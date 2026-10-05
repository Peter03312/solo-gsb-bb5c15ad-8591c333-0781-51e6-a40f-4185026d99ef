"""FastAPI 入口：上传、预览、下载、歧义人工裁决、双段结果质量筛选。

结果只在校验完全通过后保存在进程内存中（BytesIO 等价物），
任何解析错误都直接 400 返回，不创建、更不提供任何可下载文件。

- ``POST /api/liftover``：单段转换（一份 chain）；
- ``POST /api/liftover2``：双段转换（同一 BED 依次经过两份 chain，
  第二段只处理第一段实际映射到的片段，最终片段回溯到原 BED 源坐标）。

歧义裁决只作用于某一次转换结果（token 标识），绝不回写原始转换：

- 原始结果 ``GET /api/download/{token}`` 的字节内容在裁决前后完全一致；
- 裁决状态（选了哪个候选）带版本号做乐观并发控制，相同版本的并发
  修改只有一次能成功；
- 裁决版 TSV 由 ``GET /api/download/{token}/adjudicated`` 单独导出。

双段结果的质量筛选（``POST /api/result/{token}/filter``）同样只是
原结果的导出视图：按每条路径的最终覆盖率（seg2_covered/total）与保留
状态筛选，产物挂在独立的 filter token 下
（``GET /api/download/filtered/{filter_token}``），原任务的字节内容、
裁决状态与版本号在筛选前后完全不变。
"""

from __future__ import annotations

import asyncio
import math
import os
import secrets
import time
from dataclasses import dataclass, field

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles

from .bed import parse_bed
from .chain import parse_chain
from .engine import RowResult
from .engine import LiftOverEngine
from .engine2 import TwoSegmentEngine
from .filtering import FILTER_STATUSES, filter_double_rows
from .report import (
    build_adjudicated_tsv,
    build_adjudicated_tsv2,
    build_filtered_tsv2,
    build_tsv,
    build_tsv2,
    rows_to_preview,
)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
STATIC_DIR = os.path.join(BASE_DIR, "..", "static")

MAX_UPLOAD_BYTES = int(os.getenv("LIFTOVER_MAX_UPLOAD_MB", "50")) * 1024 * 1024
JOB_TTL_SECONDS = int(os.getenv("LIFTOVER_JOB_TTL_SECONDS", "3600"))
MAX_JOBS = int(os.getenv("LIFTOVER_MAX_JOBS", "100"))

app = FastAPI(title="离线基因组区间迁移 (BED liftover)", version="1.3.0")


@dataclass
class _Job:
    rows: list                   # 单段 RowResult 列表或双段 RowResult2 列表
    tsv_bytes: bytes             # 原始结果（永不改变）
    filename: str
    created_at: float
    kind: str = "single"         # single / double
    # 裁决状态：{源物理行号: 1 基候选序号}；撤销即删除键
    decisions: dict[int, int] = field(default_factory=dict)
    version: int = 0               # 每次成功修改递增
    lock: asyncio.Lock = field(default_factory=asyncio.Lock)


@dataclass
class _FilterJob:
    """一次质量筛选的产物（独立于原任务，不改动原任务任何状态）。"""

    tsv_bytes: bytes             # 筛选版 TSV（生成后不变）
    filename: str
    created_at: float
    parent_token: str            # 来源任务 token（仅溯源用）


_jobs: dict[str, _Job] = {}
_filter_jobs: dict[str, _FilterJob] = {}


def _prune_jobs() -> None:
    now = time.time()
    expired = [k for k, j in _jobs.items() if now - j.created_at > JOB_TTL_SECONDS]
    for k in expired:
        _jobs.pop(k, None)
    # 兜底：超过容量时淘汰最早的任务
    if len(_jobs) > MAX_JOBS:
        for k in sorted(_jobs, key=lambda x: _jobs[x].created_at)[: len(_jobs) - MAX_JOBS]:
            _jobs.pop(k, None)
    expired_f = [k for k, j in _filter_jobs.items()
                 if now - j.created_at > JOB_TTL_SECONDS]
    for k in expired_f:
        _filter_jobs.pop(k, None)
    if len(_filter_jobs) > MAX_JOBS:
        for k in sorted(_filter_jobs,
                        key=lambda x: _filter_jobs[x].created_at)[: len(_filter_jobs) - MAX_JOBS]:
            _filter_jobs.pop(k, None)


def _read_upload(f: UploadFile) -> str:
    data = f.file.read()
    if not data:
        raise HTTPException(status_code=400, detail={"errors": [f"文件 {f.filename} 为空"]})
    if len(data) > MAX_UPLOAD_BYTES:
        raise HTTPException(
            status_code=413,
            detail={"errors": [f"文件 {f.filename} 超过大小上限 "
                               f"{MAX_UPLOAD_BYTES // (1024 * 1024)} MB"]},
        )
    try:
        # BOM 仅允许出现在文件最前面
        return data.decode("utf-8-sig")
    except UnicodeDecodeError:
        try:
            return data.decode("latin-1")
        except UnicodeDecodeError:
            raise HTTPException(
                status_code=400,
                detail={"errors": [f"文件 {f.filename} 不是有效的 UTF-8 文本"]},
            )


@app.exception_handler(HTTPException)
async def _http_exc_handler(request: Request, exc: HTTPException) -> JSONResponse:
    return JSONResponse(status_code=exc.status_code, content=exc.detail
                        if isinstance(exc.detail, dict) else {"errors": [str(exc.detail)]})


@app.get("/healthz")
async def healthz() -> dict[str, str]:
    return {"status": "ok"}


@app.post("/api/liftover")
async def liftover(
    bed: UploadFile = File(...),
    chain: UploadFile = File(...),
    side: str = Form(...),
):
    if side not in ("target", "query"):
        raise HTTPException(status_code=400,
                            detail={"errors": ["side 必须是 target 或 query"]})

    bed_text = _read_upload(bed)
    chain_text = _read_upload(chain)

    bed_result = parse_bed(bed_text)
    chain_result = parse_chain(chain_text)

    errors = list(bed_result.errors) + list(chain_result.errors)
    if errors:
        # 校验失败：只返回错误，绝不写入 _jobs —— 不存在可下载文件
        raise HTTPException(status_code=400, detail={"errors": errors})

    engine = LiftOverEngine(chain_result.chains, side)
    rows = engine.convert_all(bed_result.records)
    preview = rows_to_preview(rows)

    tsv = build_tsv(rows)
    # Excel 友好的 UTF-8 BOM；HTTP 下载本质上仍是内存字节，不落盘
    payload = b"\xef\xbb\xbf" + tsv.encode("utf-8")
    token = secrets.token_urlsafe(24)
    stamp = time.strftime("%Y%m%d-%H%M%S", time.gmtime())
    filename = f"liftover_{side}_to_{'query' if side == 'target' else 'target'}_{stamp}.tsv"
    _prune_jobs()
    _jobs[token] = _Job(
        rows=rows, tsv_bytes=payload, filename=filename, created_at=time.time(),
    )

    return {
        "download_token": token,
        "download_name": filename,
        "kind": "single",
        "preview": preview,
    }


@app.post("/api/liftover2")
async def liftover2(
    bed: UploadFile = File(...),
    chain1: UploadFile = File(...),
    chain2: UploadFile = File(...),
    side1: str = Form(...),
    side2: str = Form(...),
):
    """双段转换：同一 BED 依次经过 chain 文件 1 与 chain 文件 2。

    第二段只处理第一段实际映射到的中间片段；最终片段回溯到原 BED 源
    坐标。任一文件（BED / chain1 / chain2）有格式错误都带行号整体拒绝，
    不生成任何可下载结果。
    """
    side_errors = [
        f"{name} 必须是 target 或 query"
        for name, val in (("side1", side1), ("side2", side2))
        if val not in ("target", "query")
    ]
    if side_errors:
        raise HTTPException(status_code=400, detail={"errors": side_errors})

    bed_text = _read_upload(bed)
    chain1_text = _read_upload(chain1)
    chain2_text = _read_upload(chain2)

    bed_result = parse_bed(bed_text)
    chain1_result = parse_chain(chain1_text)
    chain2_result = parse_chain(chain2_text)

    errors = list(bed_result.errors)
    errors += [f"chain 文件 1：{e}" for e in chain1_result.errors]
    errors += [f"chain 文件 2：{e}" for e in chain2_result.errors]
    if errors:
        # 校验失败：只返回错误，绝不写入 _jobs —— 不存在可下载文件
        raise HTTPException(status_code=400, detail={"errors": errors})

    engine = TwoSegmentEngine(
        chain1_result.chains, side1, chain2_result.chains, side2
    )
    rows = engine.convert_all(bed_result.records)
    preview = rows_to_preview(rows)

    tsv = build_tsv2(rows)
    payload = b"\xef\xbb\xbf" + tsv.encode("utf-8")
    token = secrets.token_urlsafe(24)
    stamp = time.strftime("%Y%m%d-%H%M%S", time.gmtime())
    filename = f"liftover_double_{side1}-{side2}_{stamp}.tsv"
    _prune_jobs()
    _jobs[token] = _Job(
        rows=rows, tsv_bytes=payload, filename=filename,
        created_at=time.time(), kind="double",
    )

    return {
        "download_token": token,
        "download_name": filename,
        "kind": "double",
        "preview": preview,
    }


def _get_job_or_404(token: str) -> _Job:
    _prune_jobs()
    job = _jobs.get(token)
    if job is None:
        raise HTTPException(
            status_code=404,
            detail={"errors": ["结果不存在或已过期，请重新提交转换"]},
        )
    return job


def _ambiguous_line_map(job: _Job) -> dict[int, RowResult]:
    return {r.line_no: r for r in job.rows if r.status == "ambiguous"}


@app.get("/api/download/{token}")
async def download(token: str) -> Response:
    job = _get_job_or_404(token)
    # 原始结果：裁决前后字节完全一致
    return Response(
        content=job.tsv_bytes,
        media_type="text/tab-separated-values; charset=utf-8",
        headers={"Content-Disposition": f'attachment; filename="{job.filename}"'},
    )


@app.get("/api/result/{token}/decisions")
async def list_decisions(token: str) -> dict:
    """查看某次结果当前的歧义裁决选择与结果版本。"""
    job = _get_job_or_404(token)
    ambiguous = _ambiguous_line_map(job)
    return {
        "download_token": token,
        "version": job.version,
        "decisions": [
            {
                "source_line": line_no,
                "candidate_index": job.decisions.get(line_no),
                "candidate_count": len(ambiguous[line_no].candidates),
            }
            for line_no in sorted(ambiguous)
        ],
    }


async def _read_json_body(request: Request) -> dict:
    try:
        body = await request.json()
    except Exception:
        raise HTTPException(
            status_code=400, detail={"errors": ["请求体必须是 JSON 对象"]}
        )
    if not isinstance(body, dict):
        raise HTTPException(
            status_code=400, detail={"errors": ["请求体必须是 JSON 对象"]}
        )
    return body


def _require_int(body: dict, key: str) -> int:
    if key not in body:
        raise HTTPException(status_code=400, detail={"errors": [f"缺少字段 {key}"]})
    val = body[key]
    # bool 是 int 的子类，显式拒绝，避免 True 被当成 1
    if not isinstance(val, int) or isinstance(val, bool):
        raise HTTPException(status_code=400, detail={"errors": [f"字段 {key} 必须是整数"]})
    return val


@app.post("/api/result/{token}/decisions")
async def submit_decisions(token: str, request: Request) -> dict:
    """提交歧义行裁决（选定候选 / 撤销选择）。

    请求体：``{"version": <当前结果版本>, "source_line": <BED 物理行号>,
    "candidate_index": <1 基候选序号，或 null 表示撤销>}``。

    - 结果不存在/过期 → 404；
    - source_line 不存在、该行不是歧义行、候选序号越界 → 400；
    - version 与当前版本不一致（含相同版本的并发修改竞跑失败者）→ 409；
    - 仅当选择确实发生变化时版本才 +1；与当前一致的重复提交不递增。
    """
    body = await _read_json_body(request)
    version = _require_int(body, "version")
    source_line = _require_int(body, "source_line")
    if "candidate_index" not in body:
        raise HTTPException(status_code=400,
                            detail={"errors": ["缺少字段 candidate_index"]})
    candidate_index = body["candidate_index"]
    if candidate_index is not None and (
        not isinstance(candidate_index, int) or isinstance(candidate_index, bool)
    ):
        raise HTTPException(
            status_code=400,
            detail={"errors": ["candidate_index 必须是正整数或 null"]},
        )

    job = _get_job_or_404(token)

    async with job.lock:
        if version != job.version:
            raise HTTPException(
                status_code=409,
                detail={
                    "errors": [
                        f"结果版本已过期：提交基于版本 {version}，当前为 {job.version}，"
                        "请先重新查看当前选择后再提交"
                    ],
                    "current_version": job.version,
                },
            )

        row = next((r for r in job.rows if r.line_no == source_line), None)
        if row is None:
            raise HTTPException(
                status_code=400,
                detail={"errors": [f"源文件第 {source_line} 行不存在于本次结果中"]},
            )
        if row.status != "ambiguous":
            raise HTTPException(
                status_code=400,
                detail={"errors": [
                    f"源文件第 {source_line} 行状态为 {row.status}，"
                    "只有歧义（ambiguous）行允许人工裁决"
                ]},
            )

        n_cand = len(row.candidates)
        if candidate_index is not None and not (1 <= candidate_index <= n_cand):
            raise HTTPException(
                status_code=400,
                detail={"errors": [
                    f"候选序号越界：第 {source_line} 行共有 {n_cand} 个候选，"
                    f"收到 candidate_index={candidate_index}，允许范围 1..{n_cand}"
                ]},
            )

        previous = job.decisions.get(source_line)
        modified = previous != candidate_index
        if modified:
            if candidate_index is None:
                job.decisions.pop(source_line, None)
            else:
                job.decisions[source_line] = candidate_index
            job.version += 1

        return {
            "download_token": token,
            "version": job.version,
            "modified": modified,
            "source_line": source_line,
            "candidate_index": candidate_index,
        }


def _adjudicated_name(filename: str, version: int) -> str:
    if filename.endswith(".tsv"):
        return f"{filename[:-4]}.adjudicated.v{version}.tsv"
    return f"{filename}.adjudicated.v{version}.tsv"


@app.get("/api/download/{token}/adjudicated")
async def download_adjudicated(token: str) -> Response:
    """下载裁决版 TSV。

    已选歧义行仅输出被选候选（双段为被选路径）的全部片段；未选歧义行
    输出全部候选并标明“未裁决”；其他状态行照原结果输出。原始 status
    列保持不变。
    """
    job = _get_job_or_404(token)
    if job.kind == "double":
        text = build_adjudicated_tsv2(job.rows, job.decisions)
    else:
        text = build_adjudicated_tsv(job.rows, job.decisions)
    payload = b"\xef\xbb\xbf" + text.encode("utf-8")
    filename = _adjudicated_name(job.filename, job.version)
    return Response(
        content=payload,
        media_type="text/tab-separated-values; charset=utf-8",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


# ---------------------------------------------------------------------------
# 双段结果质量筛选
# ---------------------------------------------------------------------------

def _require_min_coverage(body: dict) -> float:
    if "min_coverage" not in body:
        raise HTTPException(status_code=400,
                            detail={"errors": ["缺少字段 min_coverage"]})
    val = body["min_coverage"]
    # bool 是 int 的子类，显式拒绝，避免 true 被当成 1
    if isinstance(val, bool) or not isinstance(val, (int, float)):
        raise HTTPException(
            status_code=400,
            detail={"errors": [f"min_coverage 必须是 0~1 之间的数值，收到 {val!r}"]},
        )
    if not math.isfinite(val):
        raise HTTPException(
            status_code=400,
            detail={"errors": [f"min_coverage 必须是有限数值，收到 {val!r}"]},
        )
    if not 0.0 <= val <= 1.0:
        raise HTTPException(
            status_code=400,
            detail={"errors": [
                f"min_coverage 超出允许范围 [0, 1]：收到 {val!r}，"
                "请传入 0~1 之间（含边界）的最小最终覆盖率"
            ]},
        )
    return float(val)


def _require_keep_statuses(body: dict) -> list[str]:
    allowed = "、".join(FILTER_STATUSES)
    if "keep_statuses" not in body:
        raise HTTPException(status_code=400,
                            detail={"errors": ["缺少字段 keep_statuses"]})
    val = body["keep_statuses"]
    if not isinstance(val, list) or not val:
        raise HTTPException(
            status_code=400,
            detail={"errors": [
                f"keep_statuses 必须是非空数组，元素取自：{allowed}"
            ]},
        )
    bad = [s for s in val if not isinstance(s, str) or s not in FILTER_STATUSES]
    if bad:
        raise HTTPException(
            status_code=400,
            detail={"errors": [
                f"keep_statuses 含非法状态 {bad!r}；允许值：{allowed}"
            ]},
        )
    # 去重并保持用户给出的顺序
    return list(dict.fromkeys(val))


def _filtered_name(filename: str, min_coverage: float) -> str:
    cov = f"{min_coverage:g}"
    if filename.endswith(".tsv"):
        return f"{filename[:-4]}.filtered.min-cov-{cov}.tsv"
    return f"{filename}.filtered.min-cov-{cov}.tsv"


@app.post("/api/result/{token}/filter")
async def filter_result(token: str, request: Request) -> dict:
    """双段结果质量筛选：最小最终覆盖率 + 保留状态。

    请求体：``{"min_coverage": <0~1 数值>, "keep_statuses": ["mapped",
    "partial", "ambiguous", "unmapped"]}``。

    - 只适用于双段转换结果（``POST /api/liftover2`` 产生的 token）；
      单段结果 → 400，token 不存在/过期 → 404；
    - 覆盖率按每条路径的 seg2_covered/total 计算，达到阈值（含边界）
      的路径通过；歧义行保留通过阈值的全部路径，不自动裁决；
    - unmapped 行始终计入返回的统计（counts），但不写入筛选结果；
    - 参数非法（阈值越界、类型错误、状态集合为空或含非法值）→ 400，
      且原任务不发生任何改变；
    - 成功时返回筛选后的预览、各状态数量与 filter token，
      筛选版 TSV 由 ``GET /api/download/filtered/{filter_token}`` 下载。
    """
    body = await _read_json_body(request)
    min_coverage = _require_min_coverage(body)
    keep_statuses = _require_keep_statuses(body)

    job = _get_job_or_404(token)
    if job.kind != "double":
        raise HTTPException(
            status_code=400,
            detail={"errors": [
                "质量筛选仅适用于双段转换结果（POST /api/liftover2）；"
                "当前结果是单段转换，请重新提交双段转换后再筛选"
            ]},
        )

    kept_rows, counts = filter_double_rows(
        job.rows, min_coverage, set(keep_statuses)
    )
    tsv = build_filtered_tsv2(kept_rows)
    payload = b"\xef\xbb\xbf" + tsv.encode("utf-8")
    filter_token = secrets.token_urlsafe(24)
    filename = _filtered_name(job.filename, min_coverage)
    _prune_jobs()
    _filter_jobs[filter_token] = _FilterJob(
        tsv_bytes=payload, filename=filename,
        created_at=time.time(), parent_token=token,
    )

    # 预览与原始结果同构（total/counts/rows）；counts 为各状态统计，
    # 其中 unmapped 始终计入，但 rows/total 只含写入筛选结果的行
    preview = rows_to_preview(kept_rows)
    preview["counts"] = counts
    return {
        "download_token": token,
        "filter_token": filter_token,
        "download_name": filename,
        "download_url": f"/api/download/filtered/{filter_token}",
        "min_coverage": min_coverage,
        "keep_statuses": keep_statuses,
        "counts": counts,
        "preview": preview,
    }


@app.get("/api/download/filtered/{filter_token}")
async def download_filtered(filter_token: str) -> Response:
    """下载质量筛选版 TSV（双段 25 列 + final_coverage 列）。"""
    _prune_jobs()
    fjob = _filter_jobs.get(filter_token)
    if fjob is None:
        raise HTTPException(
            status_code=404,
            detail={"errors": ["筛选结果不存在或已过期，请重新提交筛选"]},
        )
    return Response(
        content=fjob.tsv_bytes,
        media_type="text/tab-separated-values; charset=utf-8",
        headers={"Content-Disposition": f'attachment; filename="{fjob.filename}"'},
    )


# 静态资源（纯本地，无任何外部 CDN 依赖）
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


@app.get("/", response_class=HTMLResponse)
async def index() -> FileResponse:
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))
