"""结果导出：TSV 文本与预览 JSON。

TSV 列（英文表头，避免下载文件里混入列名歧义）::

    source_line, status, source_chrom, source_start, source_end,
    candidate_index, chain_id, dest_chrom, dest_strand,
    fragment_index, src_start, src_end, dest_start, dest_end, reason

每个源行号可能对应多行（多个片段 / 多个歧义候选）。
unmapped 行没有候选，相关列写 '.'。

裁决版 TSV（``build_adjudicated_tsv``）在末尾追加一列 ``decision``：

- 已选候选的歧义行：仅输出被选候选的全部片段，decision=已裁决；
- 未选歧义行：仍输出全部候选的全部片段，decision=未裁决；
- 其他状态行：照原结果逐行输出，decision='.'。

无论是否裁决，``status`` 列始终保留引擎给出的原始状态——
人工裁决不会把 partial 改称 mapped，也不会把 ambiguous 改成别的。

双段（两份 chain 依次迁移）结果使用 ``build_tsv2`` /
``build_adjudicated_tsv2``，列定义见文末 ``TSV2_HEADER``。
双段结果的质量筛选版使用 ``build_filtered_tsv2``：在双段 25 列末尾
追加 ``final_coverage`` 列，且只包含通过筛选的行与路径（unmapped 行
绝不写入）。
"""

from __future__ import annotations

import csv
import io
from dataclasses import asdict
from typing import Any, Mapping

from .engine import Candidate, Fragment, RowResult
from .engine2 import PathCandidate, PathFragment, RowResult2, format_uncovered

STATUS_LABELS = {
    "mapped": "mapped",
    "partial": "partial",
    "ambiguous": "ambiguous",
    "unmapped": "unmapped",
}

TSV_HEADER = [
    "source_line",
    "status",
    "source_chrom",
    "source_start",
    "source_end",
    "candidate_index",
    "chain_id",
    "dest_chrom",
    "dest_strand",
    "fragment_index",
    "src_start",
    "src_end",
    "dest_start",
    "dest_end",
    "reason",
]

# 裁决版在原列尾部追加一列（不改动前面任何列与原始 status）
DECISION_HEADER = "decision"
ADJ_TSV_HEADER = TSV_HEADER + [DECISION_HEADER]

DECISION_PICKED = "已裁决"
DECISION_PENDING = "未裁决"
DECISION_NA = "."


def _fragment_values(
    row: RowResult,
    cand: Candidate,
    ci: int,
    frag: Fragment,
    fi: int,
) -> list[Any]:
    return [
        row.line_no, row.status, row.src_chrom, row.src_start, row.src_end,
        ci, cand.chain_id, cand.dest_chrom, cand.dest_strand,
        fi, frag.src_start, frag.src_end, frag.dest_start, frag.dest_end,
        row.reason,
    ]


def build_tsv(rows: list[RowResult]) -> str:
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n",
                        quoting=csv.QUOTE_MINIMAL)
    writer.writerow(TSV_HEADER)
    for row in rows:
        if not row.candidates:
            writer.writerow([
                row.line_no, row.status, row.src_chrom, row.src_start, row.src_end,
                ".", ".", ".", ".", ".", ".", ".", ".", ".", row.reason,
            ])
            continue
        for ci, cand in enumerate(row.candidates, start=1):
            for fi, frag in enumerate(cand.fragments, start=1):
                writer.writerow(_fragment_values(row, cand, ci, frag, fi))
    return buf.getvalue()


def build_adjudicated_tsv(
    rows: list[RowResult],
    decisions: Mapping[int, int],
) -> str:
    """生成裁决版 TSV。

    ``decisions`` 为 {源物理行号: 1 基候选序号}，只可能落在歧义行上
    （接口层保证；此处对非法键做容错，按未裁决处理）。

    - 已选行：只写被选候选的全部片段，decision 列写 ``已裁决``；
      candidate_index 沿用其在原始候选列表中的 1 基序号，便于与原结果对照。
    - 未选歧义行：写全部候选、全部片段，decision 列写 ``未裁决``。
    - mapped / partial / unmapped 行：前 15 列与 ``build_tsv`` 逐行一致，
      decision 列写 '.'。

    status 列在任何情况下都保持引擎给出的原始值。
    """
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n",
                        quoting=csv.QUOTE_MINIMAL)
    writer.writerow(ADJ_TSV_HEADER)
    for row in rows:
        if not row.candidates:
            writer.writerow([
                row.line_no, row.status, row.src_chrom, row.src_start, row.src_end,
                ".", ".", ".", ".", ".", ".", ".", ".", ".", row.reason,
                DECISION_NA,
            ])
            continue

        picked = decisions.get(row.line_no)
        is_picked_ambiguous = (
            row.status == "ambiguous"
            and picked is not None
            and 1 <= picked <= len(row.candidates)
        )
        if is_picked_ambiguous:
            cand = row.candidates[picked - 1]
            for fi, frag in enumerate(cand.fragments, start=1):
                values = _fragment_values(row, cand, picked, frag, fi)
                writer.writerow(values + [DECISION_PICKED])
        else:
            decision = DECISION_PENDING if row.status == "ambiguous" else DECISION_NA
            for ci, cand in enumerate(row.candidates, start=1):
                for fi, frag in enumerate(cand.fragments, start=1):
                    values = _fragment_values(row, cand, ci, frag, fi)
                    writer.writerow(values + [decision])
    return buf.getvalue()


def row_to_dict(row: RowResult) -> dict[str, Any]:
    return asdict(row)


def rows_to_preview(rows) -> dict[str, Any]:
    """预览 JSON：单段 RowResult 与双段 RowResult2 通用（状态集合一致）。"""
    counts = {"mapped": 0, "partial": 0, "ambiguous": 0, "unmapped": 0}
    for r in rows:
        counts[r.status] += 1
    return {
        "total": len(rows),
        "counts": counts,
        "rows": [asdict(r) for r in rows],
    }


# ---------------------------------------------------------------------------
# 双段（两份 chain 依次迁移）结果导出
#
# TSV 列（英文表头）::
#
#     source_line, status, failed_stage, source_chrom, source_start, source_end,
#     path_index, chain1_id, chain2_id, mid_chrom, mid_strand,
#     dest_chrom, dest_strand, fragment_index,
#     src_start, src_end, mid_start, mid_end, dest_start, dest_end,
#     seg1_covered, seg2_covered, total_len, uncovered_src, reason
#
# - src_start/src_end 是回溯到原 BED 的源坐标；mid_* 为中间坐标
#   （第一段目标 = 第二段源）；dest_* 为最终坐标；
# - 每个最终片段占一行，路径序号 path_index 从 1 开始；
# - unmapped 行保留源物理行号与失败阶段（failed_stage = segment1/segment2），
#   路径/片段列写 '.'，uncovered_src 为整个源区间；
# - 裁决版同样在末尾追加 decision 列，规则与单段一致（针对路径序号）。
# ---------------------------------------------------------------------------

TSV2_HEADER = [
    "source_line",
    "status",
    "failed_stage",
    "source_chrom",
    "source_start",
    "source_end",
    "path_index",
    "chain1_id",
    "chain2_id",
    "mid_chrom",
    "mid_strand",
    "dest_chrom",
    "dest_strand",
    "fragment_index",
    "src_start",
    "src_end",
    "mid_start",
    "mid_end",
    "dest_start",
    "dest_end",
    "seg1_covered",
    "seg2_covered",
    "total_len",
    "uncovered_src",
    "reason",
]

ADJ_TSV2_HEADER = TSV2_HEADER + [DECISION_HEADER]


def _path_fragment_values(
    row: RowResult2,
    path: PathCandidate,
    frag: PathFragment,
    fi: int,
) -> list[Any]:
    return [
        row.line_no, row.status, row.failed_stage or ".",
        row.src_chrom, row.src_start, row.src_end,
        path.path_index, path.chain1_id, path.chain2_id,
        path.mid_chrom, path.mid_strand, path.dest_chrom, path.dest_strand,
        fi, frag.src_start, frag.src_end,
        frag.mid_start, frag.mid_end, frag.dest_start, frag.dest_end,
        path.seg1_covered, path.seg2_covered, path.total,
        format_uncovered(path.uncovered_src),
        row.reason,
    ]


def _unmapped2_values(row: RowResult2) -> list[Any]:
    return [
        row.line_no, row.status, row.failed_stage or ".",
        row.src_chrom, row.src_start, row.src_end,
        ".", ".", ".", ".", ".", ".", ".",
        ".", ".", ".", ".", ".", ".", ".",
        row.seg1_covered, row.seg2_covered, row.src_end - row.src_start,
        format_uncovered([(row.src_start, row.src_end)]),
        row.reason,
    ]


def build_tsv2(rows: list[RowResult2]) -> str:
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n",
                        quoting=csv.QUOTE_MINIMAL)
    writer.writerow(TSV2_HEADER)
    for row in rows:
        if not row.candidates:
            writer.writerow(_unmapped2_values(row))
            continue
        for path in row.candidates:
            for fi, frag in enumerate(path.fragments, start=1):
                writer.writerow(_path_fragment_values(row, path, frag, fi))
    return buf.getvalue()


def build_adjudicated_tsv2(
    rows: list[RowResult2],
    decisions: Mapping[int, int],
) -> str:
    """双段裁决版 TSV：规则与单段一致，只是“候选”换成“路径”。

    - 已选歧义行：仅输出被选路径的全部片段，decision=已裁决；
    - 未选歧义行：输出全部路径的全部片段，decision=未裁决；
    - 其他状态行：照原结果输出，decision='.'；status 列保持原值。
    """
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n",
                        quoting=csv.QUOTE_MINIMAL)
    writer.writerow(ADJ_TSV2_HEADER)
    for row in rows:
        if not row.candidates:
            writer.writerow(_unmapped2_values(row) + [DECISION_NA])
            continue

        picked = decisions.get(row.line_no)
        is_picked_ambiguous = (
            row.status == "ambiguous"
            and picked is not None
            and 1 <= picked <= len(row.candidates)
        )
        if is_picked_ambiguous:
            path = row.candidates[picked - 1]
            for fi, frag in enumerate(path.fragments, start=1):
                values = _path_fragment_values(row, path, frag, fi)
                writer.writerow(values + [DECISION_PICKED])
        else:
            decision = DECISION_PENDING if row.status == "ambiguous" else DECISION_NA
            for path in row.candidates:
                for fi, frag in enumerate(path.fragments, start=1):
                    values = _path_fragment_values(row, path, frag, fi)
                    writer.writerow(values + [decision])
    return buf.getvalue()


# ---------------------------------------------------------------------------
# 双段结果的质量筛选版导出
#
# 在双段 25 列末尾追加一列 ``final_coverage``（该路径 seg2_covered/total_len，
# 即筛选时使用的最终覆盖率），只包含通过筛选的行与路径：
#
# - unmapped 行按规则绝不写入；
# - 歧义行只写出通过阈值的路径（path_index 保持原序号），status 列保持原值。
# ---------------------------------------------------------------------------

FILTERED_TSV2_HEADER = TSV2_HEADER + ["final_coverage"]


def format_ratio(value: float) -> str:
    """覆盖率比值的紧凑文本（如 ``1``、``0.833333``、``0.0625``），最多 6 位小数。"""
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return text if text else "0"


def build_filtered_tsv2(rows: list[RowResult2]) -> str:
    """质量筛选版双段 TSV。

    ``rows`` 应已通过 ``filtering.filter_double_rows`` 筛选（不含 unmapped
    行、路径均已达到阈值）；此处对无候选行同样跳过，保证 unmapped 行绝不
    写入筛选结果。
    """
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n",
                        quoting=csv.QUOTE_MINIMAL)
    writer.writerow(FILTERED_TSV2_HEADER)
    for row in rows:
        for path in row.candidates:
            ratio = format_ratio(path.final_coverage)
            for fi, frag in enumerate(path.fragments, start=1):
                values = _path_fragment_values(row, path, frag, fi)
                writer.writerow(values + [ratio])
    return buf.getvalue()
