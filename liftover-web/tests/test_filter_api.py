"""双段结果质量筛选 HTTP 接口测试。

覆盖：筛选预览与统计、歧义行按路径分别判断、unmapped 只统计不写出、
筛选版 TSV 下载、阈值边界（含边界值通过）、参数错误（越界/类型/格式）
与不改变原任务、单段结果拒绝筛选、token 缺失等。
"""

import os

import pytest
from fastapi.testclient import TestClient

from app import main as main_module
from app.main import app

client = TestClient(app)

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples")


def _post2(
    bed_name: str = "example_double.bed",
    chain1_name: str = "demo_seg1.chain",
    chain2_name: str = "demo_seg2.chain",
):
    with open(os.path.join(EXAMPLES, bed_name), "rb") as bf, \
         open(os.path.join(EXAMPLES, chain1_name), "rb") as c1, \
         open(os.path.join(EXAMPLES, chain2_name), "rb") as c2:
        return client.post(
            "/api/liftover2",
            files={
                "bed": (bed_name, bf, "text/plain"),
                "chain1": (chain1_name, c1, "text/plain"),
                "chain2": (chain2_name, c2, "text/plain"),
            },
            data={"side1": "target", "side2": "target"},
        )


@pytest.fixture()
def job2():
    resp = _post2()
    assert resp.status_code == 200, resp.text
    return resp.json()["download_token"]


def _filter(token, min_coverage, keep_statuses=("mapped", "partial", "ambiguous")):
    return client.post(
        f"/api/result/{token}/filter",
        json={"min_coverage": min_coverage, "keep_statuses": list(keep_statuses)},
    )


def _tsv_lines(resp):
    return resp.content.decode("utf-8-sig").splitlines()


# 自定义双段输入：歧义行的两条路径最终覆盖率不同（1.0 与 0.5）
_AMB_BED = b"chr1\t0\t100\n"
_AMB_CHAIN1 = (
    b"chain 1000 chr1 2000 + 0 100 chrA 500 + 0 100 1\n100\n"
    b"chain 2000 chr1 2000 + 0 100 chrB 500 + 0 100 2\n100\n"
)
_AMB_CHAIN2 = (
    b"chain 500 chrA 500 + 0 100 chrZ 1000 + 0 100 10\n100\n"
    b"chain 600 chrB 500 + 0 50 chrW 1000 + 0 50 11\n50\n"
)


@pytest.fixture()
def amb_job():
    resp = client.post(
        "/api/liftover2",
        files={
            "bed": ("amb.bed", _AMB_BED, "text/plain"),
            "chain1": ("amb1.chain", _AMB_CHAIN1, "text/plain"),
            "chain2": ("amb2.chain", _AMB_CHAIN2, "text/plain"),
        },
        data={"side1": "target", "side2": "target"},
    )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    rows = data["preview"]["rows"]
    assert len(rows) == 1 and rows[0]["status"] == "ambiguous"
    ratios = sorted(p["seg2_covered"] / p["total"] for p in rows[0]["candidates"])
    assert ratios == [0.5, 1.0]
    return data["download_token"]


# ---------- 筛选预览与统计 ----------

def test_filter_counts_and_preview(job2):
    resp = _filter(job2, 0.9)
    assert resp.status_code == 200, resp.text
    data = resp.json()

    assert data["download_token"] == job2
    assert data["min_coverage"] == 0.9
    assert data["keep_statuses"] == ["mapped", "partial", "ambiguous"]
    assert data["filter_token"]
    assert data["download_url"] == f"/api/download/filtered/{data['filter_token']}"
    assert ".filtered.min-cov-0.9.tsv" in data["download_name"]

    # mapped(1.0) 与歧义行(1.0/1.0) 通过；partial 0.833/0.0625 被滤掉；
    # unmapped 4 行始终计入统计
    assert data["counts"] == {"mapped": 1, "partial": 0, "ambiguous": 1, "unmapped": 4}
    preview = data["preview"]
    assert preview["total"] == 2
    assert preview["counts"] == data["counts"]
    # unmapped 行不写入筛选结果
    assert all(r["status"] != "unmapped" for r in preview["rows"])
    assert {r["line_no"] for r in preview["rows"]} == {3, 6}


def test_filter_keep_statuses_subset(job2):
    # 只保留 partial：mapped/ambiguous 即使覆盖率达标也不保留
    resp = _filter(job2, 0.0, keep_statuses=("partial",))
    assert resp.status_code == 200
    data = resp.json()
    assert data["counts"] == {"mapped": 0, "partial": 2, "ambiguous": 0, "unmapped": 4}
    assert {r["line_no"] for r in data["preview"]["rows"]} == {4, 5}


def test_filter_unmapped_only_in_stats_even_when_selected(job2):
    resp = _filter(job2, 0.0, keep_statuses=("unmapped",))
    assert resp.status_code == 200
    data = resp.json()
    assert data["counts"] == {"mapped": 0, "partial": 0, "ambiguous": 0, "unmapped": 4}
    assert data["preview"]["total"] == 0
    assert data["preview"]["rows"] == []
    # 筛选版 TSV 只有表头
    dl = client.get(data["download_url"])
    assert dl.status_code == 200
    assert len(_tsv_lines(dl)) == 1


def test_filter_threshold_boundaries(job2):
    # 阈值 0：全部有路径的行通过；阈值 1：只有完全覆盖的路径通过
    all_kept = _filter(job2, 0).json()
    assert all_kept["counts"] == {"mapped": 1, "partial": 2, "ambiguous": 1, "unmapped": 4}
    assert all_kept["preview"]["total"] == 4

    full_only = _filter(job2, 1).json()
    assert full_only["counts"] == {"mapped": 1, "partial": 0, "ambiguous": 1, "unmapped": 4}
    assert {r["line_no"] for r in full_only["preview"]["rows"]} == {3, 6}


def test_filter_threshold_inclusive(job2):
    # 行4 最终覆盖率 = 100/120；阈值恰好等于该值（含边界）时通过
    exact = 100 / 120
    kept = _filter(job2, exact, keep_statuses=("partial",)).json()
    assert kept["counts"]["partial"] == 1
    assert kept["preview"]["rows"][0]["line_no"] == 4

    dropped = _filter(job2, 0.84, keep_statuses=("partial",)).json()
    assert dropped["counts"]["partial"] == 0
    assert dropped["preview"]["total"] == 0


# ---------- 歧义行按路径分别判断 ----------

def test_filter_ambiguous_per_path_keeps_passing_only(amb_job):
    resp = _filter(amb_job, 0.8)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["counts"]["ambiguous"] == 1
    row = data["preview"]["rows"][0]
    # 状态保持 ambiguous，不自动替用户裁决；只保留通过阈值的路径
    assert row["status"] == "ambiguous"
    assert len(row["candidates"]) == 1
    path = row["candidates"][0]
    assert path["path_index"] == 1          # 路径序号保持原值
    assert path["dest_chrom"] == "chrZ"
    assert path["seg2_covered"] == 100 and path["total"] == 100


def test_filter_ambiguous_keeps_all_passing_paths(amb_job):
    # 两条路径都通过时全部保留（不自动挑选）；阈值 0.5 同时验证边界含等号
    both = _filter(amb_job, 0.5).json()
    row = both["preview"]["rows"][0]
    assert row["status"] == "ambiguous"
    assert len(row["candidates"]) == 2
    assert [p["path_index"] for p in row["candidates"]] == [1, 2]

    # 阈值 1.0 时只剩完全覆盖的 chrZ 路径
    strict = _filter(amb_job, 1.0).json()
    row = strict["preview"]["rows"][0]
    assert row["status"] == "ambiguous"
    assert len(row["candidates"]) == 1
    assert row["candidates"][0]["dest_chrom"] == "chrZ"


def test_filter_ambiguous_no_path_passing_drops_row():
    # 两条路径最终覆盖率都只有 0.5：阈值 0.8 时整行不保留
    chain2 = (
        b"chain 500 chrA 500 + 0 50 chrZ 1000 + 0 50 10\n50\n"
        b"chain 600 chrB 500 + 0 50 chrW 1000 + 0 50 11\n50\n"
    )
    resp = client.post(
        "/api/liftover2",
        files={
            "bed": ("amb.bed", _AMB_BED, "text/plain"),
            "chain1": ("amb1.chain", _AMB_CHAIN1, "text/plain"),
            "chain2": ("amb2.chain", chain2, "text/plain"),
        },
        data={"side1": "target", "side2": "target"},
    )
    assert resp.status_code == 200, resp.text
    token = resp.json()["download_token"]
    assert resp.json()["preview"]["rows"][0]["status"] == "ambiguous"

    data = _filter(token, 0.8).json()
    assert data["counts"]["ambiguous"] == 0
    assert data["preview"]["total"] == 0
    # 阈值降到 0.5（含边界）则两条路径都保留
    data = _filter(token, 0.5).json()
    assert len(data["preview"]["rows"][0]["candidates"]) == 2


# ---------- 筛选版 TSV 下载 ----------

def test_filtered_tsv_download(job2):
    data = _filter(job2, 0.9).json()
    dl = client.get(data["download_url"])
    assert dl.status_code == 200
    assert dl.headers["content-type"].startswith("text/tab-separated-values")
    assert data["download_name"] in dl.headers["content-disposition"]
    assert dl.content.startswith(b"\xef\xbb\xbf")

    lines = _tsv_lines(dl)
    header = lines[0].split("\t")
    assert len(header) == 26
    assert header[:3] == ["source_line", "status", "failed_stage"]
    assert header[-1] == "final_coverage"

    body = [ln.split("\t") for ln in lines[1:]]
    # 只有行3（mapped）与行6（ambiguous 两条路径）；无 partial/unmapped
    assert {f[0] for f in body} == {"3", "6"}
    assert all(f[1] != "unmapped" for f in body)
    r3 = [f for f in body if f[0] == "3"]
    assert len(r3) == 1 and r3[0][-1] == "1"
    r6 = [f for f in body if f[0] == "6"]
    assert len(r6) == 2
    assert {f[6] for f in r6} == {"1", "2"}        # path_index 保持原序号
    assert all(f[-1] == "1" for f in r6)


def test_filtered_tsv_final_coverage_values(amb_job):
    data = _filter(amb_job, 0.0).json()
    dl = client.get(data["download_url"])
    body = [ln.split("\t") for ln in _tsv_lines(dl)[1:]]
    assert len(body) == 2
    cov_by_path = {f[6]: f[-1] for f in body}
    assert cov_by_path == {"1": "1", "2": "0.5"}


def test_filtered_download_unknown_token_404():
    resp = client.get("/api/download/filtered/no-such-token")
    assert resp.status_code == 404
    assert "筛选结果" in resp.json()["errors"][0]


# ---------- 参数错误：明确报错且不改变原任务 ----------

def _assert_job_unchanged(token, jobs_before, filters_before, tsv_before, version_before):
    assert len(main_module._jobs) == jobs_before
    assert len(main_module._filter_jobs) == filters_before
    assert client.get(f"/api/download/{token}").content == tsv_before
    assert client.get(f"/api/result/{token}/decisions").json()["version"] == version_before


@pytest.mark.parametrize("bad_cov", [-0.1, 1.5, -1, 2])
def test_filter_threshold_out_of_range(job2, bad_cov):
    jobs_before = len(main_module._jobs)
    filters_before = len(main_module._filter_jobs)
    tsv_before = client.get(f"/api/download/{job2}").content
    resp = _filter(job2, bad_cov)
    assert resp.status_code == 400
    assert "min_coverage" in resp.json()["errors"][0]
    _assert_job_unchanged(job2, jobs_before, filters_before, tsv_before, 0)


@pytest.mark.parametrize("bad_cov", ["0.8", True, [0.8], {"v": 0.8}])
def test_filter_threshold_wrong_type(job2, bad_cov):
    filters_before = len(main_module._filter_jobs)
    resp = _filter(job2, bad_cov)
    assert resp.status_code == 400
    assert "min_coverage" in resp.json()["errors"][0]
    assert len(main_module._filter_jobs) == filters_before


@pytest.mark.parametrize("raw", [b'{"min_coverage": NaN, "keep_statuses": ["mapped"]}',
                                 b'{"min_coverage": Infinity, "keep_statuses": ["mapped"]}'])
def test_filter_threshold_non_finite(job2, raw):
    resp = client.post(
        f"/api/result/{job2}/filter",
        content=raw,
        headers={"Content-Type": "application/json"},
    )
    assert resp.status_code == 400
    assert "min_coverage" in resp.json()["errors"][0]


@pytest.mark.parametrize("body", [
    {"keep_statuses": ["mapped"]},                       # 缺 min_coverage
    {"min_coverage": 0.5},                                # 缺 keep_statuses
    {"min_coverage": 0.5, "keep_statuses": []},           # 空数组
    {"min_coverage": 0.5, "keep_statuses": "mapped"},     # 不是数组
    {"min_coverage": 0.5, "keep_statuses": ["bogus"]},    # 非法状态
    {"min_coverage": 0.5, "keep_statuses": [1, 2]},       # 非字符串元素
])
def test_filter_bad_request_body(job2, body):
    filters_before = len(main_module._filter_jobs)
    resp = client.post(f"/api/result/{job2}/filter", json=body)
    assert resp.status_code == 400
    assert resp.json()["errors"]
    assert len(main_module._filter_jobs) == filters_before


def test_filter_body_not_json(job2):
    resp = client.post(
        f"/api/result/{job2}/filter",
        content=b"not-json",
        headers={"Content-Type": "application/json"},
    )
    assert resp.status_code == 400


def test_filter_unknown_token_404():
    resp = _filter("no-such-token", 0.5)
    assert resp.status_code == 404
    assert "重新提交转换" in resp.json()["errors"][0]


def test_filter_single_segment_rejected():
    with open(os.path.join(EXAMPLES, "example_target.bed"), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo.chain"), "rb") as cf:
        resp = client.post(
            "/api/liftover",
            files={"bed": ("b.bed", bf, "text/plain"),
                   "chain": ("c.chain", cf, "text/plain")},
            data={"side": "target"},
        )
    assert resp.status_code == 200
    token = resp.json()["download_token"]
    filters_before = len(main_module._filter_jobs)
    resp = _filter(token, 0.5)
    assert resp.status_code == 400
    assert "双段" in resp.json()["errors"][0]
    assert len(main_module._filter_jobs) == filters_before


# ---------- 筛选不改变原任务 ----------

def test_filter_does_not_modify_original_job(job2):
    tsv_before = client.get(f"/api/download/{job2}").content
    adj_before = client.get(f"/api/download/{job2}/adjudicated").content

    resp = _filter(job2, 0.95)
    assert resp.status_code == 200

    # 原始下载与裁决版字节完全不变，版本号不变
    assert client.get(f"/api/download/{job2}").content == tsv_before
    assert client.get(f"/api/download/{job2}/adjudicated").content == adj_before
    assert client.get(f"/api/result/{job2}/decisions").json()["version"] == 0


def test_filter_ignores_adjudication_state(job2):
    # 裁决（选定行6 路径2）不影响筛选：筛选始终基于原始结果的全部路径
    decide = client.post(
        f"/api/result/{job2}/decisions",
        json={"version": 0, "source_line": 6, "candidate_index": 2},
    )
    assert decide.status_code == 200

    data = _filter(job2, 0.9).json()
    row6 = next(r for r in data["preview"]["rows"] if r["line_no"] == 6)
    assert row6["status"] == "ambiguous"
    assert len(row6["candidates"]) == 2

    dl = client.get(data["download_url"])
    r6 = [ln.split("\t") for ln in _tsv_lines(dl)[1:] if ln.startswith("6\t")]
    assert len(r6) == 2
    # 裁决状态本身也不受筛选影响
    got = client.get(f"/api/result/{job2}/decisions").json()
    assert got["version"] == 1
    assert got["decisions"][0]["candidate_index"] == 2
