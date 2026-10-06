"""双段结果冲突审计 HTTP 接口测试。

覆盖：原始/裁决两种视图、片段先合并再求交、同一源行不互比、排序、
空冲突（零汇总）、审计结果固定于发起时裁决版本、TSV 下载、
单段/过期/非法参数报错，以及大输入不做全量两两比较的扫描线行为。
"""

import os

import pytest
from fastapi.testclient import TestClient

from app import main as main_module
from app.main import app

client = TestClient(app)

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples")


def _post2(bed, chain1, chain2, names=("b.bed", "c1.chain", "c2.chain")):
    return client.post(
        "/api/liftover2",
        files={
            "bed": (names[0], bed, "text/plain"),
            "chain1": (names[1], chain1, "text/plain"),
            "chain2": (names[2], chain2, "text/plain"),
        },
        data={"side1": "target", "side2": "target"},
    )


def _post2_examples():
    with open(os.path.join(EXAMPLES, "example_double.bed"), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo_seg1.chain"), "rb") as c1, \
         open(os.path.join(EXAMPLES, "demo_seg2.chain"), "rb") as c2:
        return _post2(bf.read(), c1.read(), c2.read(),
                      ("example_double.bed", "demo_seg1.chain", "demo_seg2.chain"))


def _audit(token, view):
    return client.post(f"/api/result/{token}/audit", json={"view": view})


def _tsv_lines(resp):
    return resp.content.decode("utf-8-sig").splitlines()


# ---------- 自定义输入 ----------
#
# 视图夹具（2 行）：
#   行1 chr1:[0,100)  → 歧义：路径1 = chain1#1×chain2#10 → chrZ:[0,100)
#                        路径2 = chain1#2×chain2#11 → chrW:[0,100)
#   行2 chr1:[100,200) → chain1#3×chain2#12 → chrW:[50,150)
_VIEW_BED = b"chr1\t0\t100\nchr1\t100\t200\n"
_VIEW_CHAIN1 = (
    b"chain 1000 chr1 2000 + 0 100 chrA 500 + 0 100 1\n100\n"
    b"chain 2000 chr1 2000 + 0 100 chrB 500 + 0 100 2\n100\n"
    b"chain 3000 chr1 2000 + 100 200 chrC 500 + 0 100 3\n100\n"
)
_VIEW_CHAIN2 = (
    b"chain 500 chrA 500 + 0 100 chrZ 1000 + 0 100 10\n100\n"
    b"chain 600 chrB 500 + 0 100 chrW 1000 + 0 100 11\n100\n"
    b"chain 700 chrC 500 + 0 100 chrW 1000 + 50 150 12\n100\n"
)

# 合并夹具：行1 的路径产生两个相接的目标片段 chrZ:[0,50)+[50,100)，
# 合并后为 [0,100)；行2 → chrZ:[40,140)。交集应为单段 [40,100) 60bp。
_MERGE_BED = b"chr1\t0\t100\nchr1\t100\t200\n"
_MERGE_CHAIN1 = (
    b"chain 1000 chr1 2000 + 0 100 chrA 500 + 0 110 1\n50\t0\t10\n50\n"
    b"chain 2000 chr1 2000 + 100 200 chrB 500 + 0 100 2\n100\n"
)
_MERGE_CHAIN2 = (
    b"chain 500 chrA 500 + 0 110 chrZ 1000 + 0 100 10\n50\t10\t0\n50\n"
    b"chain 600 chrB 500 + 0 100 chrZ 1000 + 40 140 11\n100\n"
)

# 同一源行夹具：唯一的歧义行两条路径都落在 chrZ 且互相重叠——
# 同一源行的路径之间不比较，应零冲突。
_SAMELINE_BED = b"chr1\t0\t100\n"
_SAMELINE_CHAIN1 = (
    b"chain 1000 chr1 2000 + 0 100 chrA 500 + 0 100 1\n100\n"
    b"chain 2000 chr1 2000 + 0 100 chrB 500 + 0 100 2\n100\n"
)
_SAMELINE_CHAIN2 = (
    b"chain 500 chrA 500 + 0 100 chrZ 1000 + 0 100 10\n100\n"
    b"chain 600 chrB 500 + 0 100 chrZ 1000 + 10 110 11\n100\n"
)

# 多段交集夹具：行1 → chrZ:[0,10)+[100,110)，行2 → chrZ:[5,105)，
# 同一对组合产生两段交集 [5,10) 与 [100,105)。
_MULTI_BED = b"chr1\t0\t200\nchr1\t200\t300\n"
_MULTI_CHAIN1 = (
    b"chain 1000 chr1 2000 + 0 200 chrA 500 + 0 200 1\n200\n"
    b"chain 2000 chr1 2000 + 200 300 chrB 500 + 0 100 2\n100\n"
)
_MULTI_CHAIN2 = (
    b"chain 500 chrA 500 + 0 110 chrZ 1000 + 0 110 10\n10\t90\t90\n10\n"
    b"chain 600 chrB 500 + 0 100 chrZ 1000 + 5 105 11\n100\n"
)

# 排序夹具（3 行）：chrW 与 chrZ 各一段冲突，验证按染色体稳定排序。
_SORT_BED = b"chr1\t0\t100\nchr1\t100\t200\nchr1\t200\t300\n"
_SORT_CHAIN1 = _VIEW_CHAIN1 + (
    b"chain 4000 chr1 2000 + 200 300 chrD 500 + 0 100 4\n100\n"
)
_SORT_CHAIN2 = _VIEW_CHAIN2 + (
    b"chain 800 chrD 500 + 0 100 chrZ 1000 + 20 120 13\n100\n"
)


@pytest.fixture()
def view_job():
    resp = _post2(_VIEW_BED, _VIEW_CHAIN1, _VIEW_CHAIN2)
    assert resp.status_code == 200, resp.text
    rows = resp.json()["preview"]["rows"]
    assert rows[0]["status"] == "ambiguous" and rows[1]["status"] == "mapped"
    return resp.json()["download_token"]


# ---------- 示例文件：原始视图 ----------

def test_audit_examples_raw_view():
    token = _post2_examples().json()["download_token"]
    resp = _audit(token, "raw")
    assert resp.status_code == 200, resp.text
    data = resp.json()

    assert data["download_token"] == token
    assert data["view"] == "raw"
    assert data["version"] == 0
    assert data["audit_token"]
    assert data["download_url"] == f"/api/download/audit/{data['audit_token']}"
    assert ".audit.raw.v0.tsv" in data["download_name"]

    # 行3/4/5/6 有路径（4 行、5 条路径）；唯一冲突是行3×行4 在 chrZ 上
    assert data["summary"] == {
        "source_lines": 4,
        "paths": 5,
        "conflict_pairs": 1,
        "conflict_intervals": 1,
        "overlap_bases": 100,
    }
    (c,) = data["conflicts"]
    assert c == {
        "chrom": "chrZ",
        "overlap_start": 0,
        "overlap_end": 100,
        "overlap_bases": 100,
        "line_a": 3, "path_a": 1, "chain1_a": "1", "chain2_a": "10",
        "line_b": 4, "path_b": 1, "chain1_b": "1", "chain2_b": "10",
    }


# ---------- 视图语义 ----------

def test_audit_raw_view_includes_all_paths(view_job):
    data = _audit(view_job, "raw").json()
    # 行1路径2(chrW:[0,100)) × 行2(chrW:[50,150)) → [50,100) 50bp
    assert data["summary"] == {
        "source_lines": 2,
        "paths": 3,
        "conflict_pairs": 1,
        "conflict_intervals": 1,
        "overlap_bases": 50,
    }
    (c,) = data["conflicts"]
    assert (c["chrom"], c["overlap_start"], c["overlap_end"]) == ("chrW", 50, 100)
    assert (c["line_a"], c["path_a"], c["chain1_a"], c["chain2_a"]) == (1, 2, "2", "11")
    assert (c["line_b"], c["path_b"], c["chain1_b"], c["chain2_b"]) == (2, 1, "3", "12")


def test_audit_adjudicated_view_undecided_same_as_raw(view_job):
    raw = _audit(view_job, "raw").json()
    adj = _audit(view_job, "adjudicated").json()
    assert adj["version"] == 0
    assert adj["conflicts"] == raw["conflicts"]
    assert adj["summary"] == raw["summary"]


def test_audit_adjudicated_view_narrows_decided_row(view_job):
    # 裁决行1 → 路径1（chrZ）：chrW 上的冲突消失
    decide = client.post(
        f"/api/result/{view_job}/decisions",
        json={"version": 0, "source_line": 1, "candidate_index": 1},
    )
    assert decide.status_code == 200

    data = _audit(view_job, "adjudicated").json()
    assert data["version"] == 1
    assert data["conflicts"] == []
    assert data["summary"]["conflict_pairs"] == 0
    assert data["summary"]["conflict_intervals"] == 0
    assert data["summary"]["overlap_bases"] == 0
    # 参与规模：行1 只剩被选路径，行2 不受影响
    assert data["summary"]["paths"] == 2
    assert data["summary"]["source_lines"] == 2

    # 改选路径 2（chrW）：冲突恢复
    client.post(
        f"/api/result/{view_job}/decisions",
        json={"version": 1, "source_line": 1, "candidate_index": 2},
    )
    data = _audit(view_job, "adjudicated").json()
    assert data["version"] == 2
    assert len(data["conflicts"]) == 1
    assert data["conflicts"][0]["chrom"] == "chrW"


# ---------- 合并与求交语义 ----------

def test_audit_merges_own_fragments_before_intersection():
    token = _post2(_MERGE_BED, _MERGE_CHAIN1, _MERGE_CHAIN2).json()["download_token"]
    data = _audit(token, "raw").json()
    # 行1 目标片段 [0,50)+[50,100) 先合并为 [0,100)，与行2 [40,140)
    # 的交集是单段 [40,100) 60bp，而不是两段相接的小区间
    assert data["summary"]["conflict_intervals"] == 1
    (c,) = data["conflicts"]
    assert (c["overlap_start"], c["overlap_end"], c["overlap_bases"]) == (40, 100, 60)


def test_audit_same_source_line_paths_not_compared():
    token = _post2(_SAMELINE_BED, _SAMELINE_CHAIN1, _SAMELINE_CHAIN2).json()[
        "download_token"
    ]
    data = _audit(token, "raw").json()
    # 两条路径在 chrZ 上互相重叠，但属于同一源行：零冲突、零汇总
    assert data["conflicts"] == []
    assert data["summary"] == {
        "source_lines": 1,
        "paths": 2,
        "conflict_pairs": 0,
        "conflict_intervals": 0,
        "overlap_bases": 0,
    }


def test_audit_multi_interval_pair_intersection():
    token = _post2(_MULTI_BED, _MULTI_CHAIN1, _MULTI_CHAIN2).json()["download_token"]
    data = _audit(token, "raw").json()
    assert data["summary"]["conflict_pairs"] == 1
    assert data["summary"]["conflict_intervals"] == 2
    assert data["summary"]["overlap_bases"] == 10
    got = [(c["overlap_start"], c["overlap_end"], c["overlap_bases"])
           for c in data["conflicts"]]
    assert got == [(5, 10, 5), (100, 105, 5)]
    # 两段交集属于同一对组合
    assert {(c["line_a"], c["path_a"], c["line_b"], c["path_b"])
            for c in data["conflicts"]} == {(1, 1, 2, 1)}


def test_audit_results_stably_sorted():
    token = _post2(_SORT_BED, _SORT_CHAIN1, _SORT_CHAIN2).json()["download_token"]
    data = _audit(token, "raw").json()
    keys = [
        (c["chrom"], c["overlap_start"], c["overlap_end"],
         c["line_a"], c["path_a"], c["line_b"], c["path_b"])
        for c in data["conflicts"]
    ]
    assert keys == sorted(keys)
    assert keys == [
        ("chrW", 50, 100, 1, 2, 2, 1),
        ("chrZ", 20, 100, 1, 1, 3, 1),
    ]
    assert data["summary"]["overlap_bases"] == 50 + 80


# ---------- 审计结果固定于发起时的裁决版本 ----------

def test_audit_pinned_to_initiation_version(view_job):
    # v0（未裁决）发起裁决视图审计：保留全部路径，有 chrW 冲突
    first = _audit(view_job, "adjudicated").json()
    assert first["version"] == 0 and len(first["conflicts"]) == 1
    pinned_tsv = client.get(first["download_url"]).content

    # 之后裁决行1 → 路径1（chrZ），版本递增
    client.post(
        f"/api/result/{view_job}/decisions",
        json={"version": 0, "source_line": 1, "candidate_index": 1},
    )
    # 旧审计结果不变：同一 audit_token 下载字节一致
    assert client.get(first["download_url"]).content == pinned_tsv
    # 新发起的审计反映 v1 的裁决视图：冲突消失
    second = _audit(view_job, "adjudicated").json()
    assert second["version"] == 1 and second["conflicts"] == []
    # 原始视图始终包含全部路径，与裁决无关
    third = _audit(view_job, "raw").json()
    assert len(third["conflicts"]) == 1


def test_audit_does_not_modify_original_job(view_job):
    tsv_before = client.get(f"/api/download/{view_job}").content
    adj_before = client.get(f"/api/download/{view_job}/adjudicated").content
    jobs_before = len(main_module._jobs)

    resp = _audit(view_job, "raw")
    assert resp.status_code == 200

    assert len(main_module._jobs) == jobs_before
    assert client.get(f"/api/download/{view_job}").content == tsv_before
    assert client.get(f"/api/download/{view_job}/adjudicated").content == adj_before
    assert client.get(f"/api/result/{view_job}/decisions").json()["version"] == 0


# ---------- TSV 下载 ----------

def test_audit_tsv_download():
    token = _post2_examples().json()["download_token"]
    data = _audit(token, "raw").json()
    dl = client.get(data["download_url"])
    assert dl.status_code == 200
    assert dl.headers["content-type"].startswith("text/tab-separated-values")
    assert data["download_name"] in dl.headers["content-disposition"]
    assert dl.content.startswith(b"\xef\xbb\xbf")

    lines = _tsv_lines(dl)
    header = lines[0].split("\t")
    assert header == [
        "chrom", "overlap_start", "overlap_end", "overlap_bases",
        "line_a", "path_a", "chain1_a", "chain2_a",
        "line_b", "path_b", "chain1_b", "chain2_b",
    ]
    assert len(lines) == 2
    assert lines[1].split("\t") == [
        "chrZ", "0", "100", "100", "3", "1", "1", "10", "4", "1", "1", "10",
    ]


def test_audit_tsv_empty_when_no_conflict():
    token = _post2(_SAMELINE_BED, _SAMELINE_CHAIN1, _SAMELINE_CHAIN2).json()[
        "download_token"
    ]
    data = _audit(token, "raw").json()
    assert data["conflicts"] == []
    dl = client.get(data["download_url"])
    assert dl.status_code == 200
    # 无冲突：只有表头一行
    assert len(_tsv_lines(dl)) == 1


def test_audit_download_unknown_token_404():
    resp = client.get("/api/download/audit/no-such-token")
    assert resp.status_code == 404
    assert "审计结果" in resp.json()["errors"][0]


# ---------- 错误场景 ----------

def test_audit_unknown_token_404():
    resp = _audit("no-such-token", "raw")
    assert resp.status_code == 404
    assert "重新提交转换" in resp.json()["errors"][0]


def test_audit_single_segment_rejected():
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
    audits_before = len(main_module._audit_jobs)
    resp = _audit(token, "raw")
    assert resp.status_code == 400
    assert "双段" in resp.json()["errors"][0]
    assert len(main_module._audit_jobs) == audits_before


@pytest.mark.parametrize("body", [
    {},                              # 缺 view
    {"view": "bogus"},               # 非法视图
    {"view": ""},
    {"view": 1},                     # 非字符串
    {"view": True},
    {"view": ["raw"]},
    {"view": None},
])
def test_audit_bad_view_rejected(view_job, body):
    audits_before = len(main_module._audit_jobs)
    resp = client.post(f"/api/result/{view_job}/audit", json=body)
    assert resp.status_code == 400
    assert "view" in resp.json()["errors"][0]
    assert len(main_module._audit_jobs) == audits_before


def test_audit_body_not_json(view_job):
    resp = client.post(
        f"/api/result/{view_job}/audit",
        content=b"not-json",
        headers={"Content-Type": "application/json"},
    )
    assert resp.status_code == 400


# ---------- 大输入：扫描线而非全量两两比较 ----------

def test_audit_large_input_many_disjoint_rows():
    # 3000 个互不重叠的源行：扫描线一次通过即得零冲突，
    # 无需 3000×3000 的全量两两比较
    n = 3000
    bed = b"".join(f"chr1\t{i * 1000}\t{i * 1000 + 100}\n".encode()
                   for i in range(n))
    chain1 = b"chain 1000 chr1 10000000 + 0 3001000 chrA 10000000 + 0 3001000 1\n3001000\n"
    chain2 = b"chain 500 chrA 10000000 + 0 3001000 chrZ 10000000 + 0 3001000 10\n3001000\n"
    token = _post2(bed, chain1, chain2).json()["download_token"]

    data = _audit(token, "raw").json()
    assert data["summary"]["source_lines"] == n
    assert data["summary"]["paths"] == n
    assert data["summary"]["conflict_pairs"] == 0
    assert data["conflicts"] == []


def test_audit_many_overlapping_rows_all_pairs_reported():
    # 150 个源行落在同一目标区间：冲突数是组合数 C(150,2)，
    # 每对只输出一次且 line_a < line_b
    n = 150
    bed = b"chr1\t0\t100\n" * n
    chain1 = b"chain 1000 chr1 2000 + 0 100 chrA 500 + 0 100 1\n100\n"
    chain2 = b"chain 500 chrA 500 + 0 100 chrZ 1000 + 0 100 10\n100\n"
    token = _post2(bed, chain1, chain2).json()["download_token"]

    data = _audit(token, "raw").json()
    expected_pairs = n * (n - 1) // 2
    assert data["summary"]["conflict_pairs"] == expected_pairs
    assert data["summary"]["conflict_intervals"] == expected_pairs
    assert data["summary"]["overlap_bases"] == expected_pairs * 100
    first, last = data["conflicts"][0], data["conflicts"][-1]
    assert (first["line_a"], first["line_b"]) == (1, 2)
    assert (last["line_a"], last["line_b"]) == (n - 1, n)
    assert all(c["line_a"] < c["line_b"] for c in data["conflicts"])
