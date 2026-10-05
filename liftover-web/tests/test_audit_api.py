"""双段结果源行冲突审计 HTTP 接口测试。

覆盖：原始/裁决两种视图、同一源行多路径不互比、自身片段先合并再求交、
排序稳定性、无冲突的空明细零汇总、审计产物固定于发起时裁决版本、
单段/过期/非法参数报错、扫描线不做全量两两比较、审计不改变原任务。
"""

import os

import pytest
from fastapi.testclient import TestClient

from app import audit as audit_module
from app import main as main_module
from app.main import app

client = TestClient(app)

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples")


def _post2(bed, chain1, chain2):
    return client.post(
        "/api/liftover2",
        files={
            "bed": ("b.bed", bed, "text/plain"),
            "chain1": ("c1.chain", chain1, "text/plain"),
            "chain2": ("c2.chain", chain2, "text/plain"),
        },
        data={"side1": "target", "side2": "target"},
    )


def _audit(token, view=None):
    body = {} if view is None else {"view": view}
    return client.post(f"/api/result/{token}/audit", json=body)


def _decide(token, version, source_line, candidate_index):
    return client.post(
        f"/api/result/{token}/decisions",
        json={"version": version, "source_line": source_line,
              "candidate_index": candidate_index},
    )


def _tsv_lines(resp):
    return resp.content.decode("utf-8-sig").splitlines()


# ---------- 场景数据 ----------

# 恒等两段：chr1 -> chrA -> chrZ（全长单块）
_IDENT_CHAIN1 = b"chain 1000 chr1 1000000 + 0 1000000 chrA 1000000 + 0 1000000 1\n1000000\n"
_IDENT_CHAIN2 = b"chain 500 chrA 1000000 + 0 1000000 chrZ 1000000 + 0 1000000 10\n1000000\n"

# 双路径场景：chr1 同时经 chain1#1 (->chrA) 与 chain1#2 (->chrB)，
# 第二段 chrA->chrZ[0,200)、chrB->chrZ[100,300)（最终都落在 chrZ）
_TWO_PATH_CHAIN1 = (
    b"chain 1000 chr1 1000 + 0 200 chrA 1000 + 0 200 1\n200\n"
    b"chain 2000 chr1 1000 + 0 200 chrB 1000 + 0 200 2\n200\n"
)
_TWO_PATH_CHAIN2 = (
    b"chain 500 chrA 1000 + 0 200 chrZ 1000 + 0 200 10\n200\n"
    b"chain 600 chrB 1000 + 0 200 chrZ 1000 + 100 300 11\n200\n"
)

# 自身片段相接场景：第一段把源区间拆成两个中间片段（chrA [0,50) 与
# [60,110)），第二段把它们压成最终坐标系里相接的 [0,50) 与 [50,100)
_ADJ_CHAIN1 = (
    b"chain 1000 chr1 1000 + 0 100 chrA 1000 + 0 110 1\n"
    b"50 0 10\n50\n"
)
_ADJ_CHAIN2 = (
    b"chain 500 chrA 1000 + 0 110 chrZ 1000 + 0 100 10\n"
    b"50 10 0\n50\n"
)

# 多染色体场景：chrA[0,200)->chrZ、chrA[300,500)->chrW
_TWO_CHROM_CHAIN2 = (
    b"chain 500 chrA 1000 + 0 200 chrZ 1000 + 0 200 10\n200\n"
    b"chain 600 chrA 1000 + 300 500 chrW 1000 + 0 200 11\n200\n"
)


@pytest.fixture()
def ident_token():
    """恒等两段 + 三行 BED：行1/行2 在 chrZ 上冲突，行3 无冲突。"""
    resp = _post2(b"chr1\t0\t100\nchr1\t50\t150\nchr1\t300\t400\n",
                  _IDENT_CHAIN1, _IDENT_CHAIN2)
    assert resp.status_code == 200, resp.text
    return resp.json()["download_token"]


@pytest.fixture()
def two_path_token():
    """两行各有两条路径（均为歧义行），最终都落在 chrZ。"""
    resp = _post2(b"chr1\t0\t100\nchr1\t50\t150\n",
                  _TWO_PATH_CHAIN1, _TWO_PATH_CHAIN2)
    assert resp.status_code == 200, resp.text
    rows = resp.json()["preview"]["rows"]
    assert [r["status"] for r in rows] == ["ambiguous", "ambiguous"]
    return resp.json()["download_token"]


# ---------- 基本审计（原始视图） ----------

def test_audit_original_view_basic_conflict(ident_token):
    resp = _audit(ident_token, "original")
    assert resp.status_code == 200, resp.text
    data = resp.json()

    assert data["download_token"] == ident_token
    assert data["view"] == "original"
    assert data["decision_version"] == 0
    assert data["audit_token"]
    assert data["download_url"] == f"/api/download/audit/{data['audit_token']}"
    assert ".audit-original.v0.tsv" in data["download_name"]

    # 行1 chrZ[0,100) 与行2 chrZ[50,150) 冲突 [50,100) 50bp；行3 无冲突
    assert data["summary"] == {
        "source_rows": 3, "units": 3,
        "conflict_pairs": 1, "conflict_regions": 1, "conflict_bases": 50,
    }
    assert data["conflicts"] == [{
        "chrom": "chrZ", "start": 50, "end": 100, "bases": 50,
        "a_source_line": 1, "a_path_index": 1, "a_chain1_id": "1", "a_chain2_id": "10",
        "b_source_line": 2, "b_path_index": 1, "b_chain1_id": "1", "b_chain2_id": "10",
    }]


def test_audit_view_defaults_to_original(ident_token):
    data = _audit(ident_token).json()
    assert data["view"] == "original"
    assert data["summary"]["conflict_regions"] == 1


def test_audit_tsv_download(ident_token):
    data = _audit(ident_token, "original").json()
    dl = client.get(data["download_url"])
    assert dl.status_code == 200
    assert dl.headers["content-type"].startswith("text/tab-separated-values")
    assert data["download_name"] in dl.headers["content-disposition"]
    assert dl.content.startswith(b"\xef\xbb\xbf")

    lines = _tsv_lines(dl)
    header = lines[0].split("\t")
    assert header == [
        "chrom", "overlap_start", "overlap_end", "overlap_bases",
        "a_source_line", "a_path_index", "a_chain1_id", "a_chain2_id",
        "b_source_line", "b_path_index", "b_chain1_id", "b_chain2_id",
    ]
    body = [ln.split("\t") for ln in lines[1:]]
    assert body == [[
        "chrZ", "50", "100", "50", "1", "1", "1", "10", "2", "1", "1", "10",
    ]]


# ---------- 同一源行的多条路径不互相比较 ----------

def test_audit_same_row_paths_never_compared():
    # 单行双路径，两条路径最终区间在 chrZ 上重叠 [100,200)∩[100,200)；
    # 但它们是同一源行的互斥候选，不产生任何冲突记录
    resp = _post2(b"chr1\t0\t200\n", _TWO_PATH_CHAIN1, _TWO_PATH_CHAIN2)
    assert resp.status_code == 200, resp.text
    row = resp.json()["preview"]["rows"][0]
    assert row["status"] == "ambiguous"
    assert len(row["candidates"]) == 2

    data = _audit(resp.json()["download_token"], "original").json()
    assert data["summary"]["units"] == 2
    assert data["summary"]["conflict_pairs"] == 0
    assert data["summary"]["conflict_regions"] == 0
    assert data["conflicts"] == []


def test_audit_original_view_two_path_rows(two_path_token):
    data = _audit(two_path_token, "original").json()
    # 行1: p1 chrZ[0,100) p2 chrZ[100,200)；行2: p1 chrZ[50,150) p2 chrZ[150,250)
    # 跨行交集：[50,100)（1p1×2p1）、[100,150)（1p2×2p1）、[150,200)（1p2×2p2）
    assert data["summary"] == {
        "source_rows": 2, "units": 4,
        "conflict_pairs": 3, "conflict_regions": 3, "conflict_bases": 150,
    }
    got = [(c["start"], c["end"], c["a_source_line"], c["a_path_index"],
            c["b_source_line"], c["b_path_index"]) for c in data["conflicts"]]
    assert got == [
        (50, 100, 1, 1, 2, 1),
        (100, 150, 1, 2, 2, 1),
        (150, 200, 1, 2, 2, 2),
    ]
    # 双方 chain ID 随路径正确带出
    c0 = data["conflicts"][0]
    assert (c0["a_chain1_id"], c0["a_chain2_id"]) == ("1", "10")
    assert (c0["b_chain1_id"], c0["b_chain2_id"]) == ("1", "10")
    c2 = data["conflicts"][2]
    assert (c2["a_chain1_id"], c2["a_chain2_id"]) == ("2", "11")
    assert (c2["b_chain1_id"], c2["b_chain2_id"]) == ("2", "11")


# ---------- 裁决视图：仅收窄已裁决的歧义行 ----------

def test_audit_adjudicated_view_unresolved_keeps_all_paths(two_path_token):
    # 未做任何裁决时，裁决视图与原始视图一致（未裁决行保留全部路径）
    data = _audit(two_path_token, "adjudicated").json()
    assert data["view"] == "adjudicated"
    assert data["decision_version"] == 0
    assert data["summary"]["conflict_regions"] == 3


def test_audit_adjudicated_view_narrows_only_resolved_rows(two_path_token):
    # 裁决行1选路径1（version 1）；行2未裁决，仍保留全部路径
    assert _decide(two_path_token, 0, 1, 1).status_code == 200
    data = _audit(two_path_token, "adjudicated").json()
    assert data["decision_version"] == 1
    # 只剩行1p1 chrZ[0,100) × 行2p1 chrZ[50,150) = [50,100)
    assert data["summary"]["units"] == 3
    assert data["summary"]["conflict_regions"] == 1
    assert data["conflicts"][0]["start"] == 50
    assert data["conflicts"][0]["end"] == 100
    assert (data["conflicts"][0]["a_source_line"],
            data["conflicts"][0]["a_path_index"]) == (1, 1)

    # 再裁决行2选路径2（version 2）：行1p1 [0,100) 与行2p2 [150,250) 不相交
    assert _decide(two_path_token, 1, 2, 2).status_code == 200
    data = _audit(two_path_token, "adjudicated").json()
    assert data["decision_version"] == 2
    assert data["summary"]["units"] == 2
    assert data["summary"]["conflict_regions"] == 0
    assert data["conflicts"] == []


def test_audit_adjudicated_filename_carries_version(two_path_token):
    _decide(two_path_token, 0, 1, 2)
    data = _audit(two_path_token, "adjudicated").json()
    assert ".audit-adjudicated.v1.tsv" in data["download_name"]


# ---------- 审计产物固定于发起时的裁决版本 ----------

def test_audit_result_pinned_to_initiation_version(two_path_token):
    # v0（未裁决）时发起审计
    first = _audit(two_path_token, "adjudicated").json()
    assert first["decision_version"] == 0
    assert first["summary"]["conflict_regions"] == 3

    # 之后裁决两行，使冲突消失
    _decide(two_path_token, 0, 1, 1)
    _decide(two_path_token, 1, 2, 2)
    second = _audit(two_path_token, "adjudicated").json()
    assert second["decision_version"] == 2
    assert second["summary"]["conflict_regions"] == 0

    # 旧审计产物不随后续裁决改变：仍是最初的 3 条明细
    dl_first = client.get(first["download_url"])
    assert dl_first.status_code == 200
    assert len(_tsv_lines(dl_first)) == 1 + 3
    dl_second = client.get(second["download_url"])
    assert len(_tsv_lines(dl_second)) == 1

    # 继续改裁决，旧产物依然不变
    _decide(two_path_token, 2, 1, None)
    assert client.get(first["download_url"]).content == dl_first.content
    assert client.get(second["download_url"]).content == dl_second.content


def test_audit_does_not_modify_original_job(two_path_token):
    tsv_before = client.get(f"/api/download/{two_path_token}").content
    adj_before = client.get(f"/api/download/{two_path_token}/adjudicated").content

    resp = _audit(two_path_token, "adjudicated")
    assert resp.status_code == 200

    assert client.get(f"/api/download/{two_path_token}").content == tsv_before
    assert client.get(f"/api/download/{two_path_token}/adjudicated").content == adj_before
    assert client.get(f"/api/result/{two_path_token}/decisions").json()["version"] == 0


# ---------- 自身目标片段先合并再求交 ----------

def test_audit_merges_own_fragments_before_intersection():
    # 行1 -> chrZ [0,50)+[50,100)（相接）；行2 -> chrZ [25,50)+[50,75)（相接）
    resp = _post2(b"chr1\t0\t100\nchr1\t25\t75\n", _ADJ_CHAIN1, _ADJ_CHAIN2)
    assert resp.status_code == 200, resp.text
    token = resp.json()["download_token"]
    rows = resp.json()["preview"]["rows"]
    # 每行一条路径、两个最终片段（确认片段确实被拆开）
    assert [len(r["candidates"][0]["fragments"]) for r in rows] == [2, 2]

    data = _audit(token, "original").json()
    # 合并后为 [0,100) ∩ [25,75) = [25,75) 一条记录；不合并会得到两条相接记录
    assert data["summary"]["conflict_regions"] == 1
    assert data["summary"]["conflict_bases"] == 50
    c = data["conflicts"][0]
    assert (c["start"], c["end"], c["bases"]) == (25, 75, 50)


# ---------- 排序：染色体 → 起止 → 双方行号/路径序号 ----------

def test_audit_results_sorted_by_chrom_then_interval():
    bed = b"chr1\t0\t100\nchr1\t50\t150\nchr1\t300\t400\nchr1\t350\t450\n"
    resp = _post2(bed, _IDENT_CHAIN1, _TWO_CHROM_CHAIN2)
    assert resp.status_code == 200, resp.text
    token = resp.json()["download_token"]

    data = _audit(token, "original").json()
    # chrW（字典序在前）的冲突排在 chrZ 之前
    got = [(c["chrom"], c["start"], c["end"],
            c["a_source_line"], c["b_source_line"]) for c in data["conflicts"]]
    assert got == [
        ("chrW", 50, 100, 3, 4),
        ("chrZ", 50, 100, 1, 2),
    ]


# ---------- 无冲突：空明细 + 零汇总 ----------

def test_audit_no_conflicts_empty_detail_zero_summary():
    resp = _post2(b"chr1\t0\t100\nchr1\t200\t300\n", _IDENT_CHAIN1, _IDENT_CHAIN2)
    assert resp.status_code == 200, resp.text
    token = resp.json()["download_token"]

    data = _audit(token, "original").json()
    assert data["conflicts"] == []
    assert data["summary"] == {
        "source_rows": 2, "units": 2,
        "conflict_pairs": 0, "conflict_regions": 0, "conflict_bases": 0,
    }

    dl = client.get(data["download_url"])
    assert dl.status_code == 200
    assert len(_tsv_lines(dl)) == 1          # 仅表头


def test_audit_unmapped_rows_not_audited():
    # 行1 有映射；行2 第一段无映射（染色体不存在），没有最终坐标，不参与
    resp = _post2(b"chr1\t0\t100\nchrX\t0\t100\n", _IDENT_CHAIN1, _IDENT_CHAIN2)
    assert resp.status_code == 200, resp.text
    token = resp.json()["download_token"]
    data = _audit(token, "original").json()
    assert data["summary"]["source_rows"] == 1
    assert data["summary"]["units"] == 1
    assert data["summary"]["conflict_regions"] == 0


# ---------- 错误：非双段 / 过期 / 非法参数 ----------

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
    resp = _audit(token, "original")
    assert resp.status_code == 400
    assert "双段" in resp.json()["errors"][0]
    assert len(main_module._audit_jobs) == audits_before


def test_audit_unknown_token_404():
    resp = _audit("no-such-token", "original")
    assert resp.status_code == 404
    assert "重新提交转换" in resp.json()["errors"][0]


def test_audit_expired_token_404(ident_token):
    main_module._jobs[ident_token].created_at = 0   # 模拟任务已过期
    resp = _audit(ident_token, "original")
    assert resp.status_code == 404
    assert "过期" in resp.json()["errors"][0]


@pytest.mark.parametrize("view", ["bogus", "Original", "", 1, True, None, ["original"]])
def test_audit_invalid_view_rejected(ident_token, view):
    audits_before = len(main_module._audit_jobs)
    resp = client.post(f"/api/result/{ident_token}/audit", json={"view": view})
    assert resp.status_code == 400
    assert "view" in resp.json()["errors"][0]
    assert len(main_module._audit_jobs) == audits_before


def test_audit_body_not_json(ident_token):
    resp = client.post(
        f"/api/result/{ident_token}/audit",
        content=b"not-json",
        headers={"Content-Type": "application/json"},
    )
    assert resp.status_code == 400


def test_audit_download_unknown_token_404():
    resp = client.get("/api/download/audit/no-such-token")
    assert resp.status_code == 404
    assert "审计结果" in resp.json()["errors"][0]


def test_audit_download_expired_404(ident_token):
    data = _audit(ident_token, "original").json()
    main_module._audit_jobs[data["audit_token"]].created_at = 0
    resp = client.get(data["download_url"])
    assert resp.status_code == 404
    assert "审计结果" in resp.json()["errors"][0]


# ---------- 大文件：扫描线，不做全量两两比较 ----------

def test_audit_sweep_line_avoids_full_pairwise(monkeypatch):
    # 400 行链式重叠（每行只与相邻 1~2 行相交）：全量两两需 79800 次
    # 区间求交，扫描线只需与活动窗口内的区间求交（约 2N 次）
    n = 400
    bed = "".join(f"chr1\t{i * 50}\t{i * 50 + 100}\n" for i in range(n)).encode()
    resp = _post2(bed, _IDENT_CHAIN1, _IDENT_CHAIN2)
    assert resp.status_code == 200, resp.text
    token = resp.json()["download_token"]

    calls = 0
    real_overlap = audit_module._interval_overlap

    def counting_overlap(a, b):
        nonlocal calls
        calls += 1
        return real_overlap(a, b)

    monkeypatch.setattr(audit_module, "_interval_overlap", counting_overlap)
    data = _audit(token, "original").json()

    assert data["summary"]["units"] == n
    # 链式重叠的冲突明细本身约 2N 条，求交调用同阶；远小于 N*(N-1)/2
    assert data["summary"]["conflict_regions"] > 0
    assert calls < 4 * n
    assert calls < n * (n - 1) // 2


def test_audit_many_disjoint_rows_fast_and_empty():
    # 500 行互不重叠：扫描窗口始终为空，不做任何区间求交
    bed = "".join(f"chr1\t{i * 1000}\t{i * 1000 + 100}\n" for i in range(500)).encode()
    resp = _post2(bed, _IDENT_CHAIN1, _IDENT_CHAIN2)
    assert resp.status_code == 200, resp.text
    data = _audit(resp.json()["download_token"], "original").json()
    assert data["summary"]["units"] == 500
    assert data["summary"]["conflict_regions"] == 0
    assert data["conflicts"] == []
