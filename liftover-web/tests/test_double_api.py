"""双段转换 HTTP 接口测试：预览、下载、失败阶段、错误拦截、裁决。"""

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
    side1: str = "target",
    side2: str = "target",
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
            data={"side1": side1, "side2": side2},
        )


@pytest.fixture()
def job2():
    resp = _post2()
    assert resp.status_code == 200, resp.text
    data = resp.json()
    return data["download_token"], data["preview"]["rows"]


def _rows_by_line(rows):
    return {r["line_no"]: r for r in rows}


def _tsv_lines(resp):
    return resp.content.decode("utf-8-sig").splitlines()


# ---------- 预览 ----------

def test_double_preview_counts_and_kind(job2):
    token, rows = job2
    by_line = _rows_by_line(rows)
    assert len(rows) == 8

    # 行3 完全映射：chr1:0-100 -> chrA:0-100 -> chrZ:0-100
    r3 = by_line[3]
    assert r3["status"] == "mapped"
    assert r3["failed_stage"] == ""
    path = r3["candidates"][0]
    assert (path["chain1_id"], path["chain2_id"]) == ("1", "10")
    assert path["path_index"] == 1
    assert (path["mid_chrom"], path["dest_chrom"]) == ("chrA", "chrZ")
    assert (path["seg1_covered"], path["seg2_covered"], path["total"]) == (100, 100, 100)
    frag = path["fragments"][0]
    assert (frag["src_start"], frag["src_end"]) == (0, 100)
    assert (frag["mid_start"], frag["mid_end"]) == (0, 100)
    assert (frag["dest_start"], frag["dest_end"]) == (0, 100)

    # 行4 部分映射：第二段缺口，未覆盖源区间 [100,120)
    r4 = by_line[4]
    assert r4["status"] == "partial"
    p4 = r4["candidates"][0]
    assert (p4["seg1_covered"], p4["seg2_covered"]) == (120, 100)
    assert [list(iv) for iv in p4["uncovered_src"]] == [[100, 120]]

    # 行5 部分映射：第一段缺口 + 第二段未覆盖，仅源 150-160 完成两段
    r5 = by_line[5]
    assert r5["status"] == "partial"
    p5 = r5["candidates"][0]
    assert (p5["seg1_covered"], p5["seg2_covered"]) == (30, 10)
    assert [list(iv) for iv in p5["uncovered_src"]] == [[100, 150]]

    # 行6 歧义：两条可达路径全部保留
    r6 = by_line[6]
    assert r6["status"] == "ambiguous"
    assert len(r6["candidates"]) == 2
    assert {p["dest_chrom"] for p in r6["candidates"]} == {"chrW", "chrV"}
    assert [p["path_index"] for p in r6["candidates"]] == [1, 2]
    assert {p["chain2_id"] for p in r6["candidates"]} == {"11", "12"}

    # 行7/8 第二段无映射；行9/10 第一段无映射
    assert by_line[7]["failed_stage"] == "segment2"
    assert "chrM" in by_line[7]["reason"]
    assert by_line[8]["failed_stage"] == "segment2"
    assert by_line[9]["failed_stage"] == "segment1"
    assert by_line[10]["failed_stage"] == "segment1"
    assert "chrX" in by_line[10]["reason"]


def test_double_response_kind_and_name(job2):
    token, rows = job2
    resp = _post2()
    data = resp.json()
    assert data["kind"] == "double"
    assert data["download_name"].startswith("liftover_double_target-target_")
    assert data["preview"]["counts"] == {
        "mapped": 1, "partial": 2, "ambiguous": 1, "unmapped": 4,
    }


# ---------- 原始 TSV 下载 ----------

def test_double_tsv_download(job2):
    token, rows = job2
    dl = client.get(f"/api/download/{token}")
    assert dl.status_code == 200
    lines = _tsv_lines(dl)
    header = lines[0].split("\t")
    assert len(header) == 25
    assert header[:3] == ["source_line", "status", "failed_stage"]
    assert "chain1_id" in header and "chain2_id" in header
    assert "path_index" in header
    assert "seg1_covered" in header and "seg2_covered" in header
    assert "uncovered_src" in header

    body = [ln.split("\t") for ln in lines[1:]]
    # 行3 完全映射：路径与两段 chain ID 齐全
    r3 = [f for f in body if f[0] == "3"]
    assert len(r3) == 1
    f = r3[0]
    assert f[1] == "mapped" and f[2] == "."
    assert (f[6], f[7], f[8]) == ("1", "1", "10")
    assert (f[9], f[11]) == ("chrA", "chrZ")
    assert (f[14], f[15], f[16], f[17], f[18], f[19]) == ("0", "100", "0", "100", "0", "100")
    assert (f[20], f[21], f[22]) == ("100", "100", "100")
    assert f[23] == "."

    # 行4 部分映射：uncovered_src 记录未覆盖源区间
    r4 = [f for f in body if f[0] == "4"]
    assert r4[0][1] == "partial"
    assert r4[0][23] == "[100,120)"

    # 行6 歧义：两条路径各占一行
    r6 = [f for f in body if f[0] == "6"]
    assert len(r6) == 2
    assert {f[11] for f in r6} == {"chrW", "chrV"}
    assert {f[6] for f in r6} == {"1", "2"}

    # 无映射行：保留物理行号与失败阶段，未覆盖为整个源区间
    r7 = [f for f in body if f[0] == "7"]
    assert r7[0][1] == "unmapped" and r7[0][2] == "segment2"
    assert r7[0][6] == "." and r7[0][23] == "[500,600)"
    r9 = [f for f in body if f[0] == "9"]
    assert r9[0][2] == "segment1"
    assert r9[0][23] == "[900,1000)"


# ---------- 输入错误：带行号拒绝，不产生下载 ----------

def test_double_bad_bed_rejected_no_download():
    before = len(main_module._jobs)
    resp = _post2(bed_name="example_bad.bed")
    assert resp.status_code == 400
    errors = resp.json()["errors"]
    assert any("第 3 行" in e for e in errors)
    assert any("第 4 行" in e for e in errors)
    assert any("第 5 行" in e for e in errors)
    assert "download_token" not in resp.json()
    assert len(main_module._jobs) == before


def test_double_bad_chain1_rejected_with_file_label():
    bad_chain = b"chain 1000 chr1 100 + 0 10\n10\n"  # 头只有 7 列
    with open(os.path.join(EXAMPLES, "example_double.bed"), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo_seg2.chain"), "rb") as c2:
        resp = client.post(
            "/api/liftover2",
            files={
                "bed": ("b.bed", bf, "text/plain"),
                "chain1": ("bad.chain", bad_chain, "text/plain"),
                "chain2": ("c2.chain", c2, "text/plain"),
            },
            data={"side1": "target", "side2": "target"},
        )
    assert resp.status_code == 400
    errors = resp.json()["errors"]
    assert any("chain 文件 1" in e and "第 1 行" in e for e in errors)


def test_double_bad_chain2_rejected_with_file_label():
    bad_chain = b"chain 1000 chrA 500 + 0 200 chrZ 1000 + 0 200 10\n150\n"
    with open(os.path.join(EXAMPLES, "example_double.bed"), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo_seg1.chain"), "rb") as c1:
        resp = client.post(
            "/api/liftover2",
            files={
                "bed": ("b.bed", bf, "text/plain"),
                "chain1": ("c1.chain", c1, "text/plain"),
                "chain2": ("bad.chain", bad_chain, "text/plain"),
            },
            data={"side1": "target", "side2": "target"},
        )
    assert resp.status_code == 400
    errors = resp.json()["errors"]
    assert any("chain 文件 2" in e and "第 2 行" in e for e in errors)


def test_double_invalid_sides_rejected():
    for data in ({"side1": "up", "side2": "target"},
                 {"side1": "target", "side2": "down"}):
        with open(os.path.join(EXAMPLES, "example_double.bed"), "rb") as bf, \
             open(os.path.join(EXAMPLES, "demo_seg1.chain"), "rb") as c1, \
             open(os.path.join(EXAMPLES, "demo_seg2.chain"), "rb") as c2:
            resp = client.post(
                "/api/liftover2",
                files={
                    "bed": ("b.bed", bf, "text/plain"),
                    "chain1": ("c1.chain", c1, "text/plain"),
                    "chain2": ("c2.chain", c2, "text/plain"),
                },
                data=data,
            )
        assert resp.status_code == 400


# ---------- 双段结果的歧义裁决 ----------

def _ambiguous_line2(rows):
    amb = [r for r in rows if r["status"] == "ambiguous"]
    assert len(amb) == 1
    return amb[0]["line_no"]


def _decide(token, version, source_line, candidate_index):
    return client.post(
        f"/api/result/{token}/decisions",
        json={"version": version, "source_line": source_line,
              "candidate_index": candidate_index},
    )


def test_double_decisions_pick_and_unpick(job2):
    token, rows = job2
    line = _ambiguous_line2(rows)

    got = client.get(f"/api/result/{token}/decisions").json()
    assert got["version"] == 0
    assert got["decisions"] == [
        {"source_line": line, "candidate_index": None, "candidate_count": 2}
    ]

    # 选路径 2（chrV）
    resp = _decide(token, 0, line, 2)
    assert resp.status_code == 200
    assert resp.json()["version"] == 1

    dl = _tsv_lines(client.get(f"/api/download/{token}/adjudicated"))
    header = dl[0].split("\t")
    assert len(header) == 26 and header[-1] == "decision"
    amb = [ln.split("\t") for ln in dl[1:]
           if ln.startswith(f"{line}\tambiguous\t")]
    assert len(amb) == 1
    assert amb[0][6] == "2"                 # path_index 仍是原序号
    assert amb[0][11] == "chrV"
    assert amb[0][1] == "ambiguous"         # 原始 status 不变
    assert amb[0][-1] == "已裁决"

    # 撤销：恢复全部路径
    assert _decide(token, 1, line, None).json()["version"] == 2
    dl = _tsv_lines(client.get(f"/api/download/{token}/adjudicated"))
    amb = [ln for ln in dl[1:] if ln.startswith(f"{line}\tambiguous\t")]
    assert len(amb) == 2
    assert all(ln.endswith("\t未裁决") for ln in amb)


def test_double_original_download_byte_identical_after_decisions(job2):
    token, rows = job2
    line = _ambiguous_line2(rows)
    before = client.get(f"/api/download/{token}").content
    _decide(token, 0, line, 1)
    _decide(token, 1, line, None)
    after = client.get(f"/api/download/{token}").content
    assert before == after
    text = before.decode("utf-8-sig")
    assert len(text.splitlines()[0].split("\t")) == 25


def test_double_adjudicated_other_rows_identical_to_original(job2):
    token, _ = job2
    orig = [ln.split("\t") for ln in _tsv_lines(client.get(f"/api/download/{token}"))]
    adj = [ln.split("\t") for ln in _tsv_lines(
        client.get(f"/api/download/{token}/adjudicated"))]
    assert [r[:25] for r in orig] == [r[:25] for r in adj]
    partial = [r for r in adj if r[1] == "partial"]
    assert partial and all(r[-1] == "." for r in partial)


# ---------- 单段接口保持兼容 ----------

def test_single_segment_still_works_and_reports_kind():
    with open(os.path.join(EXAMPLES, "example_target.bed"), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo.chain"), "rb") as cf:
        resp = client.post(
            "/api/liftover",
            files={"bed": ("b.bed", bf, "text/plain"),
                   "chain": ("c.chain", cf, "text/plain")},
            data={"side": "target"},
        )
    assert resp.status_code == 200
    data = resp.json()
    assert data["kind"] == "single"
    assert data["preview"]["counts"] == {
        "mapped": 1, "partial": 1, "ambiguous": 1, "unmapped": 3,
    }
