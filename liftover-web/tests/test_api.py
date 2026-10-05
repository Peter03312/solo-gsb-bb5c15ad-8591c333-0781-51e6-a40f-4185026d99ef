"""HTTP 接口测试：成功转换、错误阻止导出、下载。"""

import os

import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples")


def _post(bed_name: str, side: str = "target", chain_name: str = "demo.chain"):
    with open(os.path.join(EXAMPLES, bed_name), "rb") as bf, \
         open(os.path.join(EXAMPLES, chain_name), "rb") as cf:
        return client.post(
            "/api/liftover",
            files={"bed": (bed_name, bf, "text/plain"),
                   "chain": (chain_name, cf, "text/plain")},
            data={"side": side},
        )


def test_health():
    assert client.get("/healthz").json() == {"status": "ok"}


def test_index_served():
    resp = client.get("/")
    assert resp.status_code == 200
    assert "离线" in resp.text


def test_successful_liftover_preview_and_download():
    resp = _post("example_target.bed")
    assert resp.status_code == 200, resp.text
    data = resp.json()
    preview = data["preview"]
    assert preview["total"] == 6
    assert preview["counts"] == {
        "mapped": 1, "partial": 1, "ambiguous": 1, "unmapped": 3,
    }
    token = data["download_token"]

    dl = client.get(f"/api/download/{token}")
    assert dl.status_code == 200
    body = dl.content.decode("utf-8-sig")
    header = body.splitlines()[0].split("\t")
    assert header[0] == "source_line"
    assert header[3] == "source_start"
    # 歧义行在 TSV 中有两行候选（全部候选都保留，不挑选）
    ambiguous_lines = [ln for ln in body.splitlines()[1:] if "\tambiguous\t" in ln]
    assert len(ambiguous_lines) == 2
    assert {ln.split("\t")[7] for ln in ambiguous_lines} == {"chrA", "chrB"}
    # 无映射行保留源行号与原因
    unmapped_lines = [ln for ln in body.splitlines() if "\tunmapped\t" in ln]
    assert len(unmapped_lines) == 3
    for ln in unmapped_lines:
        assert ln.split("\t")[5] == "."


def test_bad_bed_blocks_export_and_creates_no_download():
    resp = _post("example_bad.bed")
    assert resp.status_code == 400
    errors = resp.json()["errors"]
    assert any("第 3 行" in e for e in errors)
    assert any("第 4 行" in e for e in errors)
    assert any("第 5 行" in e for e in errors)


def test_query_side_liftover_minus_strand():
    resp = _post("example_query.bed", side="query")
    assert resp.status_code == 200, resp.text
    rows = resp.json()["preview"]["rows"]
    by_chrom = {(r["src_chrom"], r["src_start"]): r for r in rows}
    b = by_chrom[("chrB", 240)]
    assert b["status"] == "mapped"
    frag = b["candidates"][0]["fragments"][0]
    assert (frag["dest_start"], frag["dest_end"]) == (320, 360)
    assert b["candidates"][0]["dest_strand"] == "+"


def test_invalid_side_rejected():
    with open(os.path.join(EXAMPLES, "example_target.bed"), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo.chain"), "rb") as cf:
        resp = client.post(
            "/api/liftover",
            files={"bed": ("b.bed", bf, "text/plain"),
                   "chain": ("c.chain", cf, "text/plain")},
            data={"side": "upside"},
        )
    assert resp.status_code == 400


def test_unknown_download_token_404():
    assert client.get("/api/download/nope-not-a-token").status_code == 404
