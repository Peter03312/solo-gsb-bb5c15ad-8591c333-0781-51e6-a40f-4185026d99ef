"""歧义人工裁决测试：选择/撤销、版本并发、拒绝场景、裁决版 TSV。"""

import os

import pytest
from fastapi.testclient import TestClient

from app import main as main_module
from app.main import app

client = TestClient(app)

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples")


@pytest.fixture()
def job():
    """跑一次示例转换，返回 (token, 预览 rows)。示例歧义行是 BED 第 11 行。"""
    with open(os.path.join(EXAMPLES, "example_target.bed"), "rb") as bf, \
         open(os.path.join(EXAMPLES, "demo.chain"), "rb") as cf:
        resp = client.post(
            "/api/liftover",
            files={"bed": ("example_target.bed", bf, "text/plain"),
                   "chain": ("demo.chain", cf, "text/plain")},
            data={"side": "target"},
        )
    assert resp.status_code == 200, resp.text
    data = resp.json()
    return data["download_token"], data["preview"]["rows"]


def _ambiguous_line(rows):
    amb = [r for r in rows if r["status"] == "ambiguous"]
    assert len(amb) == 1
    assert amb[0]["line_no"] == 11
    return amb[0]["line_no"]


def _decide(token, version, source_line, candidate_index):
    return client.post(
        f"/api/result/{token}/decisions",
        json={"version": version, "source_line": source_line,
              "candidate_index": candidate_index},
    )


def _tsv_lines(resp):
    return resp.content.decode("utf-8-sig").splitlines()


# ---------- 查看当前选择 ----------

def test_initial_decisions_empty_version_zero(job):
    token, rows = job
    line = _ambiguous_line(rows)
    resp = client.get(f"/api/result/{token}/decisions")
    assert resp.status_code == 200
    data = resp.json()
    assert data["version"] == 0
    assert data["decisions"] == [
        {"source_line": line, "candidate_index": None, "candidate_count": 2}
    ]


# ---------- 裁决版 TSV：未裁决 ----------

def test_adjudicated_tsv_pending_marks_all_candidates(job):
    token, rows = job
    line = _ambiguous_line(rows)
    resp = client.get(f"/api/download/{token}/adjudicated")
    assert resp.status_code == 200
    lines = _tsv_lines(resp)
    header = lines[0].split("\t")
    assert header[-1] == "decision"
    assert len(header) == 16

    amb = [ln for ln in lines[1:] if ln.startswith(f"{line}\tambiguous\t")]
    assert len(amb) == 2                      # 两个候选都保留
    for ln in amb:
        assert ln.split("\t")[-1] == "未裁决"

    # 非歧义行照原结果输出，decision 为 '.'
    mapped = [ln for ln in lines[1:] if "\tmapped\t" in ln]
    assert mapped and all(ln.endswith("\t.") for ln in mapped)
    unmapped = [ln for ln in lines[1:] if "\tunmapped\t" in ln]
    assert len(unmapped) == 3
    assert all(ln.endswith("\t.") for ln in unmapped)


# ---------- 选定候选 ----------

def test_pick_candidate_outputs_only_selected_fragments(job):
    token, rows = job
    line = _ambiguous_line(rows)
    resp = _decide(token, 0, line, 1)
    assert resp.status_code == 200, resp.text
    data = resp.json()
    assert data["version"] == 1
    assert data["modified"] is True

    # GET 反映选择
    got = client.get(f"/api/result/{token}/decisions").json()
    assert got["version"] == 1
    assert got["decisions"][0]["candidate_index"] == 1

    dl = _tsv_lines(client.get(f"/api/download/{token}/adjudicated"))
    amb = [ln for ln in dl[1:] if ln.startswith(f"{line}\tambiguous\t")]
    assert len(amb) == 1                      # 只剩被选候选的片段
    fields = amb[0].split("\t")
    assert fields[5] == "1"                  # candidate_index 仍是原序号
    assert fields[7] == "chrA"
    assert fields[1] == "ambiguous"          # 原始 status 未被改成 mapped
    assert fields[-1] == "已裁决"


def test_pick_candidate_two_keeps_original_status_and_coords(job):
    token, rows = job
    line = _ambiguous_line(rows)
    assert _decide(token, 0, line, 2).json()["version"] == 1
    dl = _tsv_lines(client.get(f"/api/download/{token}/adjudicated"))
    amb = [ln for ln in dl[1:] if ln.startswith(f"{line}\tambiguous\t")]
    assert len(amb) == 1
    fields = amb[0].split("\t")
    assert fields[5] == "2"
    assert fields[7] == "chrB"
    # 负链候选：零基左闭右开坐标不变
    assert (fields[12], fields[13]) == ("240", "280")
    assert int(fields[12]) < int(fields[13])
    assert fields[-1] == "已裁决"


# ---------- 撤销 ----------

def test_unpick_restores_all_candidates_and_bumps_version(job):
    token, rows = job
    line = _ambiguous_line(rows)
    assert _decide(token, 0, line, 1).json()["version"] == 1
    resp = _decide(token, 1, line, None)
    assert resp.status_code == 200
    assert resp.json()["version"] == 2
    assert resp.json()["candidate_index"] is None

    dl = _tsv_lines(client.get(f"/api/download/{token}/adjudicated"))
    amb = [ln for ln in dl[1:] if ln.startswith(f"{line}\tambiguous\t")]
    assert len(amb) == 2
    assert all(ln.endswith("\t未裁决") for ln in amb)

    got = client.get(f"/api/result/{token}/decisions").json()
    assert got["decisions"][0]["candidate_index"] is None


def test_repick_different_candidate_replaces_selection(job):
    token, rows = job
    line = _ambiguous_line(rows)
    _decide(token, 0, line, 1)
    resp = _decide(token, 1, line, 2)
    assert resp.json()["version"] == 2
    dl = _tsv_lines(client.get(f"/api/download/{token}/adjudicated"))
    amb = [ln for ln in dl[1:] if ln.startswith(f"{line}\tambiguous\t")]
    assert len(amb) == 1
    assert amb[0].split("\t")[7] == "chrB"


# ---------- 版本与并发 ----------

def test_stale_version_rejected_with_409(job):
    token, rows = job
    line = _ambiguous_line(rows)
    assert _decide(token, 0, line, 1).status_code == 200
    # 仍拿旧版本 0 提交
    resp = _decide(token, 0, line, 2)
    assert resp.status_code == 409
    body = resp.json()
    assert "版本" in body["errors"][0]
    assert body["current_version"] == 1
    # 旧版本提交没有改动选择
    dl = _tsv_lines(client.get(f"/api/download/{token}/adjudicated"))
    amb = [ln for ln in dl[1:] if ln.startswith(f"{line}\tambiguous\t")]
    assert len(amb) == 1 and amb[0].split("\t")[7] == "chrA"


def test_same_version_concurrent_only_one_succeeds(job):
    token, rows = job
    line = _ambiguous_line(rows)

    async def scenario():
        import asyncio
        from httpx import ASGITransport, AsyncClient
        transport = ASGITransport(app=app)
        async with AsyncClient(transport=transport, base_url="http://t") as ac:
            payload = {"version": 0, "source_line": line, "candidate_index": 1}
            results = await asyncio.gather(
                ac.post(f"/api/result/{token}/decisions", json=payload),
                ac.post(f"/api/result/{token}/decisions", json=payload),
            )
            return [r.status_code for r in results]

    import asyncio
    codes = asyncio.run(scenario())
    assert sorted(codes) == [200, 409]
    assert client.get(f"/api/result/{token}/decisions").json()["version"] == 1


def test_idempotent_submission_keeps_version(job):
    token, rows = job
    line = _ambiguous_line(rows)
    _decide(token, 0, line, 1)
    # 相同选择重复提交：不递增版本
    resp = _decide(token, 1, line, 1)
    assert resp.status_code == 200
    body = resp.json()
    assert body["version"] == 1
    assert body["modified"] is False
    # 撤销状态下重复撤销同样不递增
    _decide(token, 1, line, None)
    again = _decide(token, 2, line, None)
    assert again.json()["version"] == 2
    assert again.json()["modified"] is False


# ---------- 明确拒绝的场景 ----------

def test_decide_non_ambiguous_line_rejected(job):
    token, rows = job
    # 第 8 行 mapped、第 9 行 partial、第 10 行 unmapped
    for line, status in ((8, "mapped"), (9, "partial"), (10, "unmapped")):
        resp = _decide(token, 0, line, 1)
        assert resp.status_code == 400, resp.text
        assert "歧义" in resp.json()["errors"][0] or status in resp.json()["errors"][0]


def test_decide_nonexistent_line_rejected(job):
    token, _ = job
    resp = _decide(token, 0, 999, 1)
    assert resp.status_code == 400
    assert "不存在" in resp.json()["errors"][0]


def test_decide_candidate_index_out_of_range(job):
    token, rows = job
    line = _ambiguous_line(rows)
    for bad in (0, 3, -1):
        resp = _decide(token, 0, line, bad)
        assert resp.status_code == 400
        assert "越界" in resp.json()["errors"][0]
    assert client.get(f"/api/result/{token}/decisions").json()["version"] == 0


def test_decide_bad_payload_shapes(job):
    token, rows = job
    line = _ambiguous_line(rows)
    base = f"/api/result/{token}/decisions"

    assert client.post(base, content="not json",
                       headers={"Content-Type": "application/json"}).status_code == 400
    assert client.post(base, json={"source_line": line,
                                  "candidate_index": 1}).status_code == 400  # 缺 version
    assert client.post(base, json={"version": 0,
                                  "candidate_index": 1}).status_code == 400  # 缺 source_line
    assert client.post(base, json={"version": 0, "source_line": line}).status_code == 400
    assert client.post(base, json={"version": "0", "source_line": line,
                                  "candidate_index": 1}).status_code == 400
    assert client.post(base, json={"version": 0, "source_line": line,
                                  "candidate_index": True}).status_code == 400
    assert client.post(base, json=[1, 2, 3]).status_code == 400


def test_expired_result_rejected(monkeypatch, job):
    token, rows = job
    line = _ambiguous_line(rows)
    monkeypatch.setattr(main_module, "JOB_TTL_SECONDS", 0)
    assert _decide(token, 0, line, 1).status_code == 404
    assert client.get(f"/api/result/{token}/decisions").status_code == 404
    assert client.get(f"/api/download/{token}/adjudicated").status_code == 404
    monkeypatch.undo()
    # TTL 恢复后，过期任务已被清理，token 仍 404
    assert client.get(f"/api/download/{token}").status_code == 404


def test_unknown_token_404():
    assert client.get("/api/result/no-such-token/decisions").status_code == 404
    assert client.post(
        "/api/result/no-such-token/decisions",
        json={"version": 0, "source_line": 1, "candidate_index": 1},
    ).status_code == 404
    assert client.get("/api/download/no-such-token/adjudicated").status_code == 404


# ---------- 原始结果兼容：裁决前后字节不变 ----------

def test_original_download_byte_identical_after_decisions(job):
    token, rows = job
    line = _ambiguous_line(rows)
    before = client.get(f"/api/download/{token}").content
    _decide(token, 0, line, 1)
    _decide(token, 1, line, None)
    after = client.get(f"/api/download/{token}").content
    assert before == after
    # 原始 TSV 仍只有 15 列、两个歧义候选、status 原样
    text = before.decode("utf-8-sig")
    header = text.splitlines()[0].split("\t")
    assert len(header) == 15 and header[-1] == "reason"
    amb = [ln for ln in text.splitlines()[1:] if "\tambiguous\t" in ln]
    assert len(amb) == 2


def test_adjudicated_filename_carries_version(job):
    token, rows = job
    line = _ambiguous_line(rows)
    resp = client.get(f"/api/download/{token}/adjudicated")
    cd = resp.headers["content-disposition"]
    assert "adjudicated.v0.tsv" in cd
    _decide(token, 0, line, 2)
    resp2 = client.get(f"/api/download/{token}/adjudicated")
    assert "adjudicated.v1.tsv" in resp2.headers["content-disposition"]


def test_adjudicated_other_rows_identical_to_original(job):
    # 除歧义行外，裁决版前 15 列与原始 TSV 完全一致；partial 不被改称 mapped
    token, _ = job
    orig = [ln.split("\t") for ln in _tsv_lines(client.get(f"/api/download/{token}"))]
    adj = [ln.split("\t") for ln in _tsv_lines(
        client.get(f"/api/download/{token}/adjudicated"))]
    assert [r[:15] for r in orig] == [r[:15] for r in adj]
    partial_rows = [r for r in adj if r[1] == "partial"]
    assert partial_rows and all(r[-1] == "." for r in partial_rows)
