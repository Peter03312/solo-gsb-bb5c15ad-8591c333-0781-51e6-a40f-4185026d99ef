"""BED / chain 解析错误定位测试。"""

from app.bed import parse_bed
from app.chain import parse_chain


def test_bad_bed_reports_each_bad_line_number():
    text = "chr1\t56\t176\nchr1\t50\nchr1\tabc\t185\nchr1\t300\t200\n"
    result = parse_bed(text)
    assert len(result.records) == 1  # 第 1 行有效
    assert len(result.errors) == 3
    assert "第 2 行" in result.errors[0]
    assert "第 3 行" in result.errors[1]
    assert "第 4 行" in result.errors[2]
    assert "start (300) 必须小于 end (200)" in result.errors[2]


def test_bed_blank_and_comment_lines_are_skipped():
    text = "# comment\n\ntrack name=x\nchr1\t0\t10\textra\tcolumns\n"
    result = parse_bed(text)
    assert result.errors == []
    assert len(result.records) == 1
    assert result.records[0].line_no == 4
    assert result.records[0].extra == "extra\tcolumns"


def test_empty_bed_is_error():
    assert parse_bed("").errors
    assert parse_bed("# only comment\n").errors


def test_chain_header_column_count_error():
    bad = "chain 1000 chr1 100 + 0 10\n10\n"
    result = parse_chain(bad)
    assert result.errors
    assert "第 1 行" in result.errors[0]
    assert "13 列" in result.errors[0]
    assert result.chains == []


def test_chain_missing_terminator_reports_line():
    text = (
        "chain 1000 chr1 100 + 0 30 chrQ 50 + 0 30 1\n"
        "10\t10\t10\n"
        "10\t0\t0\n"   # dt/dq 同时为 0 非法（应直接写终止行）
    )
    result = parse_chain(text)
    assert any("第 3 行" in e and "同时为 0" in e for e in result.errors)


def test_chain_block_running_past_end_is_caught():
    text = (
        "chain 1000 chr1 100 + 0 20 chrQ 50 + 0 20 1\n"
        "15\t10\t0\n"      # target gap 后游标 25 > tEnd 20
    )
    result = parse_chain(text)
    assert result.errors
    assert any("越过 chain 末端" in e for e in result.errors)


def test_chain_zero_gap_on_penultimate_line_is_allowed():
    # 最后一个非终止行允许一侧 gap 为 0（target 紧贴末端、query 仍有 gap）
    # tEnd=20 qEnd=30：块1 t[0,10) q[0,10)；gap dt=5 dq=10；块2 t[15,20) q[20,25)；
    # gap dt=0 dq=5；终止块 t[20,25)？target 会越界，因此 qEnd 取 30、最后块 t[20,20) 不行。
    # 合法布局：tEnd=25 qEnd=30
    #   块1 t[0,10)  q[0,10)；10 5 10
    #   块2 t[15,20) q[20,20)？重新构造：
    #   块1 t[0,10)  q[0,10)；行 "10 5 5"
    #   块2 t[15,20) q[15,20)；行 "5 0 5"   <- target gap 0
    #   块3 t[20,25) q[25,30)；终止 "5"
    text = (
        "chain 1000 chr1 100 + 0 25 chrQ 50 + 0 30 1\n"
        "10\t5\t5\n"
        "5\t0\t5\n"
        "5\n"
    )
    result = parse_chain(text)
    assert result.errors == []
    assert len(result.chains[0].blocks) == 3


def test_chain_unterminated_at_eof_reports_header():
    text = "chain 1000 chr1 100 + 0 20 chrQ 50 + 0 20 1\n10\t10\t10\n"
    result = parse_chain(text)
    assert result.errors
    assert any("缺少单字段终止行" in e and "第 1 行" in e for e in result.errors)


def test_chain_terminator_not_reaching_end():
    # 单块 size=15 终止，但 tEnd/qEnd=20
    text = "chain 1000 chr1 100 + 0 20 chrQ 50 + 0 20 1\n15\n"
    result = parse_chain(text)
    assert result.errors
    assert any("终止块未对齐 tEnd/qEnd" in e for e in result.errors)


def test_chain_data_line_outside_header():
    text = "10 3 3\n"
    result = parse_chain(text)
    assert result.errors
    assert "第 1 行" in result.errors[0]


def test_valid_demo_chain_parses_two_chains():
    import os
    demo = os.path.join(os.path.dirname(__file__), "..", "examples", "demo.chain")
    result = parse_chain(open(demo, encoding="utf-8").read())
    assert result.errors == []
    assert [c.chain_id for c in result.chains] == ["1", "2"]
    assert len(result.chains[0].blocks) == 3
