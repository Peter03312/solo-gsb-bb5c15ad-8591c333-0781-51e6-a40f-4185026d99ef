"""引擎转换正确性测试。"""

from app.bed import BedRecord
from app.chain import Block, Chain
from app.engine import LiftOverEngine


def _chain_minus_q() -> Chain:
    # target [0,30) -> query 负链，qSize=150
    # 两个对齐块（存储帧）：
    #   t[0,10)   <-> q[0,10]    （正向帧 q[140,150)）
    #   t[20,30)  <-> q[10,20]   （正向帧 q[130,140)）
    # 存储帧块间 query gap = 0（10->10，正向帧表现为两块在 140 处相接的
    # 反向排列），target gap = 10，合法。
    return Chain(
        chain_id="1", t_name="chrT", t_size=100, t_strand="+",
        t_start=0, t_end=30,
        q_name="chrQ", q_size=150, q_strand="-",
        q_start=0, q_end=20, score=1000,
        blocks=[Block(0, 10, 0, 10), Block(20, 30, 10, 20)],
        header_line=1,
    )


def _rec(start: int, end: int, chrom: str = "chrT") -> BedRecord:
    return BedRecord(0, chrom, start, end, f"{chrom}\t{start}\t{end}", "")


def test_minus_strand_target_to_query_keeps_half_open():
    eng = LiftOverEngine([_chain_minus_q()], "target")
    row = eng.convert(_rec(0, 10))
    assert row.status == "mapped"
    frag = row.candidates[0].fragments[0]
    assert (frag.dest_start, frag.dest_end) == (140, 150)
    assert frag.dest_start < frag.dest_end          # 仍然左闭右开
    assert row.candidates[0].dest_strand == "-"


def test_minus_strand_gap_splits_into_ordered_fragments_descending_dest():
    eng = LiftOverEngine([_chain_minus_q()], "target")
    row = eng.convert(_rec(0, 30))
    assert row.status == "partial"
    assert "2 个有序片段" in row.reason
    assert len(row.candidates[0].fragments) == 2
    f1, f2 = row.candidates[0].fragments
    # 片段按源坐标有序
    assert (f1.src_start, f1.src_end) == (0, 10)
    assert (f2.src_start, f2.src_end) == (20, 30)
    # 目标侧负链：片段按目标坐标降序，每个仍 start < end
    assert (f1.dest_start, f1.dest_end) == (140, 150)
    assert (f2.dest_start, f2.dest_end) == (130, 140)
    assert f1.dest_start > f2.dest_start
    assert row.candidates[0].covered == 20


def test_query_to_target_on_minus_strand_maps_forward_half_open():
    eng = LiftOverEngine([_chain_minus_q()], "query")
    # 输入 query 正向帧 [140,150) -> target [0,10)，目标方向为 +
    row = eng.convert(_rec(140, 150, chrom="chrQ"))
    assert row.status == "mapped"
    frag = row.candidates[0].fragments[0]
    assert (frag.dest_start, frag.dest_end) == (0, 10)
    assert frag.dest_start < frag.dest_end
    assert row.candidates[0].dest_strand == "+"
    assert row.candidates[0].dest_chrom == "chrT"


def _chain_minus_q_offset() -> Chain:
    # 与 _chain_minus_q 相同的对齐，但 qStart=50（非 0）。
    # size 枢轴与 (qStart+qEnd) 枢轴会给出不同结果，用于锁定 UCSC 语义。
    # query 长 160，存储帧 [50,70)；归一化正向帧 [90,110)。
    #   t[0,10)  <-> stored q[50,60]  -> forward [100,110)
    #   t[20,30) <-> stored q[60,70]  -> forward [90,100)
    return Chain(
        chain_id="9", t_name="chrT", t_size=100, t_strand="+",
        t_start=0, t_end=30,
        q_name="chrQ2", q_size=160, q_strand="-",
        q_start=50, q_end=70, score=1000,
        blocks=[Block(0, 10, 50, 60), Block(20, 30, 60, 70)],
        header_line=1,
    )


def test_negative_strand_pivot_is_chromosome_size_not_span():
    eng = LiftOverEngine([_chain_minus_q_offset()], "target")
    row = eng.convert(_rec(0, 10, chrom="chrT"))
    assert row.status == "mapped"
    frag = row.candidates[0].fragments[0]
    # 若错误地用 (qStart+qEnd) 枢轴会得到 [150,160)，必须是 size 枢轴的 [100,110)
    assert (frag.dest_start, frag.dest_end) == (100, 110)

    eng_rev = LiftOverEngine([_chain_minus_q_offset()], "query")
    back = eng_rev.convert(BedRecord(0, "chrQ2", 100, 110, "chrQ2\t100\t110", ""))
    assert back.status == "mapped"
    assert (back.candidates[0].fragments[0].dest_start,
            back.candidates[0].fragments[0].dest_end) == (0, 10)


def test_double_minus_strand_round_trip():
    # target 与 query 均为负链：相对同向，块为平移
    c = Chain(
        "1", "chrT", 100, "-", 10, 30,
        "chrQ", 80, "-", 20, 40, 1000,
        [Block(10, 30, 20, 40)], 1,
    )
    # 存储帧 t[10,30)/q[20,40) -> 正向帧 t[70,90)/q[40,60)
    eng = LiftOverEngine([c], "target")
    r = eng.convert(BedRecord(1, "chrT", 70, 90, "", ""))
    assert r.status == "mapped"
    f = r.candidates[0].fragments[0]
    assert (f.dest_start, f.dest_end) == (40, 60)
    assert r.candidates[0].dest_strand == "-"

    back = LiftOverEngine([c], "query").convert(
        BedRecord(1, "chrQ", 40, 60, "", ""))
    bf = back.candidates[0].fragments[0]
    assert (bf.dest_start, bf.dest_end) == (70, 90)


def test_source_minus_dest_plus_intrablock_reverse():
    # 源（target）负链、目标（query）正链：块内反向映射
    c = Chain(
        "2", "chrT", 100, "-", 10, 30,
        "chrQ", 80, "+", 40, 60, 1000,
        [Block(10, 30, 40, 60)], 1,
    )
    # t 正向 [70,80) -> q [50,60)
    eng = LiftOverEngine([c], "target")
    f = eng.convert(BedRecord(1, "chrT", 70, 80, "", "")).candidates[0].fragments[0]
    assert (f.dest_start, f.dest_end) == (50, 60)
    assert f.dest_start < f.dest_end

    back = LiftOverEngine([c], "query").convert(
        BedRecord(1, "chrQ", 50, 60, "", "")).candidates[0].fragments[0]
    assert (back.dest_start, back.dest_end) == (70, 80)


def test_demo_chain_full_mapping(tmp_path):
    from app.chain import parse_chain
    import os
    demo = os.path.join(os.path.dirname(__file__), "..", "examples", "demo.chain")
    result = parse_chain(open(demo, encoding="utf-8").read())
    assert result.errors == []
    eng = LiftOverEngine(result.chains, "target")

    row = eng.convert(_rec(56, 176, "chr1"))
    assert row.status == "mapped"
    f = row.candidates[0]
    assert f.dest_chrom == "chrA"
    assert (f.fragments[0].dest_start, f.fragments[0].dest_end) == (100, 220)


def test_demo_chain_partial_split():
    from app.chain import parse_chain
    import os
    demo = os.path.join(os.path.dirname(__file__), "..", "examples", "demo.chain")
    chains = parse_chain(open(demo, encoding="utf-8").read()).chains
    eng = LiftOverEngine(chains, "target")
    row = eng.convert(_rec(50, 185, "chr1"))
    assert row.status == "partial"
    frags = row.candidates[0].fragments
    assert len(frags) == 2
    # 源片段 [50,56) 在缺口内 -> 被丢弃；第一块命中 [56,176)
    assert (frags[0].src_start, frags[0].src_end) == (56, 176)
    assert (frags[0].dest_start, frags[0].dest_end) == (100, 220)
    # gap [176,179) 丢弃；第二块命中 [179,185)
    assert (frags[1].src_start, frags[1].src_end) == (179, 185)
    assert (frags[1].dest_start, frags[1].dest_end) == (223, 229)


def test_demo_chain_unmapped_in_gap_and_outside():
    from app.chain import parse_chain
    import os
    demo = os.path.join(os.path.dirname(__file__), "..", "examples", "demo.chain")
    chains = parse_chain(open(demo, encoding="utf-8").read()).chains
    eng = LiftOverEngine(chains, "target")

    gap = eng.convert(_rec(16, 56, "chr1"))
    assert gap.status == "unmapped"
    assert "缺口" in gap.reason or "对齐范围" in gap.reason

    outside = eng.convert(_rec(500, 600, "chr1"))
    assert outside.status == "unmapped"

    nochrom = eng.convert(_rec(0, 10, chrom="chrX"))
    assert nochrom.status == "unmapped"
    assert "chrX" in nochrom.reason


def test_demo_chain_ambiguous_lists_all_candidates_without_picking():
    from app.chain import parse_chain
    import os
    demo = os.path.join(os.path.dirname(__file__), "..", "examples", "demo.chain")
    chains = parse_chain(open(demo, encoding="utf-8").read()).chains
    eng = LiftOverEngine(chains, "target")
    row = eng.convert(_rec(320, 360, "chr1"))
    assert row.status == "ambiguous"
    assert len(row.candidates) == 2
    dests = sorted(c.dest_chrom for c in row.candidates)
    assert dests == ["chrA", "chrB"]
    assert "不自动挑选" in row.reason
    # 负链候选坐标仍为左闭右开
    cand_b = next(c for c in row.candidates if c.dest_chrom == "chrB")
    f = cand_b.fragments[0]
    assert (f.dest_start, f.dest_end) == (240, 280)
    assert f.dest_start < f.dest_end
    assert cand_b.dest_strand == "-"
