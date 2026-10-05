"""双段（两份 chain 依次）迁移引擎测试。"""

from app.bed import BedRecord
from app.chain import Block, Chain
from app.engine2 import TwoSegmentEngine


def _rec(start: int, end: int, chrom: str = "chr1", line_no: int = 1) -> BedRecord:
    return BedRecord(line_no, chrom, start, end, f"{chrom}\t{start}\t{end}", "")


def _seg1_chain() -> Chain:
    # chr1(target) -> chrA(query)，正链；块 [0,120) 与 [150,200)，缺口 [120,150)
    return Chain(
        chain_id="1", t_name="chr1", t_size=2000, t_strand="+",
        t_start=0, t_end=200,
        q_name="chrA", q_size=500, q_strand="+",
        q_start=0, q_end=200, score=1000,
        blocks=[Block(0, 120, 0, 120), Block(150, 200, 150, 200)],
        header_line=1,
    )


def _seg2_chain() -> Chain:
    # chrA(target) -> chrZ(query)，正链；块 [0,100) 与 [120,200)，缺口 [100,120)
    return Chain(
        chain_id="10", t_name="chrA", t_size=500, t_strand="+",
        t_start=0, t_end=200,
        q_name="chrZ", q_size=1000, q_strand="+",
        q_start=0, q_end=200, score=500,
        blocks=[Block(0, 100, 0, 100), Block(120, 200, 120, 200)],
        header_line=1,
    )


def _engine(chains1=None, chains2=None, side1="target", side2="target"):
    return TwoSegmentEngine(
        chains1 if chains1 is not None else [_seg1_chain()], side1,
        chains2 if chains2 is not None else [_seg2_chain()], side2,
    )


# ---------- 基本映射 / 两段缺口 ----------

def test_full_two_segment_mapping():
    row = _engine().convert(_rec(0, 100))
    assert row.status == "mapped"
    assert row.failed_stage == ""
    assert len(row.candidates) == 1
    path = row.candidates[0]
    assert (path.chain1_id, path.chain2_id) == ("1", "10")
    assert path.path_index == 1
    assert (path.mid_chrom, path.dest_chrom) == ("chrA", "chrZ")
    assert (path.seg1_covered, path.seg2_covered, path.total) == (100, 100, 100)
    assert path.uncovered_src == []
    frag = path.fragments[0]
    assert (frag.src_start, frag.src_end) == (0, 100)
    assert (frag.mid_start, frag.mid_end) == (0, 100)
    assert (frag.dest_start, frag.dest_end) == (0, 100)


def test_segment2_gap_makes_partial_and_keeps_uncovered():
    # 第一段完整 [0,120)->chrA[0,120)；第二段 chrA[100,120) 落在缺口
    row = _engine().convert(_rec(0, 120))
    assert row.status == "partial"
    assert "完全映射" not in row.reason
    path = row.candidates[0]
    assert (path.seg1_covered, path.seg2_covered, path.total) == (120, 100, 120)
    assert path.uncovered_src == [(100, 120)]
    assert "[100,120)" in row.reason
    frag = path.fragments[0]
    assert (frag.src_start, frag.src_end) == (0, 100)
    assert (frag.dest_start, frag.dest_end) == (0, 100)


def test_segment1_gap_plus_segment2_gap_partial():
    # 源 [100,160)：第一段命中 [100,120)+[150,160)（缺口 120-150）；
    # 中间 chrA[100,120) 又落在第二段缺口，只有 chrA[150,160) 完成两段
    row = _engine().convert(_rec(100, 160))
    assert row.status == "partial"
    path = row.candidates[0]
    assert (path.seg1_covered, path.seg2_covered, path.total) == (30, 10, 60)
    assert path.uncovered_src == [(100, 150)]
    assert len(path.fragments) == 1
    frag = path.fragments[0]
    assert (frag.src_start, frag.src_end) == (150, 160)
    assert (frag.mid_start, frag.mid_end) == (150, 160)
    assert (frag.dest_start, frag.dest_end) == (150, 160)


# ---------- 失败阶段 ----------

def test_unmapped_segment1_missing_chrom():
    row = _engine().convert(_rec(0, 10, chrom="chrX"))
    assert row.status == "unmapped"
    assert row.failed_stage == "segment1"
    assert "chain 文件 1" in row.reason and "chrX" in row.reason
    assert row.candidates == []
    assert (row.seg1_covered, row.seg2_covered) == (0, 0)


def test_unmapped_segment1_in_gap():
    row = _engine().convert(_rec(120, 150))
    assert row.status == "unmapped"
    assert row.failed_stage == "segment1"
    assert "chain 文件 1" in row.reason


def test_unmapped_segment2_missing_chrom():
    c1 = Chain(
        "3", "chr1", 2000, "+", 500, 600,
        "chrM", 500, "+", 0, 100, 1000,
        [Block(500, 600, 0, 100)], 1,
    )
    row = _engine(chains1=[c1]).convert(_rec(500, 600))
    assert row.status == "unmapped"
    assert row.failed_stage == "segment2"
    assert "chain 文件 2" in row.reason and "chrM" in row.reason
    # 第一段实际映射了 100 bp，但没有任何碱基完成两段
    assert (row.seg1_covered, row.seg2_covered) == (100, 0)


def test_unmapped_segment2_outside_alignment():
    c1 = Chain(
        "4", "chr1", 2000, "+", 700, 800,
        "chrA", 500, "+", 300, 400, 1000,
        [Block(700, 800, 300, 400)], 1,
    )
    # chrA 在 chain 文件 2 中存在，但 [300,400) 不在任何对齐块内
    row = _engine(chains1=[c1]).convert(_rec(700, 800))
    assert row.status == "unmapped"
    assert row.failed_stage == "segment2"
    assert "chain 文件 2" in row.reason
    assert (row.seg1_covered, row.seg2_covered) == (100, 0)


# ---------- 回溯方向 ----------

def test_traceback_through_minus_strand_segment1():
    # 第一段目标侧为负链：chr1[0,100) <-> chrA 存储帧 [400,500)（正向帧 [0,100)），
    # 块内反向。第二段只覆盖 chrA[30,60)。
    c1 = Chain(
        "1", "chr1", 1000, "+", 0, 100,
        "chrA", 500, "-", 400, 500, 1000,
        [Block(0, 100, 400, 500)], 1,
    )
    c2 = Chain(
        "10", "chrA", 500, "+", 30, 60,
        "chrZ", 1000, "+", 300, 330, 500,
        [Block(30, 60, 300, 330)], 1,
    )
    row = _engine(chains1=[c1], chains2=[c2]).convert(_rec(0, 100))
    assert row.status == "partial"
    path = row.candidates[0]
    assert path.seg2_covered == 30
    frag = path.fragments[0]
    # 中间 [30,60) 经反向回溯对应源 [40,70)，而不是 [30,60)
    assert (frag.src_start, frag.src_end) == (40, 70)
    assert (frag.mid_start, frag.mid_end) == (30, 60)
    assert (frag.dest_start, frag.dest_end) == (300, 330)
    assert path.uncovered_src == [(0, 40), (70, 100)]


def test_traceback_through_offset_same_orientation():
    # 第一段平移：chr1[100,200) -> chrA[50,150)；第二段覆盖 chrA[60,90)
    c1 = Chain(
        "1", "chr1", 1000, "+", 100, 200,
        "chrA", 500, "+", 50, 150, 1000,
        [Block(100, 200, 50, 150)], 1,
    )
    c2 = Chain(
        "10", "chrA", 500, "+", 60, 90,
        "chrZ", 1000, "+", 600, 630, 500,
        [Block(60, 90, 600, 630)], 1,
    )
    row = _engine(chains1=[c1], chains2=[c2]).convert(_rec(100, 200))
    assert row.status == "partial"
    frag = row.candidates[0].fragments[0]
    assert (frag.src_start, frag.src_end) == (110, 140)
    assert (frag.mid_start, frag.mid_end) == (60, 90)
    assert (frag.dest_start, frag.dest_end) == (600, 630)


def test_fragments_sorted_by_source_after_minus_strand_traceback():
    # 第一段反向；第二段两个块命中中间 [10,20) 与 [40,50)，
    # 回溯后源区间顺序翻转，输出必须按源坐标升序。
    c1 = Chain(
        "1", "chr1", 1000, "+", 0, 100,
        "chrA", 500, "-", 400, 500, 1000,
        [Block(0, 100, 400, 500)], 1,
    )
    c2 = Chain(
        "10", "chrA", 500, "+", 10, 50,
        "chrZ", 1000, "+", 0, 50, 500,
        [Block(10, 20, 0, 10), Block(40, 50, 40, 50)], 1,
    )
    row = _engine(chains1=[c1], chains2=[c2]).convert(_rec(0, 100))
    frags = row.candidates[0].fragments
    assert [(f.src_start, f.src_end) for f in frags] == [(50, 60), (80, 90)]
    assert [(f.mid_start, f.mid_end) for f in frags] == [(40, 50), (10, 20)]


# ---------- 多路径歧义 ----------

def _ambiguous_engine() -> TwoSegmentEngine:
    # 第一段 chr1[300,400) -> chrB[0,100)；第二段 chrB[0,100) 同时到 chrW 与 chrV
    c1 = Chain(
        "2", "chr1", 2000, "+", 300, 400,
        "chrB", 400, "+", 0, 100, 1000,
        [Block(300, 400, 0, 100)], 1,
    )
    c2a = Chain(
        "11", "chrB", 400, "+", 0, 100,
        "chrW", 800, "+", 0, 100, 500,
        [Block(0, 100, 0, 100)], 1,
    )
    c2b = Chain(
        "12", "chrB", 400, "+", 0, 100,
        "chrV", 900, "+", 50, 150, 500,
        [Block(0, 100, 50, 150)], 1,
    )
    return TwoSegmentEngine([c1], "target", [c2a, c2b], "target")


def test_multiple_reachable_paths_all_kept_as_ambiguous():
    row = _ambiguous_engine().convert(_rec(300, 400))
    assert row.status == "ambiguous"
    assert len(row.candidates) == 2
    assert "不自动选择" in row.reason
    dests = {p.dest_chrom for p in row.candidates}
    assert dests == {"chrW", "chrV"}
    assert [p.path_index for p in row.candidates] == [1, 2]
    for p in row.candidates:
        assert (p.seg1_covered, p.seg2_covered, p.total) == (100, 100, 100)
        assert p.uncovered_src == []


def test_paths_via_two_segment1_chains_all_kept():
    c1a = Chain(
        "2", "chr1", 2000, "+", 300, 400,
        "chrB", 400, "+", 0, 100, 1000,
        [Block(300, 400, 0, 100)], 1,
    )
    c1b = Chain(
        "3", "chr1", 2000, "+", 300, 400,
        "chrC", 400, "+", 0, 100, 1000,
        [Block(300, 400, 0, 100)], 1,
    )
    c2a = Chain(
        "11", "chrB", 400, "+", 0, 100,
        "chrW", 800, "+", 0, 100, 500,
        [Block(0, 100, 0, 100)], 1,
    )
    c2b = Chain(
        "12", "chrC", 400, "+", 0, 100,
        "chrV", 900, "+", 0, 100, 500,
        [Block(0, 100, 0, 100)], 1,
    )
    eng = TwoSegmentEngine([c1a, c1b], "target", [c2a, c2b], "target")
    row = eng.convert(_rec(300, 400))
    assert row.status == "ambiguous"
    assert {(p.chain1_id, p.chain2_id) for p in row.candidates} == {
        ("2", "11"), ("3", "12"),
    }


def test_dead_end_segment1_hit_noted_but_not_a_path():
    # 第一段两个命中：chrB（可达 chrW）与 chrM（第二段不存在该染色体）。
    # 只有一条可达路径 -> 不判歧义；死路只在 reason 中说明。
    c1a = Chain(
        "2", "chr1", 2000, "+", 300, 400,
        "chrB", 400, "+", 0, 100, 1000,
        [Block(300, 400, 0, 100)], 1,
    )
    c1b = Chain(
        "3", "chr1", 2000, "+", 300, 400,
        "chrM", 500, "+", 0, 100, 1000,
        [Block(300, 400, 0, 100)], 1,
    )
    c2 = Chain(
        "11", "chrB", 400, "+", 0, 100,
        "chrW", 800, "+", 0, 100, 500,
        [Block(0, 100, 0, 100)], 1,
    )
    eng = TwoSegmentEngine([c1a, c1b], "target", [c2], "target")
    row = eng.convert(_rec(300, 400))
    assert row.status == "mapped"
    assert len(row.candidates) == 1
    assert "未形成可达路径" in row.reason
    assert "chain 3" in row.reason


def test_same_chain2_hit_by_multiple_mid_fragments_forms_one_path():
    # 第一段两个片段都落在同一条第二段 chain 上 -> 同一路径两个片段
    row = _engine().convert(_rec(50, 160))
    assert row.status == "partial"
    assert len(row.candidates) == 1
    path = row.candidates[0]
    # 第一段片段 chr1[50,120)+[150,160)；第二段只覆盖 chrA[50,100) 与 [150,160)
    assert [(f.src_start, f.src_end) for f in path.fragments] == [(50, 100), (150, 160)]
    assert path.seg1_covered == 80
    assert path.seg2_covered == 60
    assert path.uncovered_src == [(100, 150)]


# ---------- 侧别组合 ----------

def test_query_side_segment1():
    # side1=query：BED 在 chain 文件 1 的 query 侧（chrA），第一段产物在
    # chrZ 上；chain 文件 2 的源侧没有 chrZ -> 第二段无映射
    c1 = Chain(
        "1", "chrZ", 1000, "+", 0, 100,
        "chrA", 500, "+", 0, 100, 1000,
        [Block(0, 100, 0, 100)], 1,
    )
    eng = TwoSegmentEngine([c1], "query", [_seg2_chain()], "target")
    row = eng.convert(_rec(0, 100, chrom="chrA"))
    assert row.status == "unmapped"
    assert row.failed_stage == "segment2"
    assert "chrZ" in row.reason


def test_query_side_roundtrip_consistent_with_single_segment():
    # 双段 chain2 是 chain1 的逆（同一份文件、侧别互换）时应回到原坐标
    c = Chain(
        "1", "chr1", 1000, "+", 0, 100,
        "chrA", 500, "+", 50, 150, 1000,
        [Block(0, 100, 50, 150)], 1,
    )
    eng = TwoSegmentEngine([c], "target", [c], "query")
    row = eng.convert(_rec(0, 100))
    assert row.status == "mapped"
    frag = row.candidates[0].fragments[0]
    assert (frag.src_start, frag.src_end) == (0, 100)
    assert (frag.mid_start, frag.mid_end) == (50, 150)
    assert (frag.dest_start, frag.dest_end) == (0, 100)
