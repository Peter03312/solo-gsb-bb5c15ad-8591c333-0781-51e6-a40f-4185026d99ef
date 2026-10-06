"""冲突审计核心逻辑单元测试：区间合并、视图选路、扫描线求交。"""

from app.audit import (
    find_conflicts,
    merge_intervals,
    select_audit_paths,
    summarize,
)
from app.engine2 import PathCandidate, PathFragment, RowResult2


def _path(path_index, dest_chrom, dest_intervals, chain1="1", chain2="10"):
    return PathCandidate(
        path_index=path_index,
        chain1_id=chain1,
        chain2_id=chain2,
        mid_chrom="chrA",
        mid_strand="+",
        dest_chrom=dest_chrom,
        dest_strand="+",
        fragments=[
            PathFragment(
                src_start=0, src_end=0, mid_start=0, mid_end=0,
                dest_start=s, dest_end=e,
            )
            for s, e in dest_intervals
        ],
        seg1_covered=0,
        seg2_covered=0,
        total=100,
    )


def _row(line_no, paths, status=None):
    if status is None:
        status = "ambiguous" if len(paths) > 1 else "mapped"
    return RowResult2(
        line_no=line_no,
        src_chrom="chr1",
        src_start=0,
        src_end=100,
        raw="",
        status=status,
        candidates=paths,
    )


# ---------- merge_intervals ----------

def test_merge_intervals_unsorted_and_overlapping():
    assert merge_intervals([(50, 80), (0, 30), (20, 60)]) == ((0, 80),)


def test_merge_intervals_touching_merged():
    # 相接（左闭右开的端点相接）也并为一个区间
    assert merge_intervals([(0, 50), (50, 100)]) == ((0, 100),)


def test_merge_intervals_nested_and_disjoint():
    assert merge_intervals([(0, 100), (10, 40), (200, 300)]) == ((0, 100), (200, 300))


# ---------- select_audit_paths ----------

def test_select_raw_keeps_all_paths_and_merges_fragments():
    rows = [
        _row(1, [_path(1, "chrZ", [(0, 50), (50, 100)]),
                 _path(2, "chrW", [(0, 100)])]),
        _row(2, [_path(1, "chrZ", [(40, 140)])]),
        _row(3, [], status="unmapped"),   # 无映射行不参与
    ]
    combos = select_audit_paths(rows, "raw", decisions={})
    assert [(c.line_no, c.path_index) for c in combos] == [(1, 1), (1, 2), (2, 1)]
    by_key = {(c.line_no, c.path_index): c for c in combos}
    # 自身目标片段先合并：相接片段 [0,50)+[50,100) 并成 [0,100)
    assert by_key[(1, 1)].intervals == ((0, 100),)
    assert by_key[(1, 1)].dest_chrom == "chrZ"
    assert by_key[(1, 2)].chain1_id == "1"


def test_select_adjudicated_narrows_only_decided_ambiguous_rows():
    rows = [
        _row(1, [_path(1, "chrZ", [(0, 100)]),
                 _path(2, "chrW", [(0, 100)])]),          # 歧义，已裁决路径 2
        _row(2, [_path(1, "chrZ", [(0, 100)]),
                 _path(2, "chrW", [(0, 100)])]),          # 歧义，未裁决
        _row(3, [_path(1, "chrV", [(0, 100)])], status="mapped"),
    ]
    combos = select_audit_paths(rows, "adjudicated", decisions={1: 2})
    got = [(c.line_no, c.path_index, c.dest_chrom) for c in combos]
    # 行1 只保留被选路径 2；行2 未裁决保留全部路径；行3 非歧义不受影响
    assert got == [(1, 2, "chrW"), (2, 1, "chrZ"), (2, 2, "chrW"), (3, 1, "chrV")]


def test_select_adjudicated_tolerates_bogus_decision():
    rows = [_row(1, [_path(1, "chrZ", [(0, 100)]), _path(2, "chrW", [(0, 100)])])]
    # 越界的选择按未裁决处理（接口层本就不会产生，此处容错）
    combos = select_audit_paths(rows, "adjudicated", decisions={1: 99})
    assert [(c.line_no, c.path_index) for c in combos] == [(1, 1), (1, 2)]


# ---------- find_conflicts ----------

def test_same_line_paths_never_compared():
    combos = select_audit_paths(
        [_row(1, [_path(1, "chrZ", [(0, 100)]), _path(2, "chrZ", [(10, 60)])])],
        "raw", {},
    )
    assert find_conflicts(combos) == []


def test_touching_intervals_do_not_conflict():
    combos = select_audit_paths(
        [_row(1, [_path(1, "chrZ", [(0, 50)])]),
         _row(2, [_path(1, "chrZ", [(50, 100)])])],
        "raw", {},
    )
    assert find_conflicts(combos) == []


def test_different_chroms_not_compared():
    combos = select_audit_paths(
        [_row(1, [_path(1, "chrZ", [(0, 100)])]),
         _row(2, [_path(1, "chrW", [(0, 100)])])],
        "raw", {},
    )
    assert find_conflicts(combos) == []


def test_nested_interval_conflict_canonical_order():
    combos = select_audit_paths(
        [_row(5, [_path(2, "chrZ", [(0, 100)], chain1="7", chain2="8")]),
         _row(2, [_path(1, "chrZ", [(10, 40)], chain1="3", chain2="4")])],
        "raw", {},
    )
    (c,) = find_conflicts(combos)
    # 源行号较小的一侧为 A，与组合传入顺序无关
    assert (c.line_a, c.path_a, c.chain1_a, c.chain2_a) == (2, 1, "3", "4")
    assert (c.line_b, c.path_b, c.chain1_b, c.chain2_b) == (5, 2, "7", "8")
    assert (c.overlap_start, c.overlap_end, c.overlap_bases) == (10, 40, 30)


def test_multi_interval_pair_intersection_no_double_count():
    # A = [0,10)+[100,110)，B = [5,105)：两段交集分别输出，互不重叠
    combos = select_audit_paths(
        [_row(1, [_path(1, "chrZ", [(0, 10), (100, 110)])]),
         _row(2, [_path(1, "chrZ", [(5, 105)])])],
        "raw", {},
    )
    conflicts = find_conflicts(combos)
    assert [(c.overlap_start, c.overlap_end, c.overlap_bases) for c in conflicts] == [
        (5, 10, 5),
        (100, 105, 5),
    ]
    summary = summarize(combos, conflicts)
    assert summary["conflict_pairs"] == 1
    assert summary["conflict_intervals"] == 2
    assert summary["overlap_bases"] == 10


def test_output_sorted_by_chrom_start_end_lines_paths():
    combos = select_audit_paths(
        [_row(1, [_path(1, "chrZ", [(0, 100)]), _path(2, "chrW", [(0, 100)])]),
         _row(2, [_path(1, "chrZ", [(20, 120)])]),
         _row(3, [_path(1, "chrW", [(50, 150)])]),
         _row(4, [_path(1, "chrZ", [(10, 30)])])],
        "raw", {},
    )
    conflicts = find_conflicts(combos)
    keys = [
        (c.chrom, c.overlap_start, c.overlap_end, c.line_a, c.path_a, c.line_b, c.path_b)
        for c in conflicts
    ]
    assert keys == sorted(keys)
    # chrW 只有一对；chrZ 上 (1,1)×(2,1)、(1,1)×(4,1)、(2,1)×(4,1) 三对
    assert keys[0][0] == "chrW"
    assert [k for k in keys if k[0] == "chrZ"] == [
        ("chrZ", 10, 30, 1, 1, 4, 1),
        ("chrZ", 20, 30, 2, 1, 4, 1),
        ("chrZ", 20, 100, 1, 1, 2, 1),
    ]


def test_summarize_zero_when_no_conflict():
    combos = select_audit_paths([_row(1, [_path(1, "chrZ", [(0, 100)])])], "raw", {})
    conflicts = find_conflicts(combos)
    assert conflicts == []
    assert summarize(combos, conflicts) == {
        "source_lines": 1,
        "paths": 1,
        "conflict_pairs": 0,
        "conflict_intervals": 0,
        "overlap_bases": 0,
    }
