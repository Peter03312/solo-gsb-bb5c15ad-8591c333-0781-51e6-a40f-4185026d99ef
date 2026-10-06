"""双段迁移结果的冲突审计：不同 BED 源行在最终坐标系上的重叠检测。

审计对象是某一次**已完成的双段转换**（``POST /api/liftover2`` 产生），
可选两种视图：

- ``raw``：原始全部路径——每个源行的每条可达路径都参与比较；
- ``adjudicated``：发起审计**那一刻**的人工裁决视图——已裁决的歧义行
  只取被选路径，未裁决的歧义行仍保留全部路径，其余行不变。

比较规则：

- 只比较**同一目标染色体**（最终坐标系）上、**不同源物理行**的零基
  左闭右开区间；同一源行的不同路径之间不比较；
- 每个（源行, 路径）组合先把自身的目标片段合并为互不相交的区间集，
  再对两个组合求交集——同一碱基不会被重复计数；
- 交集按区间逐条输出，按（染色体, 交集起点, 交集终点, 双方源行号,
  双方路径序号）稳定排序。

性能：按目标染色体分组后在组内做扫描线（起点排序 + 终点小根堆活动集），
只与真正重叠的活动区间求交，不做全量两两比较；复杂度
O(n log n + k)，k 为实际重叠的区间对数。
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

from .engine2 import RowResult2

# 允许的审计视图：raw=原始全部路径；adjudicated=发起时的人工裁决视图
AUDIT_VIEWS = ("raw", "adjudicated")


@dataclass(frozen=True)
class AuditPath:
    """参与审计的一个（源行, 路径）组合：目标片段已合并为互不相交的区间。"""

    line_no: int                 # BED 源物理行号
    path_index: int              # 1 基路径序号（与双段结果 TSV 一致）
    chain1_id: str
    chain2_id: str
    dest_chrom: str              # 最终坐标系染色体
    intervals: tuple[tuple[int, int], ...]   # 合并后的目标区间（有序、互不相交）


@dataclass(frozen=True)
class Conflict:
    """两个（源行, 路径）组合在最终坐标系上的一段交集（零基左闭右开）。"""

    chrom: str
    overlap_start: int
    overlap_end: int
    overlap_bases: int
    line_a: int                  # 源行号较小的一侧为 A
    path_a: int
    chain1_a: str
    chain2_a: str
    line_b: int
    path_b: int
    chain1_b: str
    chain2_b: str


def merge_intervals(
    intervals: Iterable[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    """合并为互不相交的有序区间（相邻相接的区间也并为一个）。

    与 ``engine2._union_cov`` 的并集口径一致：``s <= 当前终点`` 即合并。
    """
    merged: list[list[int]] = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1]:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return tuple((s, e) for s, e in merged)


def select_audit_paths(
    rows: Sequence[RowResult2],
    view: str,
    decisions: Mapping[int, int],
) -> list[AuditPath]:
    """按视图挑选参与审计的（源行, 路径）组合，并合并各自的目标片段。

    - ``raw``：全部路径参与；
    - ``adjudicated``：已裁决的歧义行只取被选路径，未裁决的歧义行保留
      全部路径（decisions 只可能落在歧义行上；此处对非法键容错，按未
      裁决处理，与裁决版 TSV 的口径一致）。

    无映射行没有路径，天然不参与比较。
    """
    combos: list[AuditPath] = []
    for row in rows:
        paths = row.candidates
        if not paths:
            continue
        if view == "adjudicated" and row.status == "ambiguous":
            picked = decisions.get(row.line_no)
            if picked is not None and 1 <= picked <= len(paths):
                paths = [paths[picked - 1]]
        for path in paths:
            merged = merge_intervals(
                (f.dest_start, f.dest_end) for f in path.fragments
            )
            combos.append(
                AuditPath(
                    line_no=row.line_no,
                    path_index=path.path_index,
                    chain1_id=path.chain1_id,
                    chain2_id=path.chain2_id,
                    dest_chrom=path.dest_chrom,
                    intervals=merged,
                )
            )
    return combos


def find_conflicts(combos: Sequence[AuditPath]) -> list[Conflict]:
    """检测不同源行的组合在最终坐标系上的交集（扫描线，非全量两两比较）。

    按目标染色体分组；组内把各组合已合并的区间按起点排序扫描，只对
    终点越过当前起点的活动区间求交。同一组合内区间互不相交，因此同一
    对组合产生的各段交集互不重叠，同一碱基不会被重复计数。
    """
    by_chrom: dict[str, list[int]] = {}
    for ci, combo in enumerate(combos):
        by_chrom.setdefault(combo.dest_chrom, []).append(ci)

    conflicts: list[Conflict] = []
    for chrom, idxs in by_chrom.items():
        # （起点, 终点, 组合序号）按起点排序；同起点时按终点、序号定序，
        # 保证扫描顺序确定
        entries = sorted(
            (s, e, ci)
            for ci in idxs
            for (s, e) in combos[ci].intervals
        )
        active: list[tuple[int, int, int]] = []  # （终点, 起点, 组合序号）小根堆
        for s, e, ci in entries:
            while active and active[0][0] <= s:
                heapq.heappop(active)
            for a_end, a_start, a_ci in active:
                other = combos[a_ci]
                cur = combos[ci]
                if other.line_no == cur.line_no:
                    continue  # 只审计不同源物理行之间的重叠
                ov_s = max(a_start, s)
                ov_e = min(a_end, e)
                if ov_s >= ov_e:
                    continue
                # 每对组合只输出一次：源行号（再路径序号）较小的一侧为 A
                lo, hi = (other, cur) if (
                    (other.line_no, other.path_index)
                    <= (cur.line_no, cur.path_index)
                ) else (cur, other)
                conflicts.append(
                    Conflict(
                        chrom=chrom,
                        overlap_start=ov_s,
                        overlap_end=ov_e,
                        overlap_bases=ov_e - ov_s,
                        line_a=lo.line_no,
                        path_a=lo.path_index,
                        chain1_a=lo.chain1_id,
                        chain2_a=lo.chain2_id,
                        line_b=hi.line_no,
                        path_b=hi.path_index,
                        chain1_b=hi.chain1_id,
                        chain2_b=hi.chain2_id,
                    )
                )
            heapq.heappush(active, (e, s, ci))

    conflicts.sort(
        key=lambda c: (
            c.chrom, c.overlap_start, c.overlap_end,
            c.line_a, c.path_a, c.line_b, c.path_b,
        )
    )
    return conflicts


def summarize(
    combos: Sequence[AuditPath],
    conflicts: Sequence[Conflict],
) -> dict[str, int]:
    """审计汇总：参与规模与冲突规模（无冲突时冲突三项全为 0）。"""
    pairs = {(c.line_a, c.path_a, c.line_b, c.path_b) for c in conflicts}
    return {
        "source_lines": len({c.line_no for c in combos}),
        "paths": len(combos),
        "conflict_pairs": len(pairs),
        "conflict_intervals": len(conflicts),
        "overlap_bases": sum(c.overlap_bases for c in conflicts),
    }
