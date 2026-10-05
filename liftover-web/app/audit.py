"""双段迁移结果的源行冲突审计（最终坐标系重叠检测）。

对一次**已完成的双段迁移任务**（``POST /api/liftover2`` 的结果），审计
不同 BED 源物理行在最终坐标系（第二段目标侧）上的区间冲突：

- 审计单元是（源行 × 路径）组合，可选两种视图：

  - ``original``：原始全部路径（口径与原始下载一致）；
  - ``adjudicated``：发起时的人工裁决视图——已裁决的歧义行只保留被选
    路径，未裁决的歧义行仍保留全部路径（口径与裁决版 TSV 一致）。
    裁决状态在发起时快照，之后的裁决改动不影响已生成的审计结果，
    原始下载字节也始终不变。

- 只比较**同一目标染色体**上、**不同源物理行**的组合；同一源行的多条
  路径是互斥候选，彼此不比较；
- 每个组合先把自己的目标片段按染色体**合并**（重叠或相接的区间并为
  一个），再对两个组合求交集——同一碱基不会被重复计数；
- 区间一律零基左闭右开；仅相接（``end == start``）不算重叠；
- 不做全量两两比较：按目标染色体分组后用扫描线，只对扫描窗口内实际
  相交的区间对求交，代价随输出规模增长而不是随组合数平方增长。

结果按（染色体、交集起、交集止、双方源行号、双方路径序号）稳定排序；
每条明细列出双方的两段 chain ID、交集坐标与重叠碱基数。无冲突时输出
空明细与零汇总（不是错误）。
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from typing import Mapping

from .engine2 import RowResult2

VIEW_ORIGINAL = "original"
VIEW_ADJUDICATED = "adjudicated"
AUDIT_VIEWS = (VIEW_ORIGINAL, VIEW_ADJUDICATED)


@dataclass(frozen=True)
class AuditUnit:
    """参与审计的（源行 × 路径）组合。

    ``merged`` 按目标染色体给出该组合合并后的目标片段（零基左闭右开、
    互不重叠也不相接、按起点升序）。
    """

    line_no: int                 # BED 源物理行号（1 基）
    path_index: int              # 1 基路径序号（与原结果一致）
    chain1_id: str
    chain2_id: str
    merged: dict[str, tuple[tuple[int, int], ...]]


@dataclass(frozen=True)
class Conflict:
    """一条冲突明细：两个组合在同一目标染色体上的一段交集。"""

    chrom: str
    start: int                   # 交集区间 [start, end)，零基左闭右开
    end: int
    a_line: int                  # 键序（源行号, 路径序号）较小的一方
    a_path: int
    a_chain1: str
    a_chain2: str
    b_line: int
    b_path: int
    b_chain1: str
    b_chain2: str

    @property
    def bases(self) -> int:
        """重叠碱基数（双方均已先合并，交集内同一碱基只计一次）。"""
        return self.end - self.start


def merge_intervals(ivs) -> tuple[tuple[int, int], ...]:
    """合并重叠或相接的 [s, e) 区间，返回升序、互不相接的区间元组。"""
    out: list[list[int]] = []
    for s, e in sorted(ivs):
        if out and s <= out[-1][1]:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return tuple((s, e) for s, e in out)


def build_audit_units(
    rows: list[RowResult2],
    view: str,
    decisions: Mapping[int, int],
) -> list[AuditUnit]:
    """按视图从双段结果中挑选参与审计的（源行 × 路径）组合。

    - ``original``：所有有路径的行的全部路径；
    - ``adjudicated``：歧义行若已裁决（``decisions`` 给出有效路径序号）
      只保留被选路径，未裁决的歧义行仍保留全部路径；其余行不变；
    - 无映射行没有最终坐标，不参与审计。
    """
    units: list[AuditUnit] = []
    for row in rows:
        if not row.candidates:
            continue
        paths = row.candidates
        if view == VIEW_ADJUDICATED and row.status == "ambiguous":
            picked = decisions.get(row.line_no)
            if picked is not None and 1 <= picked <= len(row.candidates):
                paths = [row.candidates[picked - 1]]
        for path in paths:
            by_chrom: dict[str, list[tuple[int, int]]] = {}
            for frag in path.fragments:
                by_chrom.setdefault(path.dest_chrom, []).append(
                    (frag.dest_start, frag.dest_end)
                )
            merged = {chrom: merge_intervals(ivs) for chrom, ivs in by_chrom.items()}
            units.append(
                AuditUnit(
                    line_no=row.line_no,
                    path_index=path.path_index,
                    chain1_id=path.chain1_id,
                    chain2_id=path.chain2_id,
                    merged=merged,
                )
            )
    return units


def _interval_overlap(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int] | None:
    """两个零基左闭右开区间的交集；无交集（含仅相接）返回 None。"""
    lo = a[0] if a[0] > b[0] else b[0]
    hi = a[1] if a[1] < b[1] else b[1]
    return (lo, hi) if lo < hi else None


def _make_conflict(
    chrom: str, lo: int, hi: int, u1: AuditUnit, u2: AuditUnit
) -> Conflict:
    """构造一条冲突明细；双方按（源行号, 路径序号）键序规范为 a/b。"""
    if (u2.line_no, u2.path_index) < (u1.line_no, u1.path_index):
        u1, u2 = u2, u1
    return Conflict(
        chrom=chrom, start=lo, end=hi,
        a_line=u1.line_no, a_path=u1.path_index,
        a_chain1=u1.chain1_id, a_chain2=u1.chain2_id,
        b_line=u2.line_no, b_path=u2.path_index,
        b_chain1=u2.chain1_id, b_chain2=u2.chain2_id,
    )


def find_conflicts(units: list[AuditUnit]) -> list[Conflict]:
    """找出所有冲突明细（扫描线算法，不做全量两两比较）。

    按目标染色体分组后，各区间的起点扫描线只让当前区间与"尚未结束"
    的区间求交；同一源物理行的组合（互斥路径）直接跳过。每个组合的
    区间已合并，故同一对组合在同一染色体上的各交集互不相交，同一
    碱基不会被重复计数。
    """
    by_chrom: dict[str, list[tuple[int, int, int]]] = {}
    for ui, unit in enumerate(units):
        for chrom, ivs in unit.merged.items():
            for s, e in ivs:
                by_chrom.setdefault(chrom, []).append((s, e, ui))

    conflicts: list[Conflict] = []
    for chrom in sorted(by_chrom):
        entries = sorted(by_chrom[chrom])
        active: list[tuple[int, int, int]] = []   # (end, start, unit_index)
        for s, e, ui in entries:
            # 零基左闭右开：end <= s 的区间与当前区间不重叠，移出窗口
            active = [a for a in active if a[0] > s]
            for a_end, a_start, a_ui in active:
                ua, ub = units[a_ui], units[ui]
                if ua.line_no == ub.line_no:
                    continue   # 同一源行的多条路径互斥，不互相比较
                ov = _interval_overlap((a_start, a_end), (s, e))
                if ov is not None:
                    conflicts.append(_make_conflict(chrom, ov[0], ov[1], ua, ub))
            active.append((e, s, ui))

    conflicts.sort(key=lambda c: (
        c.chrom, c.start, c.end, c.a_line, c.a_path, c.b_line, c.b_path,
    ))
    return conflicts


def summarize(units: list[AuditUnit], conflicts: list[Conflict]) -> dict[str, int]:
    """审计汇总：参与规模与冲突规模（无冲突时冲突项全为 0）。"""
    pairs = {(c.a_line, c.a_path, c.b_line, c.b_path) for c in conflicts}
    return {
        "source_rows": len({u.line_no for u in units}),
        "units": len(units),
        "conflict_pairs": len(pairs),
        "conflict_regions": len(conflicts),
        "conflict_bases": sum(c.bases for c in conflicts),
    }


def conflict_to_dict(c: Conflict) -> dict:
    """冲突明细的 JSON 表示（与 TSV 列一一对应）。"""
    return {
        "chrom": c.chrom,
        "start": c.start,
        "end": c.end,
        "bases": c.bases,
        "a_source_line": c.a_line,
        "a_path_index": c.a_path,
        "a_chain1_id": c.a_chain1,
        "a_chain2_id": c.a_chain2,
        "b_source_line": c.b_line,
        "b_path_index": c.b_path,
        "b_chain1_id": c.b_chain1,
        "b_chain2_id": c.b_chain2,
    }


AUDIT_TSV_HEADER = [
    "chrom",
    "overlap_start",
    "overlap_end",
    "overlap_bases",
    "a_source_line",
    "a_path_index",
    "a_chain1_id",
    "a_chain2_id",
    "b_source_line",
    "b_path_index",
    "b_chain1_id",
    "b_chain2_id",
]


def build_audit_tsv(conflicts: list[Conflict]) -> str:
    """审计明细 TSV（无冲突时只有表头）。"""
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n",
                        quoting=csv.QUOTE_MINIMAL)
    writer.writerow(AUDIT_TSV_HEADER)
    for c in conflicts:
        writer.writerow([
            c.chrom, c.start, c.end, c.bases,
            c.a_line, c.a_path, c.a_chain1, c.a_chain2,
            c.b_line, c.b_path, c.b_chain1, c.b_chain2,
        ])
    return buf.getvalue()
