"""双段（两份 chain 依次）区间迁移引擎。

第一段把 BED 源区间迁到中间坐标系，第二段只处理第一段**实际映射到**
的片段，再把它们迁到最终坐标系；每个最终片段都按映射方向（第一段
片段的块内仿射逆映射）精确回溯到原 BED 源坐标，坐标一律零基左闭右开。

规则（与单段一致的严格语义）：

- 任一段的缺口/未覆盖都体现为源区间的未覆盖信息，残余片段绝不标成
  完全映射；
- 多条可达路径（第一段 chain × 第二段 chain 的组合）全部保留为歧义，
  不自动选择、不合并；
- 没有任何最终映射时记录失败阶段（``segment1`` / ``segment2``），
  第二段源染色体不存在同样按第二段无映射处理。
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .bed import BedRecord
from .chain import Chain
from .engine import LiftOverEngine


@dataclass(frozen=True)
class PathFragment:
    src_start: int               # 原 BED 源坐标（经两段回溯）
    src_end: int
    mid_start: int               # 中间坐标（第一段目标 = 第二段源）
    mid_end: int
    dest_start: int              # 最终坐标（第二段目标）
    dest_end: int


@dataclass
class PathCandidate:
    path_index: int              # 1 基路径序号（按发现顺序，稳定）
    chain1_id: str
    chain2_id: str
    mid_chrom: str               # 中间染色体（chain 文件 1 的目标侧）
    mid_strand: str
    dest_chrom: str              # 最终染色体（chain 文件 2 的目标侧）
    dest_strand: str
    fragments: list[PathFragment]
    seg1_covered: int            # 第一段覆盖的源碱基数（该路径的第一段 chain）
    seg2_covered: int            # 第二段覆盖的碱基数 = 完成两段迁移的源碱基数
    total: int                   # 源区间长度
    uncovered_src: list[tuple[int, int]] = field(default_factory=list)

    @property
    def final_coverage(self) -> float:
        """最终覆盖率：完成两段迁移的源碱基占源区间长度的比例（0~1）。"""
        return self.seg2_covered / self.total if self.total else 0.0


@dataclass
class RowResult2:
    line_no: int
    src_chrom: str
    src_start: int
    src_end: int
    raw: str
    status: str                          # mapped / partial / ambiguous / unmapped
    failed_stage: str = ""               # "" / "segment1" / "segment2"
    reason: str = ""
    candidates: list[PathCandidate] = field(default_factory=list)
    # 行级覆盖（供无映射行展示每段覆盖情况；有路径时以各路径自己的覆盖为准）
    seg1_covered: int = 0                # 第一段所有 chain 覆盖的并集
    seg2_covered: int = 0                # 所有可达路径最终覆盖的并集


def _union_cov(ranges) -> int:
    """若干 [s, e) 区间的并集长度。"""
    total = 0
    cur_s = cur_e = None
    for s, e in sorted(ranges):
        if cur_s is None:
            cur_s, cur_e = s, e
        elif s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
    if cur_s is not None:
        total += cur_e - cur_s
    return total


def _uncovered(lo: int, hi: int, frags: list[PathFragment]) -> list[tuple[int, int]]:
    """源区间 [lo, hi) 内未被任何最终片段覆盖的子区间（frags 已按源排序）。"""
    out: list[tuple[int, int]] = []
    cursor = lo
    for f in frags:
        if f.src_start > cursor:
            out.append((cursor, f.src_start))
        if f.src_end > cursor:
            cursor = f.src_end
    if cursor < hi:
        out.append((cursor, hi))
    return out


def format_uncovered(unc: list[tuple[int, int]]) -> str:
    """未覆盖源区间的紧凑文本表示，如 ``[100,150);[300,320)``。"""
    if not unc:
        return "."
    return ";".join(f"[{a},{b})" for a, b in unc)


def _summarize(notes: list[str], cap: int = 4) -> str:
    dedup = list(dict.fromkeys(notes))
    if len(dedup) > cap:
        return "；".join(dedup[:cap]) + f"；等 {len(dedup)} 处"
    return "；".join(dedup)


class TwoSegmentEngine:
    """同一 BED 区间依次经过两份 chain 的迁移。"""

    def __init__(
        self,
        chains1: list[Chain],
        side1: str,
        chains2: list[Chain],
        side2: str,
    ) -> None:
        self.engine1 = LiftOverEngine(chains1, side1)
        self.engine2 = LiftOverEngine(chains2, side2)
        self.side1 = side1
        self.side2 = side2

    def convert(self, rec: BedRecord) -> RowResult2:
        total = rec.end - rec.start
        row = RowResult2(
            line_no=rec.line_no,
            src_chrom=rec.chrom,
            src_start=rec.start,
            src_end=rec.end,
            raw=rec.raw,
            status="unmapped",
        )

        # ---- 第一段：源区间 → 中间坐标系 ----
        cands1 = self.engine1.candidates_for(rec.chrom, rec.start, rec.end)
        if not cands1:
            row.failed_stage = "segment1"
            if not self.engine1.has_chrom(rec.chrom):
                row.reason = (
                    f"第一段无映射（chain 文件 1）：染色体 {rec.chrom} "
                    f"在 chain 文件 1 的源侧不存在"
                )
            else:
                row.reason = (
                    f"第一段无映射（chain 文件 1）：[{rec.start},{rec.end}) 完全落在 "
                    f"chain 文件 1 所有 chain 的对齐范围之外（链内缺口或未覆盖区域），"
                    f"没有任何碱基可映射"
                )
            return row

        row.seg1_covered = _union_cov(
            (f.src_start, f.src_end) for c in cands1 for f in c.fragments
        )

        # ---- 第二段：只处理第一段实际映射到的中间片段 ----
        # 路径 = （第一段候选, 第二段 chain）组合；同一组合下不同中间片段
        # 命中的最终片段归入同一路径。
        paths: list[PathCandidate] = []
        path_keys: dict[tuple[int, int], PathCandidate] = {}
        dead_notes: list[str] = []
        n_mid = 0

        for ci1, cand1 in enumerate(cands1):
            for frag1 in cand1.fragments:
                n_mid += 1
                cands2 = self.engine2.candidates_for(
                    cand1.dest_chrom, frag1.dest_start, frag1.dest_end
                )
                if not cands2:
                    if not self.engine2.has_chrom(cand1.dest_chrom):
                        dead_notes.append(
                            f"中间染色体 {cand1.dest_chrom} 在 chain 文件 2 的源侧不存在"
                        )
                    else:
                        dead_notes.append(
                            f"中间片段 {cand1.dest_chrom}:[{frag1.dest_start},{frag1.dest_end}) "
                            f"落在 chain 文件 2 的对齐范围之外/链内缺口"
                        )
                    continue
                for cand2 in cands2:
                    key = (ci1, cand2.chain_uid)
                    path = path_keys.get(key)
                    if path is None:
                        path = PathCandidate(
                            path_index=len(paths) + 1,
                            chain1_id=cand1.chain_id,
                            chain2_id=cand2.chain_id,
                            mid_chrom=cand1.dest_chrom,
                            mid_strand=cand1.dest_strand,
                            dest_chrom=cand2.dest_chrom,
                            dest_strand=cand2.dest_strand,
                            fragments=[],
                            seg1_covered=cand1.covered,
                            seg2_covered=0,
                            total=total,
                        )
                        path_keys[key] = path
                        paths.append(path)
                    for frag2 in cand2.fragments:
                        # frag2.src_* 是中间坐标系子区间，按第一段片段的
                        # 块内仿射逆映射回溯到原 BED 源坐标（长度不变）。
                        lo = frag1.back_a * frag2.src_start + frag1.back_b
                        hi = frag1.back_a * frag2.src_end + frag1.back_b
                        s_lo, s_hi = (lo, hi) if lo <= hi else (hi, lo)
                        path.fragments.append(
                            PathFragment(
                                src_start=s_lo,
                                src_end=s_hi,
                                mid_start=frag2.src_start,
                                mid_end=frag2.src_end,
                                dest_start=frag2.dest_start,
                                dest_end=frag2.dest_end,
                            )
                        )

        if not paths:
            row.failed_stage = "segment2"
            row.reason = (
                f"第二段无映射（chain 文件 2）：第一段实际映射到的 {n_mid} 个中间片段"
                f"均未在 chain 文件 2 中获得对齐（{_summarize(dead_notes)}），"
                f"源区间没有任何碱基完成两段迁移"
            )
            return row

        alive = {k[0] for k in path_keys}
        dead1 = [c for i, c in enumerate(cands1) if i not in alive]
        survivor_note = ""
        if dead1:
            ids = ", ".join(c.chain_id for c in dead1)
            survivor_note = (
                f"；另有 {len(dead1)} 个第一段命中（chain {ids}）在第二段无映射，"
                f"未形成可达路径"
            )

        for path in paths:
            path.fragments.sort(key=lambda f: (f.src_start, f.src_end))
            path.seg2_covered = sum(f.src_end - f.src_start for f in path.fragments)
            path.uncovered_src = _uncovered(rec.start, rec.end, path.fragments)

        row.candidates = paths
        row.seg2_covered = _union_cov(
            (f.src_start, f.src_end) for p in paths for f in p.fragments
        )

        # 多条可达路径 => 歧义，全部保留，不自动选择或合并
        if len(paths) > 1:
            row.status = "ambiguous"
            row.reason = (
                f"存在 {len(paths)} 条可达路径（第一段 chain × 第二段 chain 组合），"
                f"按要求全部保留为歧义，不自动选择或合并；"
                f"各路径两段覆盖情况见 seg1_covered/seg2_covered 列" + survivor_note
            )
            return row

        path = paths[0]
        if path.seg2_covered == total:
            row.status = "mapped"
            row.reason = survivor_note.lstrip("；")
        else:
            row.status = "partial"
            row.reason = (
                f"{total - path.seg2_covered} bp 源序列未能完成两段迁移"
                f"（第一段覆盖 {path.seg1_covered}/{total}，"
                f"第二段覆盖 {path.seg2_covered}/{path.seg1_covered}）；"
                f"未覆盖源区间 {format_uncovered(path.uncovered_src)} 被丢弃，"
                f"残余片段为部分映射" + survivor_note
            )
        return row

    def convert_all(self, records: list[BedRecord]) -> list[RowResult2]:
        return [self.convert(r) for r in records]
