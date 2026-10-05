"""区间转换引擎（零基、左闭右开）。

支持把区间从 chain 的 target 侧转到 query 侧，或反向转换。

核心做法：

1. 把每块在负链帧下书写的坐标按 UCSC 语义镜像为正向帧；
2. 对每个与源区间相交的对齐块做仿射映射，源区间被链内缺口
   自然切成有序片段；
3. 跨多个 chain 命中时，所有候选都保留并标记“歧义”，
   绝不自动挑选。
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass, field

from .bed import BedRecord
from .chain import Chain


@dataclass(frozen=True)
class Fragment:
    src_start: int
    src_end: int
    dest_start: int
    dest_end: int
    # 块内仿射逆映射参数：src = back_a * dest + back_b（back_a ∈ {+1, -1}）。
    # 每个片段都落在单个对齐块内，因此目标侧子区间可以精确回溯到源坐标；
    # 双段串联迁移靠它把第二段命中的中间子区间折算回原 BED 源坐标。
    back_a: int = 1
    back_b: int = 0


@dataclass
class Candidate:
    chain_id: str
    dest_chrom: str
    dest_strand: str           # 该 chain 在目标坐标系一侧的方向
    fragments: list[Fragment]
    covered: int
    total: int
    # 内部 prepared chain 的稳定标识（同一 chain 多次命中间隔不变），
    # 双段串联时用于把不同中间片段命中同一条 chain 的片段归入同一路径。
    chain_uid: int = 0


@dataclass
class RowResult:
    line_no: int
    src_chrom: str
    src_start: int
    src_end: int
    raw: str
    status: str                # mapped / partial / ambiguous / unmapped
    reason: str = ""
    candidates: list[Candidate] = field(default_factory=list)


@dataclass(frozen=True)
class _NormBlock:
    s_start: int               # 源侧正向帧
    s_end: int
    d_lo: int                  # 目标侧正向帧（d_lo < d_hi）
    d_hi: int
    same_orientation: bool     # 源与目标在该 chain 上同向（同为 + 或同为 -）


@dataclass
class _PreparedChain:
    chain_id: str
    src_chrom: str
    dest_chrom: str
    dest_strand: str
    s_lo: int
    s_hi: int
    blocks: tuple[_NormBlock, ...]
    starts: tuple[int, ...]

    @classmethod
    def prepare(cls, chain: Chain, side: str) -> "_PreparedChain":
        if side == "target":
            src_strand, src_lo, src_hi = chain.t_strand, chain.t_start, chain.t_end
            dst_strand = chain.q_strand
            src_chrom, dest_chrom = chain.t_name, chain.q_name
        else:
            src_strand, src_lo, src_hi = chain.q_strand, chain.q_start, chain.q_end
            dst_strand = chain.t_strand
            src_chrom, dest_chrom = chain.q_name, chain.t_name

        norms: list[_NormBlock] = []
        # 负链帧镜像的枢轴是该侧染色体长度（UCSC 语义：
        # forward_coord = size - stored_coord），不是 (start+end)。
        if side == "target":
            src_size = chain.t_size
            dst_size = chain.q_size
            dst_neg = chain.q_strand == "-"
        else:
            src_size = chain.q_size
            dst_size = chain.t_size
            dst_neg = chain.t_strand == "-"
        src_neg = src_strand == "-"
        # 源侧覆盖范围也要归一化到正向帧，供快速排除使用
        norm_src_lo, norm_src_hi = src_lo, src_hi
        if src_neg:
            norm_src_lo, norm_src_hi = src_size - src_hi, src_size - src_lo

        for b in chain.blocks:
            ss, se, ds, de = b.t_start, b.t_end, b.q_start, b.q_end
            if side == "query":
                ss, se, ds, de = b.q_start, b.q_end, b.t_start, b.t_end
            if src_neg:
                ss, se = src_size - se, src_size - ss
            if dst_neg:
                ds, de = dst_size - de, dst_size - ds
            d_lo, d_hi = min(ds, de), max(ds, de)
            norms.append(
                _NormBlock(ss, se, d_lo, d_hi, src_neg == dst_neg)
            )

        # 源帧可能是倒序书写的（目标侧为负链时），统一按源坐标排序
        norms.sort(key=lambda nb: nb.s_start)
        return cls(
            chain_id=chain.chain_id,
            src_chrom=src_chrom,
            dest_chrom=dest_chrom,
            dest_strand=dst_strand,
            s_lo=norm_src_lo,
            s_hi=norm_src_hi,
            blocks=tuple(norms),
            starts=tuple(nb.s_start for nb in norms),
        )

    def map_interval(self, start: int, end: int) -> list[Fragment]:
        """把源区间与对齐块求交并映射。

        块已统一到两侧正向帧。源与目标相对方向为同向时块内是平移；
        反向时块内镜像：``dest = d_hi + s_start - src``。
        无论哪种情况，返回片段的目标坐标都满足 d_start < d_end
        （零基左闭右开）。
        """
        frags: list[Fragment] = []
        # 块按源坐标有序且互不相交。起点 <= start 的最后一块可能跨过 start；
        # 若它在 start 前已结束，则起点等于 start 的块（在下一个位置）才是
        # 第一候选块。
        i = bisect_right(self.starts, start) - 1
        if i < 0 or self.blocks[i].s_end <= start:
            i += 1
        while i < len(self.blocks):
            nb = self.blocks[i]
            if nb.s_start >= end:
                break
            ovl_lo = max(start, nb.s_start)
            ovl_hi = min(end, nb.s_end)
            if ovl_lo < ovl_hi:
                if nb.same_orientation:
                    delta = nb.d_lo - nb.s_start
                    d_start, d_end = ovl_lo + delta, ovl_hi + delta
                    back_a, back_b = 1, nb.s_start - nb.d_lo
                else:
                    # 反向：归一化前存储帧 dst_stored = d_bstart + (src - s_bstart)，
                    # 归一化后 dst = d_hi + s_start - src。
                    # 区间 [ovl_lo, ovl_hi) 因此映射到：
                    pivot = nb.d_hi + nb.s_start
                    d_start = pivot - ovl_hi
                    d_end = pivot - ovl_lo
                    back_a, back_b = -1, pivot
                frags.append(
                    Fragment(
                        src_start=ovl_lo,
                        src_end=ovl_hi,
                        dest_start=d_start,
                        dest_end=d_end,
                        back_a=back_a,
                        back_b=back_b,
                    )
                )
            i += 1
        return frags


class LiftOverEngine:
    def __init__(self, chains: list[Chain], side: str) -> None:
        if side not in ("target", "query"):
            raise ValueError("side 必须是 'target' 或 'query'")
        self.side = side
        prepared = [_PreparedChain.prepare(c, side) for c in chains]
        self._index: dict[str, list[_PreparedChain]] = {}
        for pc in prepared:
            self._index.setdefault(pc.src_chrom, []).append(pc)
        for lst in self._index.values():
            lst.sort(key=lambda p: (p.s_lo, p.chain_id))

    def has_chrom(self, chrom: str) -> bool:
        """源侧是否存在该染色体（区分“无此染色体”与“落在缺口/未覆盖区”）。"""
        return chrom in self._index

    def candidates_for(self, chrom: str, start: int, end: int) -> list[Candidate]:
        """对任意源侧区间给出全部命中候选，不做状态判定。

        单段 ``convert`` 与双段串联引擎共用同一套命中逻辑，保证行为一致；
        命中多个 chain 时全部返回，由调用方决定如何标记歧义。
        """
        length = end - start
        candidates: list[Candidate] = []
        for pc in self._index.get(chrom, []):
            if pc.s_hi <= start or pc.s_lo >= end:
                continue
            frags = pc.map_interval(start, end)
            if not frags:
                continue
            covered = sum(f.src_end - f.src_start for f in frags)
            candidates.append(
                Candidate(
                    chain_id=pc.chain_id,
                    dest_chrom=pc.dest_chrom,
                    dest_strand=pc.dest_strand,
                    fragments=frags,
                    covered=covered,
                    total=length,
                    chain_uid=id(pc),
                )
            )
        return candidates

    def convert(self, rec: BedRecord) -> RowResult:
        row = RowResult(
            line_no=rec.line_no,
            src_chrom=rec.chrom,
            src_start=rec.start,
            src_end=rec.end,
            raw=rec.raw,
            status="unmapped",
        )
        length = rec.end - rec.start
        candidates = self.candidates_for(rec.chrom, rec.start, rec.end)

        if not candidates:
            if not self.has_chrom(rec.chrom):
                row.reason = f"染色体 {rec.chrom} 在 chain 文件的源侧不存在"
            else:
                row.reason = (
                    f"[{rec.start},{rec.end}) 完全落在所有 chain 的对齐范围之外"
                    f"（链内缺口或未覆盖区域），没有任何碱基可映射"
                )
            return row

        # 多个 chain 各自给出映射 => 歧义，保留全部候选，不挑选
        if len(candidates) > 1:
            row.status = "ambiguous"
            row.reason = (
                f"命中 {len(candidates)} 个 chain（id: "
                f"{', '.join(c.chain_id for c in candidates)}），存在多个目标位置；"
                f"按要求不自动挑选，全部候选列于下方"
            )
            row.candidates = candidates
            return row

        cand = candidates[0]
        row.candidates = [cand]
        if cand.covered == length:
            row.status = "mapped"
        else:
            row.status = "partial"
            n_gap = len(cand.fragments) - 1
            bits = [f"{length - cand.covered} bp 落在链内缺口/未对齐区域"]
            if n_gap >= 1:
                bits.append(f"区间被缺口拆成 {len(cand.fragments)} 个有序片段")
            bits.append("部分映射，未映射部分被丢弃")
            row.reason = "；".join(bits)
        return row

    def convert_all(self, records: list[BedRecord]) -> list[RowResult]:
        return [self.convert(r) for r in records]
