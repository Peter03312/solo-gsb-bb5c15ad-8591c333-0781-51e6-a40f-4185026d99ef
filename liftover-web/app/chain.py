"""UCSC chain 文件解析器。

chain 头示例::

    chain 5527 score 1 139107393 + 9518082 139107393 chr1 248956422 + 10000 248946422 1

    <tName> <tSize> <tStrand> <tSize?> …

13 个字段（含开头的 'chain'），随后是若干对齐数据行::

    12       # size 大小的对齐块，其后 target/query 各有 12 与 3 的缺口
    3 3
    5        # 最后一个对齐块（只有 size）

负链坐标说明（UCSC 语义）：任一侧 strand 为 '-' 时，该侧块在文件里
以反向帧书写——size 不变，而 start 是到序列末端的距离。归一化方法::

    forward_start = start + end - stored_end
    forward_end   = start + end - stored_start
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Block:
    """存储坐标帧下的对齐块（即文件中的书写坐标）。"""

    t_start: int
    t_end: int
    q_start: int
    q_end: int


@dataclass
class Chain:
    chain_id: str
    t_name: str
    t_size: int
    t_strand: str
    t_start: int
    t_end: int
    q_name: str
    q_size: int
    q_strand: str
    q_start: int
    q_end: int
    score: int
    blocks: list[Block]
    header_line: int


@dataclass
class ChainParseResult:
    chains: list[Chain]
    errors: list[str]


def _uint(tok: str) -> int:
    if not tok.isdigit():
        raise ValueError
    return int(tok)


def parse_chain(text: str) -> ChainParseResult:
    chains: list[Chain] = []
    errors: list[str] = []

    current: Chain | None = None          # 已经成功起始、等待块/终止的 chain
    poisoned = False                       # 头损坏：吞掉后续数据行直到下一个 chain 头
    t_cursor = q_cursor = 0
    header_line = 0

    def finish(open_chain: Chain, tc: int, qc: int, line_no: int) -> None:
        if tc != open_chain.t_end or qc != open_chain.q_end:
            errors.append(
                f"chain 第 {line_no} 行: chain id={open_chain.chain_id} 的对齐块未到达 "
                f"tEnd/qEnd（累计 target {tc}/{open_chain.t_end}, query {qc}/{open_chain.q_end}），"
                f"缺少终止行或块长错误"
            )
        if not open_chain.blocks:
            errors.append(
                f"chain 第 {open_chain.header_line} 行: chain id={open_chain.chain_id} 没有任何对齐块"
            )

    lines = text.splitlines()
    for line_no, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        tokens = line.split()

        if tokens[0] == "chain":
            # 开始新 chain 前，先收口上一个
            if current is not None:
                if not poisoned:
                    finish(current, t_cursor, q_cursor, line_no)
                else:
                    errors.append(
                        f"chain 第 {header_line} 行: chain 头损坏且未读取到任何有效块"
                    )
            current, poisoned = None, False

            if len(tokens) != 13:
                errors.append(
                    f"chain 第 {line_no} 行: chain 头必须是 13 列，实际 {len(tokens)} 列"
                )
                poisoned = True
                header_line = line_no
                continue

            try:
                score = _uint(tokens[1])
                t_size, t_start, t_end = _uint(tokens[3]), _uint(tokens[5]), _uint(tokens[6])
                q_size, q_start, q_end = _uint(tokens[8]), _uint(tokens[10]), _uint(tokens[11])
            except ValueError:
                errors.append(f"chain 第 {line_no} 行: chain 头中的长度/坐标必须是非负整数")
                poisoned = True
                header_line = line_no
                continue

            t_name, t_strand = tokens[2], tokens[4]
            q_name, q_strand, chain_id = tokens[7], tokens[9], tokens[12]

            def head_err(msg: str) -> None:
                errors.append(f"chain 第 {line_no} 行: {msg}")

            if t_strand not in "+-" or q_strand not in "+-":
                head_err(f"strand 只能是 '+' 或 '-'，得到 tStrand={t_strand!r} qStrand={q_strand!r}")
                poisoned = True
                header_line = line_no
                continue
            if not (0 <= t_start < t_end <= t_size):
                head_err(f"target 范围非法，要求 0 <= tStart < tEnd <= tSize（{t_start},{t_end},{t_size}）")
                poisoned = True
                header_line = line_no
                continue
            if not (0 <= q_start < q_end <= q_size):
                head_err(f"query 范围非法，要求 0 <= qStart < qEnd <= qSize（{q_start},{q_end},{q_size}）")
                poisoned = True
                header_line = line_no
                continue
            if not chain_id.isdigit() or chain_id == "0":
                head_err(f"chain id 必须是正整数，得到 {chain_id!r}")
                poisoned = True
                header_line = line_no
                continue

            current = Chain(
                chain_id=chain_id,
                t_name=t_name,
                t_size=t_size,
                t_strand=t_strand,
                t_start=t_start,
                t_end=t_end,
                q_name=q_name,
                q_size=q_size,
                q_strand=q_strand,
                q_start=q_start,
                q_end=q_end,
                score=score,
                blocks=[],
                header_line=line_no,
            )
            t_cursor, q_cursor = t_start, q_start
            header_line = line_no
            continue

        # 数据行
        if current is None or poisoned:
            errors.append(
                f"chain 第 {line_no} 行: 数据行出现在有效 chain 头之外：{line!r}"
            )
            continue

        try:
            nums = [_uint(t) for t in tokens]
        except ValueError:
            errors.append(f"chain 第 {line_no} 行: 数据行只能是 1 列(size)或 3 列(size dt dq)非负整数")
            continue

        if len(nums) not in (1, 3):
            errors.append(f"chain 第 {line_no} 行: 数据行列数必须是 1 或 3，实际 {len(nums)} 列")
            continue

        size = nums[0]
        if size == 0:
            errors.append(f"chain 第 {line_no} 行: 对齐块 size 不能为 0")
            continue

        tb_start, qb_start = t_cursor, q_cursor
        tb_end, qb_end = t_cursor + size, q_cursor + size

        if tb_end > current.t_end or qb_end > current.q_end:
            errors.append(
                f"chain 第 {line_no} 行: 对齐块超出 chain 范围"
                f"（target [{tb_start},{tb_end})/{current.t_end}, "
                f"query [{qb_start},{qb_end})/{current.q_end}）"
                f"，chain id={current.chain_id} 作废"
            )
            current, poisoned = None, True
            continue

        current.blocks.append(Block(tb_start, tb_end, qb_start, qb_end))
        t_cursor, q_cursor = tb_end, qb_end

        if len(nums) == 1:
            # 终止行：块必须恰好顶到 chain 末端，否则整条 chain 作废
            if t_cursor != current.t_end or q_cursor != current.q_end:
                errors.append(
                    f"chain 第 {line_no} 行: 终止块未对齐 tEnd/qEnd"
                    f"（累计 target {t_cursor}/{current.t_end}, query {q_cursor}/{current.q_end}）"
                    f"，chain id={current.chain_id} 作废"
                )
            else:
                chains.append(current)
            current, poisoned = None, False
        else:
            dt, dq = nums[1], nums[2]
            # UCSC 允许最后一个非终止行的 gap 为 0（一侧最后一块紧贴链末端，
            # 另一侧仍有间隙）；但两个 gap 同时为 0 时该行没有存在意义，拒绝。
            if dt == 0 and dq == 0:
                errors.append(
                    f"chain 第 {line_no} 行: dt 与 dq 不能同时为 0"
                    f"（无间隙的数据行应直接写成终止行）"
                )
                continue
            t_cursor += dt
            q_cursor += dq
            if t_cursor > current.t_end or q_cursor > current.q_end:
                errors.append(
                    f"chain 第 {line_no} 行: gap 之后越过 chain 末端"
                    f"（target 游标 {t_cursor}/{current.t_end}, "
                    f"query 游标 {q_cursor}/{current.q_end}），后续块无法解析，"
                    f"chain id={current.chain_id} 作废"
                )
                poisoned = True
                current = None

    # 文件结束收口
    if current is not None:
        if poisoned:
            errors.append(f"chain 第 {header_line} 行: chain 头损坏且未读取到任何有效块")
        else:
            errors.append(
                f"chain 第 {len(lines) or 1} 行: chain id={current.chain_id}（头在第 "
                f"{current.header_line} 行）缺少单字段终止行（文件提前结束）"
            )

    if not chains:
        errors.append("chain 文件不包含任何完整的 chain（没有可用于转换的对齐块）")

    return ChainParseResult(chains=chains, errors=errors)
