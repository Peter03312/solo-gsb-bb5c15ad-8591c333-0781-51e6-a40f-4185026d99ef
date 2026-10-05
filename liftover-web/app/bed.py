"""BED3 解析器。

坐标按零基、左闭右开解释。只取前 3 列（chrom/start/end），
其余列保留在 ``extra`` 中但不参与转换。

解析采用“收集全部错误”策略：坏行不会立即中断，调用方可以一次性
向用户展示所有问题；只有错误列表为空时结果才可用。
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class BedRecord:
    line_no: int          # 物理文件行号，从 1 开始
    chrom: str
    start: int
    end: int
    raw: str              # 原始行（去掉换行符）
    extra: str            # 第 4 列起的原始文本


@dataclass
class ParseResult:
    records: list[BedRecord]
    errors: list[str]


def _parse_uint(text: str, line_no: int, field: str) -> int:
    """解析非负十进制整数；显式拒绝 '+5'、'1_000'、小数和空白。"""
    t = text.strip()
    if not t or not t.isdigit():
        raise ValueError(f"BED 第 {line_no} 行: {field} 必须是非负整数，得到 {text!r}")
    return int(t)


def parse_bed(text: str) -> ParseResult:
    records: list[BedRecord] = []
    errors: list[str] = []

    for idx, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        # browser/track 指令行与 # 注释行允许出现但不参与转换
        if line.startswith(("#", "browser", "track")):
            continue

        fields = line.split()
        if len(fields) < 3:
            errors.append(
                f"BED 第 {idx} 行: 至少需要 3 列 (chrom start end)，实际 {len(fields)} 列"
            )
            continue

        chrom = fields[0]
        try:
            start = _parse_uint(fields[1], idx, "start")
            end = _parse_uint(fields[2], idx, "end")
        except ValueError as exc:
            errors.append(str(exc))
            continue

        if start >= end:
            errors.append(
                f"BED 第 {idx} 行: start ({start}) 必须小于 end ({end})"
            )
            continue

        records.append(
            BedRecord(
                line_no=idx,
                chrom=chrom,
                start=start,
                end=end,
                raw=line,
                extra="\t".join(fields[3:]),
            )
        )

    if not records and not errors:
        errors.append("BED 文件为空（没有任何有效数据行）")

    return ParseResult(records=records, errors=errors)
