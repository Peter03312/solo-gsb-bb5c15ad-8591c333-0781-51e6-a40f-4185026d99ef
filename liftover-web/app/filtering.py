"""双段结果的质量筛选（最小最终覆盖率 + 保留状态）。

筛选只作用于某一次双段转换结果（``POST /api/liftover2`` 产生）的导出
视图，绝不改动原始任务：原始下载字节、裁决状态与版本号在筛选前后完全
一致。规则：

- 覆盖率按**每条路径**计算：``seg2_covered / total``（完成两段迁移的源
  碱基占源区间长度的比例），达到 ``min_coverage``（含边界）即通过；
- 歧义行按路径分别判断，保留通过阈值的全部路径，**不自动替用户裁决**——
  即使只剩一条路径通过，行的状态仍保持 ``ambiguous``，路径序号
  （``path_index``）保持原值；
- ``unmapped`` 行始终计入统计（返回的 ``counts``），但绝不写入筛选结果
  （预览行与 TSV）；
- 筛选产物挂在独立的 filter token 下，与原任务互不影响。
"""

from __future__ import annotations

from dataclasses import replace

from .engine2 import RowResult2

# 允许勾选的状态集合。unmapped 是合法取值，但只影响统计口径说明：
# 无映射行始终计入统计、绝不写入筛选结果，与是否勾选无关。
FILTER_STATUSES = ("mapped", "partial", "ambiguous", "unmapped")


def filter_double_rows(
    rows: list[RowResult2],
    min_coverage: float,
    keep_statuses: set[str],
) -> tuple[list[RowResult2], dict[str, int]]:
    """按最小最终覆盖率与保留状态筛选双段结果。

    返回 ``(筛选后的行, 各状态统计)``：

    - 筛选后的行是原行或其拷贝（歧义行可能只保留部分路径），原任务的
      行对象不被改动；
    - 统计中 ``unmapped`` 恒为原结果的无映射行数（始终计入），其余状态
      为通过筛选的行数。
    """
    counts = {"mapped": 0, "partial": 0, "ambiguous": 0, "unmapped": 0}
    kept: list[RowResult2] = []
    for row in rows:
        if row.status == "unmapped":
            # 始终保留在统计中，但不写入筛选结果
            counts["unmapped"] += 1
            continue
        if row.status not in keep_statuses:
            continue
        passing = [p for p in row.candidates if p.final_coverage >= min_coverage]
        if not passing:
            continue
        counts[row.status] += 1
        if len(passing) == len(row.candidates):
            kept.append(row)
        else:
            # 歧义行只保留通过阈值的路径；状态与路径序号保持原值
            kept.append(replace(row, candidates=passing))
    return kept, counts
