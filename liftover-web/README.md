# 离线基因组区间迁移 Web 应用（BED Liftover）

供实验室**内网/离线环境**使用的基因组坐标区间迁移工具。用户在浏览器中
上传 BED 文件与对应的 UCSC `chain` 文件，指定 BED 区间属于链的
**target 侧**还是**query 侧**，先在页面上逐条预览转换状态，确认后再
下载 TSV 结果。支持两种模式：

- **单段**：一份 chain 的常规迁移；
- **双段**：同一 BED 依次经过**两份 chain** 迁移（如 hg19→hg38→T2T），
  第二段只处理第一段实际映射到的片段，最终片段精确回溯到原 BED 源坐标。

- 后端：Python 3.11+ / FastAPI / Uvicorn
- 前端：原生 HTML + CSS + JavaScript（无任何 CDN / 外网依赖）
- 运行时数据不落盘：上传内容与转换结果只保存在进程内存中；校验失败时
  **不会**生成任何可下载文件

---

## 目录结构

```
<交付外层目录>/
├── Dockerfile            # 方式 A 构建入口：以外层目录为上下文（COPY liftover-web/...）
├── .dockerignore         # 外层构建上下文排除规则
└── liftover-web/
    ├── app/                  # 后端
    │   ├── bed.py            # BED 解析（零基左闭右开，错误带行号）
    │   ├── chain.py          # UCSC chain 解析（严格校验，错误带行号）
    │   ├── engine.py         # 单段区间转换引擎（缺口拆分 / 歧义 / 负链）
    │   ├── engine2.py        # 双段串联引擎（路径 / 回溯 / 失败阶段）
    │   ├── filtering.py      # 双段结果质量筛选（最小最终覆盖率 + 保留状态）
    │   ├── report.py         # 单段与双段 TSV 导出、筛选版 TSV、预览 JSON
    │   └── main.py           # FastAPI 入口（含歧义人工裁决与质量筛选接口）
    ├── static/               # 前端（纯本地资源）
    ├── examples/             # 可直接试用的 BED / chain 示例（含双段示例）
    ├── tests/                # pytest 测试（105 个用例）
    ├── scripts/init.sh       # 本地开发依赖初始化
    ├── Dockerfile            # 方式 B 构建入口：以本目录为上下文
    ├── docker-compose.yml    # 两种启动方式共用
    ├── .env.example
    ├── requirements.txt      # 运行依赖
    └── requirements-dev.txt  # 运行 + 测试依赖
```

---

## 一、Docker / Compose 启动（推荐，完全离线）

镜像在构建时安装全部 Python 依赖，运行阶段不需要访问外网。下面两种
入口构建的是**同一个镜像** `liftover-web:1.0.0`：应用代码与静态资源
全部打进镜像，缺少任一构建必需文件（`requirements.txt`、`app/`、
`static/`）时构建直接失败，不会启动缺少资源却显示健康的容器。

### 方式 A：在外层目录一键启动

在交付外层目录直接指定 `liftover-web/docker-compose.yml`，并以外层
目录为 Compose 项目目录：

```bash
# 构建并后台启动
docker compose --project-directory . -f liftover-web/docker-compose.yml up -d --build

# 查看状态 / 日志
docker compose --project-directory . -f liftover-web/docker-compose.yml ps
docker compose --project-directory . -f liftover-web/docker-compose.yml logs -f

# 停止
docker compose --project-directory . -f liftover-web/docker-compose.yml down
```

此时构建上下文是外层目录，使用外层目录下的 `Dockerfile`（COPY 源路径
带 `liftover-web/` 前缀）。

### 方式 B：进入 liftover-web 目录启动

```bash
cd liftover-web

# 构建并后台启动
docker compose up -d --build

# 查看状态 / 日志
docker compose ps
docker compose logs -f

# 停止
docker compose down
```

此时构建上下文是 `liftover-web/` 目录，使用该目录下的 `Dockerfile`。
两份 Dockerfile 仅 COPY 源路径前缀不同，构建产物一致。

访问地址（两种方式相同）：

- **http://localhost:8080** （应用页面）
- 健康检查：`http://localhost:8080/healthz`

### 端口与配置

复制环境变量模板并按需修改（方式 B 在 `liftover-web/` 内执行；方式 A
可将 `.env` 放在外层目录，或直接在命令行注入变量）：

```bash
cp .env.example .env
```

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `LIFTOVER_PORT` | `8080` | 宿主机映射端口（容器内固定监听 8080） |
| `LIFTOVER_MAX_UPLOAD_MB` | `50` | 单个上传文件大小上限（MB） |
| `LIFTOVER_JOB_TTL_SECONDS` | `3600` | 转换结果在内存中的保留时长（秒） |
| `LIFTOVER_MAX_JOBS` | `100` | 内存中同时保留的结果任务数量上限 |

改端口示例（两种入口写法都列出来）：

```bash
# 方式 A（外层目录）
LIFTOVER_PORT=9000 docker compose --project-directory . -f liftover-web/docker-compose.yml up -d

# 方式 B（liftover-web/ 内）
LIFTOVER_PORT=9000 docker compose up -d
```

随后访问 `http://localhost:9000`。

不用 Compose 时：

```bash
# 在 liftover-web/ 内构建
docker build -t liftover-web:1.0.0 .
# 或在外层目录构建（等价）
docker build -f Dockerfile -t liftover-web:1.0.0 .

docker run -d --name liftover-web -p 8080:8080 \
  -e LIFTOVER_MAX_UPLOAD_MB=100 \
  liftover-web:1.0.0
```

---

## 二、本地 Python 直接启动（开发 / 调试）

需要 Python 3.11+。

```bash
cd liftover-web

# 1) 初始化：创建虚拟环境并安装依赖（含测试依赖）
bash scripts/init.sh

# 2) 启动
.venv/bin/uvicorn app.main:app --host 0.0.0.0 --port 8080
```

访问 **http://localhost:8080**。

运行测试：

```bash
.venv/bin/python -m pytest -q
```

---

## 三、使用方法

### 单段转换

1. 转换模式选择「单段（一份 chain）」。
2. 选择 BED 文件（至少前 3 列 `chrom start end`，坐标为**零基左闭右开**，
   支持第 4 列起的额外列但不参与转换；允许 `#`、`browser`、`track` 行）。
3. 选择对应的 UCSC `chain` 文件。
4. 选择 BED 区间位于 chain 的哪一侧：
   - **target 侧**：target → query；
   - **query 侧**：query → target。
5. 点击「校验并转换」。页面显示每条源区间的状态。
6. 若存在**歧义**行，可先在页面上按 BED 物理行号逐行**人工裁决**
   （见下一节），也可以直接下载。
7. 点击「下载原始结果 TSV」获取未改动的原始结果；对歧义行做过裁决后，
   可点击「下载裁决版 TSV」获取只保留所选候选的裁决版。

### 双段转换（两份 chain 依次迁移）

1. 转换模式选择「双段（同一 BED 依次经过两份 chain）」。
2. 选择 BED 文件，再**按顺序**选择 chain 文件 1（第一段）与
   chain 文件 2（第二段），并分别指定每段的源侧：
   - `side1`：BED 区间属于 chain 文件 1 的哪一侧；
   - `side2`：第一段产物（中间坐标）属于 chain 文件 2 的哪一侧。
3. 点击「校验并转换」。语义规则：

   - 第一段把源区间迁到中间坐标系；**第二段只处理第一段实际映射到的
     片段**（第一段未覆盖的部分不参与第二段）；
   - 每个最终片段都按映射方向**回溯到原 BED 源坐标**（零基左闭右开），
     预览与 TSV 中同时给出源、中间、最终三套坐标；
   - 每条结果标明**两段 chain ID、路径序号与每段覆盖情况**
     （`seg1_covered` / `seg2_covered` / `total_len`）；
   - 任一段的缺口都会留下**源区间未覆盖信息**（`uncovered_src` 列），
     残余片段只会标为 `partial`，**绝不**标成完全映射；
   - 多条可达路径（第一段 chain × 第二段 chain 的组合）**全部保留为
     歧义**，不自动选择、不合并；在第二段走不通的第一段命中不构成
     可达路径，只在说明中注明；
   - 没有任何最终映射的源行仍保留**原物理行号**与**失败阶段**
     （`failed_stage` = `segment1` / `segment2`）；第二段找不到中间
     染色体同样按第二段无映射处理；
   - BED、chain 文件 1、chain 文件 2 中**任一**文件有格式错误，都会
     带行号列出并拒绝整次转换，不产生任何下载结果。

4. 歧义行的裁决与单段相同（按路径序号选择），下载按钮、
   原始/裁决版 TSV 的用法完全一致。
5. 双段结果还可做**质量筛选**（见下一节）：按每条路径的最终覆盖率
   与状态过滤后单独导出 TSV，不影响原始结果与裁决。

### 双段结果质量筛选

双段转换完成后，预览页底部出现「质量筛选」面板（仅双段结果可用）：

1. 输入**最小最终覆盖率**（0~1，含边界）：筛选按**每条路径**的
   `seg2_covered / total_len`（完成两段迁移的源碱基占源区间长度的比例）
   计算，达到阈值的路径通过；
2. 勾选要保留的**状态**（完全映射 / 部分映射 / 歧义；无映射行见下）；
3. 点击「应用筛选」，页面显示筛选后的预览与**各状态数量**，并可
   「下载筛选结果 TSV」（双段 25 列 + 末尾 `final_coverage` 列）。

筛选语义：

- **歧义行按路径分别判断**：通过阈值的路径全部保留（路径序号保持
  原值），系统**不自动替用户裁决**——即使只剩一条路径通过，行的状态
  仍保持 `ambiguous`；
- **无映射（unmapped）行始终计入统计**（各状态数量中恒为原结果的无
  映射行数，与是否勾选无关），但**绝不写入筛选结果**（预览行与 TSV）；
- 筛选只是原结果的导出视图：原始下载、裁决状态与版本号在筛选前后
  **完全不变**；参数非法（阈值越界、类型错误、状态集合为空或含非法
  值）时返回明确错误，原任务不受任何影响。

### 歧义区间人工裁决

歧义区间（`ambiguous`）不能自动进入下游分析，必须经人工裁决。
裁决**只作用于这一次转换结果**，不会改动原始转换：

- 在预览表的每个歧义行，候选按候选序号（1、2、…，与 TSV 的
  `candidate_index` 一致）列出，点候选前的单选框即提交裁决；
  点「撤销选择」可恢复为未裁决（重新输出全部候选）；
- 页面上方实时显示「已选/未裁决行数」与**结果版本号**；
- 同一个结果被连续修改时版本号递增；若两人/两个标签页同时拿着同一
  版本提交，只有一次成功，另一次收到 409 并自动刷新当前选择后重试；
- 原始下载（`GET /api/download/{token}`）的内容在裁决前后**字节不变**，
  行的原始 `status`（如 `ambiguous`/`partial`）也始终保留——裁决不会把
  部分映射改称完全映射。

### 状态定义

| 状态 | 含义 |
| --- | --- |
| 完全映射 `mapped` | 区间全部落在对齐块内（双段：全部源碱基经唯一可达路径完成两段迁移） |
| 部分映射 `partial` | 一部分落在链内缺口/未对齐区域；命中部分被拆成**按源坐标有序的片段**，未映射部分丢弃 |
| 歧义 `ambiguous` | 同一源区间命中**多个 chain / 多个目标位置**（双段：多条可达路径）；系统列出全部候选，**绝不自动挑选**，需人工逐行裁决后再用于下游分析 |
| 无映射 `unmapped` | 没有任何碱基可映射；结果中保留源行号、源区间与原因（双段另附 `failed_stage` 标明失败发生在第一段还是第二段） |

### 坐标与负链

- 输入与输出一律是**零基、左闭右开**（`start < end`）。
- chain 负链侧的存储坐标按 UCSC 语义（以染色体长度为枢轴）归一化后
  再做区间求交，因此**负链转换后的坐标仍然是 start < end 的零基左闭
  右开区间**；结果的 `dest_strand` 列标明该候选在目标坐标系中的方向。

### 下载文件格式（TSV）

UTF-8（带 BOM，Excel 友好）。

**原始结果**（`下载原始结果 TSV` / `GET /api/download/{token}`）列：

```
source_line, status, source_chrom, source_start, source_end,
candidate_index, chain_id, dest_chrom, dest_strand,
fragment_index, src_start, src_end, dest_start, dest_end, reason
```

- `source_line` 是源文件**物理行号**（从 1 开始）；
- 跨缺口区间每个片段占一行；歧义区间的每个候选、每个片段各占一行；
- `unmapped` 行保留源行号与原因，候选/片段列写 `.`。

**裁决版**（`下载裁决版 TSV` / `GET /api/download/{token}/adjudicated`）
在原始 15 列末尾追加一列 `decision`，取值：

| decision | 出现在 | 输出内容 |
| --- | --- | --- |
| `已裁决` | 已选定候选的歧义行 | **仅输出被选候选的全部片段**（`candidate_index` 沿用原序号） |
| `未裁决` | 尚未选择的歧义行 | 仍输出**全部候选**的全部片段 |
| `.` | mapped / partial / unmapped 行 | 与原始结果对应行**逐列一致** |

裁决版的 `status` 列与原始结果完全一致：歧义行即使选了候选仍写作
`ambiguous`，`partial` 也不会因为人工选择而变成 `mapped`。
文件名形如 `liftover_target_to_query_<时间戳>.adjudicated.v<版本>.tsv`。

**双段结果**（`POST /api/liftover2` 产生）列：

```
source_line, status, failed_stage, source_chrom, source_start, source_end,
path_index, chain1_id, chain2_id, mid_chrom, mid_strand,
dest_chrom, dest_strand, fragment_index,
src_start, src_end, mid_start, mid_end, dest_start, dest_end,
seg1_covered, seg2_covered, total_len, uncovered_src, reason
```

- `src_start`/`src_end` 是**回溯到原 BED 的源坐标**；`mid_*` 是中间坐标
  （第一段目标 = 第二段源）；`dest_*` 是最终坐标；
- `path_index` 为 1 基路径序号，`chain1_id`/`chain2_id` 是两段各自
  命中的 chain ID；`seg1_covered`/`seg2_covered`/`total_len` 给出该路径
  每段覆盖的碱基数与源区间总长；
- `uncovered_src` 列出该路径未覆盖的源子区间（如 `[100,120);[300,320)`），
  完全覆盖时为 `.`；
- `failed_stage` 仅对 `unmapped` 行有意义：`segment1`（第一段无映射）
  或 `segment2`（第二段无映射，含中间染色体在 chain 文件 2 中不存在），
  其余行写 `.`；无映射行的路径/片段列写 `.`，`uncovered_src` 为整个
  源区间，`seg1_covered`/`seg2_covered` 为行级覆盖（第一段各 chain 的
  并集 / 恒为 0）；
- 裁决版同样在末尾追加 `decision` 列（`已裁决`/`未裁决`/`.`），
  `candidate_index` 即路径序号，规则与单段一致。

**双段筛选版**（`POST /api/result/{token}/filter` 后
`GET /api/download/filtered/{filter_token}`）在双段 25 列末尾追加一列
`final_coverage`（该路径 `seg2_covered/total_len`，即筛选时使用的最终
覆盖率，如 `1`、`0.833333`），只包含通过筛选的行与路径；unmapped 行
不写入。

### 失败即阻止导出

- BED 坏行（列数不足、坐标非整数、`start >= end` 等）或 chain 格式错误
  （头列数/字段非法、块越界、缺少终止行、gap 越界等）都会在页面上
  **带行号列出全部问题**，接口返回 400，且**不会生成或保留任何可下载
  文件**；修正后需重新提交。
- 双段模式下 BED、chain 文件 1、chain 文件 2 任一文件出错同样整体
  拒绝，错误分别以 `BED 第 N 行`、`chain 文件 1：…`、`chain 文件 2：…`
  标明来源。

---

## 四、HTTP API 调用示例

页面的全部功能都通过以下 HTTP 接口完成（端口以 `LIFTOVER_PORT` 为准，
默认 8080；示例中的文件路径按方式 A 的外层目录给出，方式 B 下去掉
`liftover-web/` 前缀即可）。

### 健康检查

```bash
curl http://localhost:8080/healthz
# {"status":"ok"}
```

### 上传并转换（预览）

`POST /api/liftover`，`multipart/form-data`，字段：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `bed` | 文件 | BED 文件（零基左闭右开） |
| `chain` | 文件 | UCSC chain 文件 |
| `side` | 文本 | `target`（target → query）或 `query`（query → target） |

```bash
curl -s -X POST http://localhost:8080/api/liftover \
  -F "bed=@liftover-web/examples/example_target.bed" \
  -F "chain=@liftover-web/examples/demo.chain" \
  -F "side=target"
```

成功（200）返回预览与下载令牌：

```json
{
  "download_token": "……",
  "download_name": "liftover_target_to_query_20261005-120000.tsv",
  "preview": {
    "total": 6,
    "counts": {"mapped": 1, "partial": 1, "ambiguous": 1, "unmapped": 3},
    "rows": []
  }
}
```

输入校验失败（400）只返回带行号的错误列表，**不会**产生任何可下载
文件：

```json
{"errors": ["BED 第 3 行：……", "……"]}
```

### 双段转换（两份 chain 依次迁移）

`POST /api/liftover2`，`multipart/form-data`，字段：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `bed` | 文件 | BED 文件（零基左闭右开） |
| `chain1` | 文件 | 第一段 UCSC chain 文件 |
| `chain2` | 文件 | 第二段 UCSC chain 文件 |
| `side1` | 文本 | BED 区间属于 chain 文件 1 的哪一侧（`target`/`query`） |
| `side2` | 文本 | 第一段产物属于 chain 文件 2 的哪一侧（`target`/`query`） |

```bash
curl -s -X POST http://localhost:8080/api/liftover2 \
  -F "bed=@liftover-web/examples/example_double.bed" \
  -F "chain1=@liftover-web/examples/demo_seg1.chain" \
  -F "chain2=@liftover-web/examples/demo_seg2.chain" \
  -F "side1=target" -F "side2=target"
```

成功（200）返回预览与下载令牌（`kind` 为 `double`；单段接口返回
`single`）：

```json
{
  "download_token": "……",
  "download_name": "liftover_double_target-target_20261005-120000.tsv",
  "kind": "double",
  "preview": {
    "total": 8,
    "counts": {"mapped": 1, "partial": 2, "ambiguous": 1, "unmapped": 4},
    "rows": [
      {
        "line_no": 6, "status": "ambiguous", "failed_stage": "",
        "src_chrom": "chr1", "src_start": 300, "src_end": 400,
        "candidates": [
          {"path_index": 1, "chain1_id": "2", "chain2_id": "11",
           "mid_chrom": "chrB", "mid_strand": "-",
           "dest_chrom": "chrW", "dest_strand": "-",
           "seg1_covered": 100, "seg2_covered": 100, "total": 100,
           "uncovered_src": [],
           "fragments": [
             {"src_start": 300, "src_end": 400,
              "mid_start": 300, "mid_end": 400,
              "dest_start": 0, "dest_end": 100}
           ]}
        ]
      }
    ]
  }
}
```

- 预览行中 `candidates` 即**可达路径**（含两段 chain ID、路径序号、
  每段覆盖与未覆盖源区间）；`failed_stage` 标明无映射行的失败阶段；
- BED、chain 文件 1、chain 文件 2 任一解析失败都返回 400，错误分别
  以 `BED 第 N 行`、`chain 文件 1：…`、`chain 文件 2：…` 开头，整次
  转换被拒绝，不产生任何下载令牌；
- 下载、裁决接口与单段完全共用：`GET /api/download/{token}` 下载双段
  原始 TSV（25 列），`GET/POST /api/result/{token}/decisions` 按路径
  序号裁决，`GET /api/download/{token}/adjudicated` 下载裁决版
  （26 列）。

### 下载结果 TSV

`GET /api/download/{download_token}`。token 来自上一步响应，只在转换
完全成功后才存在，且结果只保存在服务进程内存中：

```bash
TOKEN=<上一步返回的 download_token>
curl -OJ "http://localhost:8080/api/download/${TOKEN}"
```

- 成功（200）：`text/tab-separated-values`（UTF-8 带 BOM），
  `Content-Disposition` 携带上一步返回的 `download_name` 文件名；
- token 不存在或已过期（404）：
  `{"errors": ["结果不存在或已过期，请重新提交转换"]}`。

### 查看当前裁决选择

`GET /api/result/{download_token}/decisions`。返回某次结果的版本号与
全部歧义行的当前选择（`candidate_index` 为 `null` 表示未裁决）：

```bash
curl -s "http://localhost:8080/api/result/${TOKEN}/decisions"
```

```json
{
  "download_token": "……",
  "version": 0,
  "decisions": [
    {"source_line": 11, "candidate_index": null, "candidate_count": 2}
  ]
}
```

### 提交人工裁决（选定候选 / 撤销）

`POST /api/result/{download_token}/decisions`，`application/json`：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `version` | 整数 | 客户端所基于的**当前结果版本**（首次为 0） |
| `source_line` | 整数 | BED 源文件**物理行号**，该行必须是歧义行 |
| `candidate_index` | 整数或 `null` | 选定的 1 基候选序号；`null` 表示撤销选择 |

```bash
# 为第 11 行选定候选 2（基于版本 0）
curl -s -X POST "http://localhost:8080/api/result/${TOKEN}/decisions" \
  -H 'Content-Type: application/json' \
  -d '{"version":0,"source_line":11,"candidate_index":2}'
```

成功修改（200）后版本递增：

```json
{"download_token":"……","version":1,"modified":true,
 "source_line":11,"candidate_index":2}
```

提交与当前选择相同（无实际变化）时同样返回 200，但 `modified` 为
`false`、`version` 不变。随后可在版本 1 上改选候选 1 或撤销：

```bash
# 撤销第 11 行的选择（基于版本 1）
curl -s -X POST "http://localhost:8080/api/result/${TOKEN}/decisions" \
  -H 'Content-Type: application/json' \
  -d '{"version":1,"source_line":11,"candidate_index":null}'
```

明确拒绝的情形：

| 场景 | 状态码 | 说明 |
| --- | --- | --- |
| token 不存在 / 结果已过期 | 404 | 需重新提交转换 |
| `version` 与服务端当前版本不一致 | **409** | 相同版本的并发修改只有一次成功；响应体带 `current_version`，查看最新选择后用新版本重试 |
| `source_line` 不存在于本次结果 | 400 | 错误信息带行号 |
| 该行不是歧义行（mapped/partial/unmapped） | 400 | 只有 `ambiguous` 行允许裁决 |
| `candidate_index` 越界（0、负数或超过候选数） | 400 | 错误信息标明允许范围 1..N |
| 请求体不是 JSON 对象 / 缺字段 / 字段类型错误 | 400 | 如 `version` 传字符串、`candidate_index` 传 `true` |

### 下载裁决版 TSV

`GET /api/download/{download_token}/adjudicated`，无需额外参数，
始终按服务端当前裁决状态即时生成：

```bash
curl -OJ "http://localhost:8080/api/download/${TOKEN}/adjudicated"
# 文件名：liftover_target_to_query_<时间戳>.adjudicated.v<版本>.tsv
```

- 已选歧义行只输出被选候选的全部片段（`decision=已裁决`）；
- 未选歧义行输出全部候选并标明 `decision=未裁决`；
- 其他状态行照原结果输出（`decision=.`），`status` 列保持原值；
- 与原始下载互不影响：`GET /api/download/{token}` 仍是裁决前的内容。

### 双段结果质量筛选

`POST /api/result/{download_token}/filter`，`application/json`，仅适用于
双段转换结果（`POST /api/liftover2` 返回的 token）：

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `min_coverage` | 数值 | 最小最终覆盖率，0~1（含边界）；按每条路径的 `seg2_covered/total_len` 判定 |
| `keep_statuses` | 字符串数组 | 要保留的状态，元素取自 `mapped`/`partial`/`ambiguous`/`unmapped`，至少一个 |

```bash
# 保留最终覆盖率 ≥ 0.9 的完全映射/部分映射/歧义行
curl -s -X POST "http://localhost:8080/api/result/${TOKEN}/filter" \
  -H 'Content-Type: application/json' \
  -d '{"min_coverage":0.9,"keep_statuses":["mapped","partial","ambiguous"]}'
```

成功（200）返回筛选后的预览、各状态数量与筛选下载令牌：

```json
{
  "download_token": "……",
  "filter_token": "……",
  "download_name": "liftover_double_target-target_20261005-120000.filtered.min-cov-0.9.tsv",
  "download_url": "/api/download/filtered/……",
  "min_coverage": 0.9,
  "keep_statuses": ["mapped", "partial", "ambiguous"],
  "counts": {"mapped": 1, "partial": 0, "ambiguous": 1, "unmapped": 4},
  "preview": {
    "total": 2,
    "counts": {"mapped": 1, "partial": 0, "ambiguous": 1, "unmapped": 4},
    "rows": []
  }
}
```

- `counts` 为各状态统计：`unmapped` **始终计入**（原结果的无映射行数，
  与 `keep_statuses` 是否包含它无关）；`preview.rows` / `preview.total`
  只含写入筛选结果的行，**unmapped 行绝不写入**；
- 歧义行按路径分别判断，通过阈值的路径全部保留（`path_index` 保持原
  序号），状态仍为 `ambiguous`，不自动裁决；
- 筛选不影响原任务：原始/裁决版下载字节与裁决版本号在筛选前后不变。

下载筛选版 TSV（双段 25 列 + 末尾 `final_coverage` 列，UTF-8 带 BOM）：

```bash
FILTER_TOKEN=<上一步返回的 filter_token>
curl -OJ "http://localhost:8080/api/download/filtered/${FILTER_TOKEN}"
```

明确拒绝的情形（均不改变原任务）：

| 场景 | 状态码 | 说明 |
| --- | --- | --- |
| token 不存在 / 结果已过期 | 404 | 需重新提交转换 |
| 该 token 是单段转换结果 | 400 | 质量筛选仅适用于双段结果 |
| `min_coverage` 越界（<0 或 >1） | 400 | 错误信息标明允许范围 [0, 1] |
| `min_coverage` 类型非法（字符串、`true`、NaN 等） | 400 | 必须是 0~1 的有限数值 |
| `keep_statuses` 缺失 / 非数组 / 为空 / 含非法状态 | 400 | 错误信息列出允许值 |
| 请求体不是 JSON 对象 | 400 | |
| filter_token 不存在 / 筛选结果已过期 | 404 | 需重新提交筛选 |

---

## 五、快速试用

```bash
# 页面上分别选择：
#   BED:   examples/example_target.bed
#   chain: examples/demo.chain
#   侧别:  target 侧
```

该示例覆盖全部四种状态（完全映射 / 跨缺口拆分的部分映射 / 歧义 /
无映射）。`examples/example_query.bed` 演示 query → target，
包含负链 query 的反向转换；`examples/example_bad.bed` 演示错误定位
与导出拦截。

歧义行（示例 BED 第 11 行 `chr1 320 360`，候选 1=chrA、候选 2=chrB）
可直接在预览页裁决，或用上面的 `POST /api/result/{token}/decisions`
接口按物理行号选定 / 撤销，再下载 `…/adjudicated` 裁决版 TSV。

### 双段快速试用

```bash
# 页面上选择「双段」模式，然后：
#   BED:        examples/example_double.bed
#   chain 文件 1: examples/demo_seg1.chain   （chr1 → chrA/chrB/chrM）
#   chain 文件 2: examples/demo_seg2.chain   （chrA → chrZ；chrB → chrW/chrV）
#   side1/side2: 均为 target
```

或用 API：

```bash
TOKEN=$(curl -s -X POST http://localhost:8080/api/liftover2 \
  -F "bed=@liftover-web/examples/example_double.bed" \
  -F "chain1=@liftover-web/examples/demo_seg1.chain" \
  -F "chain2=@liftover-web/examples/demo_seg2.chain" \
  -F "side1=target" -F "side2=target" | python3 -c "import json,sys; print(json.load(sys.stdin)['download_token'])")

# 下载双段原始 TSV（25 列，含两段 chain ID、路径序号、每段覆盖与未覆盖源区间）
curl -OJ "http://localhost:8080/api/download/${TOKEN}"

# 歧义行（BED 第 6 行，路径 1=chrW、路径 2=chrV）选定路径 2 后下载裁决版
curl -s -X POST "http://localhost:8080/api/result/${TOKEN}/decisions" \
  -H 'Content-Type: application/json' \
  -d '{"version":0,"source_line":6,"candidate_index":2}'
curl -OJ "http://localhost:8080/api/download/${TOKEN}/adjudicated"

# 质量筛选：保留最终覆盖率 ≥ 0.9 的 mapped/partial/ambiguous 行并下载筛选版
FILTER_TOKEN=$(curl -s -X POST "http://localhost:8080/api/result/${TOKEN}/filter" \
  -H 'Content-Type: application/json' \
  -d '{"min_coverage":0.9,"keep_statuses":["mapped","partial","ambiguous"]}' \
  | python3 -c "import json,sys; print(json.load(sys.stdin)['filter_token'])")
curl -OJ "http://localhost:8080/api/download/filtered/${FILTER_TOKEN}"
```

该示例 8 行数据覆盖：完全映射（行 3）、第二段缺口导致的部分映射
（行 4）、第一段+第二段缺口叠加的部分映射（行 5）、双路径歧义
（行 6）、第二段无映射（行 7 中间染色体不存在、行 8 落在对齐范围外）、
第一段无映射（行 9 范围外、行 10 无此染色体）。

---

## 六、离线部署说明

- 前端资源全部由本服务在 `/static/` 下提供，页面不含任何外部字体、
  脚本或样式表；
- Docker 构建需在有网络的机器上 `docker build`（或导出镜像），之后把
  镜像 `docker save | gzip > liftover-web.tgz` 拷贝到内网，
  `docker load` 后即可运行；
- 容器以非 root 用户运行，只暴露 8080 端口，不挂载任何数据卷；
- 上传文件与结果仅存内存，进程重启即清除，适合处理受限数据。
