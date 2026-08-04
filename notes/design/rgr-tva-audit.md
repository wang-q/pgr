# rgr 功能梳理：tva 替代关系与 pgr 候选命令（2026-08-04）

> 目的：梳理 intspan 项目 `rgr` 二进制的 14 个子命令，与专门做 TSV 的
> `tva` 项目（27 子命令）核对功能，明确哪些已被 tva 替代、哪些仍有独立
> 价值，并评估哪些值得移植进 pgr。源码：`~/Scripts/intspan/src/cmd_rgr/`
> 与 `~/Scripts/tva/src/cmd_tva/`。

## 1. 结论摘要

* rgr 是"杂烩"：6 个通用 TSV 工具 + 8 个 range 专用命令，靠 `extract_rg`
  （行内取第一个合法 range）把 `.rg` 单列文件和 TSV 任意列统一起来。
* **6 个通用 TSV 命令全部被 tva 替代**（`dedup→uniq`、`keep→keep-header`、
  `md→to md`、`filter→filter`、`select→select`；`replace` 是唯一缺口，
  tva 无直接等价命令，但属于 tva 域）。
* rgr 的剩余价值全部在 **8 个 range 专用命令**上。按 pgr 一级命令以
  **输入格式**命名的规则，这批行级 range 命令应归属**新建的 `pgr rg`
  家族**（**只操作 `.rg` 单列文件**，决策见 §3.0），而不是并入 `runlist`。
  其中 3 个高价值（`count`、`prop`、runlist 过滤），1 个中价值（行级
  `span`），4 个低价值/不建议（`sort`、`field`、`merge`、`pl-2rmp`——
  merge 是找片段重复的旧实现，已被 pgr `sd` 家族替代，不再移植）。
  `runlist` 家族保持"输入 runlist JSON、集合级"的定位不变。

## 2. rgr × tva 覆盖矩阵

### 2.1 通用 TSV 命令（6 个，均归 tva）

| rgr 命令 | 功能 | tva 对应 | 关系 |
| :--- | :--- | :--- | :--- |
| `dedup` | 整行/字段去重（单遍哈希） | `uniq` | **被替代**（tva 更强：ignore-case/repeated/at-least/max/equiv） |
| `keep` | 保留头部 N 行，其余送外部命令 | `keep-header` | **被替代**（几乎同功能，tva 支持多行头） |
| `md` | TSV → Markdown 表格 | `to md` | **被替代**（tva 还有 `to csv` / `to xlsx`） |
| `replace` | 两列映射文件替换字段值 | 无直接命令（libs 有 `read_replaces`；`extend`/`mutate` 表达式可近似） | **部分替代 / 缺口** |
| `filter` | 字段 str/num 比较过滤 | `filter` | **被替代**（tva 表达式语言更系统） |
| `select` | 字段选择/重排（号或名） | `select` | **被替代** |

### 2.2 range 专用命令（8 个，tva 不覆盖）

| rgr 命令 | 功能 | pgr 落点（新 `rg` 家族） | 移植价值 |
| :--- | :--- | :--- | :--- |
| `field` | 从 chr/start/end 字段建 range | `gff rg` 覆盖思路 | 低（不建议） |
| `sort` | 按 range（chr/start/strand）排序 | `rg sort` | 中低（外部 `sort` 不可完整替代；若做则原生做） |
| `count` | 每条 range 与一组 range 的重叠数（lapper） | `rg count` | **高** |
| `prop` | 每条 range 与 runlist 的交集比例（--full 加 length/size） | `rg prop` | **高** |
| `runlist` | 按 runlist 过滤行（overlap/non-overlap/superset） | `rg runlist`（沿用 rgr 命令名） | **高** |
| `span` | 行级 trim/pad/shift/flank/excise（5p/3p，保行） | `rg span` | 中 |
| `merge` | 覆盖度阈值重叠图合并（petgraph，O(n²)） | 不做 | **不做**（找片段重复的旧实现，`pgr sd` 家族已覆盖该目标） |
| `pl-2rmp` | 两轮 merge+replace 管道（分片降 O(n²)） | 不做 | 不做（merge 的工程包装，随 merge 一起弃） |

## 3. 新家族 `pgr rg` 的候选命令设计

> `rg` 家族定位：**输入 `.rg` 单列文件（每行整行是一个 range），行级处理**，
> 输出 `.rg` 行、统计值或 runlist JSON。与 `runlist`（输入 runlist JSON、
> 集合级）严格区分。一级命令按输入格式命名，与 `fa`/`fas`/`psl`/`gff`
> 一致。

### 3.0 家族输入契约决策（2026-08-04，待定稿）

**决策 A：`pgr rg` 只操作 `.rg` 单列文件**（每行整行是一个 range），不处理
"TSV 行里含 range 字段"的混合输入，不引入 `-f/--field`、`-H/--header`、
`extract_rg` 双模式。

理由：

1. 家族 = 单格式契约（`fa`/`fas`/`psl`/`gff`/`runlist` 皆然）；TSV+range
   是混合格式，放进来会重蹈 rgr 杂烩。
2. TSV+range 的老工作流（links / `gars tsv`）由 **rgr 继续承担**（与 spanr
   同样作为外部工具保留）；通用 TSV 归 **tva**。分工：
   `pgr rg`（.rg）→ `rgr`（TSV+range）→ `tva`（通用 TSV）。
3. 每个命令省掉双模式（-f/-H/-s + 行内找 range），代码、测试、边界减半。
4. 若未来出现"给 TSV 行标注 range 统计"的硬需求，落点是 tva 侧加 range
   感知（tva 是 TSV 正主），而不是把 pgr rg 做成混合体。

后果：要用 `rg count` 标注一张基因表，先 `tva select` 抽出 range 列转
.rg 再进 pgr rg——显式一步，换来家族边界干净。

本决策影响 §2.2 表中各命令的输入契约（均改为 .rg 单列）与 §3.1 设计，
标记为**待定稿**：命令清单、命名、输出格式后续可能调整。

### 3.1 边界与迁移问题（待决策）

现有 `runlist cover`（.rg → runlist JSON）与 `runlist coverage`（.rg →
深度）的**输入就是 .rg**，按新体系逻辑上应属 `rg` 家族。`rg runlist`
这个名字预留给"按 runlist 过滤行"（rgr 同款命令名），因此转换命令不用
`rg runlist` 命名。

**已迁移（2026-08-04）**：`runlist cover/coverage` → **`rg cover`** /
**`rg coverage`**（操作名，与 rgr/spanr 的 cover 语义一致、无撞名）；
`runlist` 家族从此只收 runlist JSON 输入。迁移内容：新建
`cmd_pgr/rg/`（mod/cover/coverage）、注册 `pgr rg`、删除 runlist 的
cover/coverage 命令、测试迁至 `tests/cli_rg.rs`、新增 `docs/rg.md` 并
更新 `docs/runlist.md` 与内部引用（psl to-range 帮助、notes）。

### 3.2 高价值

**`pgr rg count`**——逐区间重叠计数（待定稿）

* 输入：target `.rg` 文件 + 一个或多个 `.rg` 区间文件；输出：每行追加
  `count`（该 range 被多少区间覆盖），即 `range<TAB>count`。
* 实现：pgr 已有 `coitrees` 依赖（paf/pbit 在用），每染色体建
  `BasicCOITree`，对每条 target range 查询计数，O(n log n)，比 rgr 的
  rust-lapper 更稳（病态长区间有保证）。也直接服务 mosdepth 讨论里
  "per-region 计数 ≈ thresholds.bed.gz" 的缺口。
* 与 `rg coverage` 互补：`coverage` 是"按位置"的深度，`count` 是"按区间"的
  重叠数。

**`pgr rg prop`**——交集比例（待定稿）

* 输入：runlist JSON + `.rg` 文件；输出：每行追加 `prop`（与 runlist 交集
  占比），`--full` 追加 `length`、`size`，即 `range<TAB>prop[<TAB>length
  <TAB>size]`。
* 实现：`IntSpan::intersect` + `cardinality` 现成（pgr 的 `json_to_set`
  已有）；`.rg` 行解析复用 `rg_to_set` 同款 `Range::from_str`。
* 场景：评估一组区间（如基因、链）被重复库/比对区间覆盖的比例。

**`pgr rg runlist`**（沿用 rgr 命令名，待定稿）

* 按 runlist 过滤 `.rg` 行：`overlap`（与 runlist 相交）/ `non-overlap` /
  `superset`（被 runlist 包含）。与 rgr 的 `runlist` 子命令同名同义，
  降低从 rgr 迁移用户的理解成本。
* 实现：与 prop 同一套 intersect 逻辑；行级输出保留原行结构。

### 3.3 中价值

**不做 `pgr rg merge` / `pl-2rmp`**

* 背景：`rgr merge` 是老的"片段重复（SD）发现"实现——把互相重叠 ≥ 阈值
  （0.95/0.98）的区间（links 文件中的 part）归并成代表区间，靠
  `replace` 应用映射。它是 SD 检测流程的一个中间环节，本身不是 SD 检测。
* 现状：pgr 已有完整的 `sd` 家族（`sd search/align/cluster/decompose/
  cover/cross/run`，见 notes/design/sd.md 与 notes/references/biser.md），
  覆盖了"找片段重复"的目标，且输出/语义比 rgr merge 更适合 pgr。
* 结论：**不移植**。`rg merge` 的方案讨论（单命令直接替换 vs 两命令
  镜像 rgr、链向处理、coitrees+union-find 优化）全部归档，不再推进；
  若未来出现"links 类多 part 归一"的真实需求，再按当时场景重新评估。

**`pgr rg span`（行级 span，待定稿）**

* 对 `.rg` 每行做 trim/pad/shift/flank/excise（5p/3p 方向），输出变换后的
  `.rg` 行（`--append` 可追加新 range 字段）。pgr 现有 `runlist span` 是
  JSON 集合级，二者语义不同。
* 场景：处理带链向的行（如 PSL/链文件），pgr 的 `Range` 已有
  `trim_5p/3p`、`shift_5p/3p`、`flank_5p/3p`（rgr span 直接复用这些）。

### 3.4 低价值 / 不建议

* `rg sort`（待定稿）：**外部 `sort` 不可完整替代**——排序键是解析出的
  `(chr, start, strand)` 三元组（`-k` 无法处理 `name.chr(strand)` 变长
  组合），且非法行置尾、`-H` 头保留、`-f` 指定字段、`--group` 分组均无
  对应（实测：带链向时 `chr1(-)` 被外部 sort 整体排到 `chr1(+)` 前、
  物种前缀 `S288c.I` 被当成主键、非法行被排到最前而非末尾）。对 pgr
  内部无加速价值（`coverage` 的事件排序自带），但若"输出有序 .rg 供
  用户/下游工具查看"的场景成立，pgr 原生实现成本很低（复用 `Range`
  解析 + `sort_by_key`，stable 排序顺带修正 rgr 的整行去重副作用），
  属中低优先级、可选做。
* `rg field`：`gff rg` 已覆盖"字段 → range"；通用 TSV 字段化属于 tva 域。
* `rg merge` / `pl-2rmp`：找片段重复的旧实现，`pgr sd` 家族已替代（§3.2）。
* 通用 TSV 六件套：不进 pgr（tva 是正主，重复实现无意义）。

## 4. 与命名体系 / 既有讨论的衔接

* 按此前归纳的 pgr 命名体系，**一级命令按输入格式命名**：新家族
  `rg`（输入 `.rg`/含 range 的 TSV）与 `fa`/`fas`/`psl`/`gff` 同规则；
  `count`/`prop`/`runlist`/`span`/`sort`/`merge` 是家族内**操作名**
  （不改变格式，行→行；`runlist` 沿用 rgr 的命令名，表示"按 runlist
  过滤行"）；若迁移 cover/coverage，用 `rg cover`/`rg coverage`（操作名，
  不用 `rg runlist`，避免与过滤命令撞名）。
* rgr 的教训：把"通用 TSV"和"range 专用"混在一个二进制里，靠
  `extract_rg` 启发式统一，命令语义变得模糊。pgr 新建 `rg` 家族后，
  输入域明确（`.rg`/含 range 字段的 TSV），通用 TSV 交给 tva，runlist
  JSON 交给 `runlist`，三者各司其职，避免重蹈"杂烩"覆辙。
* 与 mosdepth 讨论衔接：`coverage --per-base`（逐位点，待实现）是"按
  位置"维度；`rg count`（逐区间）是"按区间"维度；`count` 的 per-region
  输出即 mosdepth `thresholds.bed.gz` 的对应物。

## 5. 建议实施顺序

1. 搭建 `rg` 家族骨架（`make_subcommand`/`execute` + 注册 + `docs/rg.md`）。
2. `rg count`——coitrees 现成、需求直接（重叠计数）。
3. `rg prop`——IntSpan 现成，与 count 共用 `-f/--header` 输入层。
4. `rg runlist`（overlap 过滤）——与 prop 同基础，几行代码。
5. 行级 `rg span`、`rg sort`——看需求（5p/3p 方向变换、有序输出是否真有场景）。
6. `runlist cover/coverage` 迁移到 `rg cover`/`rg coverage`（§3.0 决策后）。
