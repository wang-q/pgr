# rgr 功能梳理：tva 替代关系与 pgr 候选命令（2026-08-04）

> 目的：梳理 intspan 项目 `rgr` 二进制的 14 个子命令，与专门做 TSV 的
> `tva` 项目（27 子命令）核对功能，明确哪些已被 tva 替代、哪些仍有独立
> 价值，并评估哪些值得移植进 pgr。源码：`~/Scripts/intspan/src/cmd_rgr/`
> 与 `~/Scripts/tva/src/cmd_tva/`。

## 0. 当前进度（2026-08-04）

* `pgr rg` 家族已落地：`cover` / `coverage`（自 `runlist` 迁入）、`count`
  （自 rgr 移植，coitrees 索引）、`prop`（自 rgr 移植，IntSpan 交集）。
* `runlist` 家族已移除 cover/coverage，只收 runlist JSON 输入。
* 测试迁移：`tests/cli_rg.rs`（14 个用例，含 rgr `command_count` 与
  `command_prop` 原始测试的 .rg 部分；TSV `-f` 部分按决策 A 有意不迁）。
* 文档：`docs/rg.md` 新建，`docs/runlist.md` 更新为 JSON-only。
* 基准：`pgr rg count` vs `rgr count` 快 ~3.4×、输出逐行一致、内存 1.6×
  （见 notes/benchmarks/bench-rg-count.md）。
* `pgr rg sort` 已落地（2026-08-04）：按 (chr, start, strand) 稳定排序，
  非法行置尾。理由：有序 `.rg` 是下游工具与人工查看的常见需求，且外部
  `sort` 无法复刻解析键语义（见 §3.4）。
* `pgr rg runlist` / `rg span` / `rg merge` 已落地（2026-08-04）：
  runlist 过滤（overlap/non-overlap/superset）、行级 span（trim/pad/
  shift/flank/excise，5p/3p）、merge（互覆盖聚类出映射）。merge 反转了
  此前"被 sd 家族替代、不做"的决定（用户指示补做 .rg 版）。
* **rg 家族全部命令已落地**，无待办。

## 1. 结论摘要

* rgr 是"杂烩"：6 个通用 TSV 工具 + 8 个 range 专用命令，靠 `extract_rg`
  （行内取第一个合法 range）把 `.rg` 单列文件和 TSV 任意列统一起来。
* **6 个通用 TSV 命令全部被 tva 替代**（`dedup→uniq`、`keep→keep-header`、
  `md→to md`、`filter→filter`、`select→select`；`replace` 是唯一缺口，
  tva 无直接等价命令，但属于 tva 域）。
* rgr 的剩余价值全部在 **8 个 range 专用命令**上。按 pgr 一级命令以
  **输入格式**命名的规则，这批行级 range 命令应归属**新建的 `pgr rg`
  家族**（**只操作 `.rg` 单列文件**，决策见 §3.0），而不是并入 `runlist`。
  其中 3 个高价值（`count`、`prop`、runlist 过滤）与 1 个中价值（行级
  `span`）均已移植；`sort`、`merge` 后续按需补做（sort 论证见 §3.4，
  merge 反转"不做"决定见 §3.3）；仅 `field`、`pl-2rmp` 不做。`runlist`
  家族保持"输入 runlist JSON、集合级"的定位不变。

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

| rgr 命令 | 功能 | pgr 落点（新 `rg` 家族） | 状态 |
| :--- | :--- | :--- | :--- |
| `field` | 从 chr/start/end 字段建 range | `gff rg` 覆盖思路 | 不做 |
| `sort` | 按 range（chr/start/strand）排序 | `rg sort` | **✅ 已实现**（稳定排序，非法行置尾） |
| `count` | 每条 range 与一组 range 的重叠数（lapper） | `rg count` | **✅ 已实现**（coitrees，快 ~3.4×） |
| `prop` | 每条 range 与 runlist 的交集比例（--full 加 length/size） | `rg prop` | **✅ 已实现** |
| `runlist` | 按 runlist 过滤行（overlap/non-overlap/superset） | `rg runlist`（沿用 rgr 命令名） | **✅ 已实现** |
| `span` | 行级 trim/pad/shift/flank/excise（5p/3p，保行） | `rg span` | **✅ 已实现** |
| `merge` | 覆盖度阈值重叠图合并（petgraph，O(n²)） | `rg merge`（.rg 版：coitrees + union-find） | **✅ 已实现**（反转此前"不做"决定，见 §3.3） |
| `pl-2rmp` | 两轮 merge+replace 管道（分片降 O(n²)） | 不做 | 不做（merge 的工程包装，.rg 版无需分片） |

## 3. 新家族 `pgr rg` 的候选命令设计

> `rg` 家族定位：**输入 `.rg` 单列文件（每行整行是一个 range），行级处理**，
> 输出 `.rg` 行、统计值或 runlist JSON。与 `runlist`（输入 runlist JSON、
> 集合级）严格区分。一级命令按输入格式命名，与 `fa`/`fas`/`psl`/`gff`
> 一致。

### 3.0 家族输入契约决策（2026-08-04，已实施）

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

本决策已落地于现有命令（cover/coverage/count 均为 .rg 单列输入、无
-f/-H 双模式）；命令清单、命名、输出格式已定稿，后续新增的
prop/runlist/span/sort 均遵守同一契约。

### 3.1 边界与迁移问题（已迁移）

现有 `runlist cover`（.rg → runlist JSON）与 `runlist coverage`（.rg →
深度）的**输入就是 .rg**，按新体系逻辑上应属 `rg` 家族。`rg runlist`
这个名字预留给"按 runlist 过滤行"（rgr 同款命令名），因此转换命令不用
`rg runlist` 命名。

**已迁移（2026-08-04）**：`runlist cover/coverage` → **`rg cover`** /
**`rg coverage`**（操作名，与 rgr/spanr 的 cover 语义一致、无撞名）；
`runlist` 家族从此只收 runlist JSON 输入。迁移内容：新建
`cmd_pgr/rg/`（mod/cover/coverage）、注册 `pgr rg`、删除 runlist 的
cover/coverage 命令、测试迁至 `tests/cli_rg.rs`、新增 `docs/rg.md` 并
更新 `docs/runlist.md` 与内部引用（psl to-rg 帮助、notes）。

### 3.2 高价值

**`pgr rg count`**——逐区间重叠计数（✅ 已实现，2026-08-04）

* 输入：target `.rg` 文件 + 一个或多个 `.rg` 区间文件；输出：每行追加
  `count`（该 range 被多少区间覆盖），即 `range<TAB>count`。
* 实现：`libs/runlist::RgIndex`（每染色体 `BasicCOITree`，inclusive 坐标），
  对每条 target range 查询计数，O(n log n)；非法行跳过（rgr parity）。
  也直接服务 mosdepth 讨论里"per-region 计数 ≈ thresholds.bed.gz"的缺口。
* 与 `rg coverage` 互补：`coverage` 是"按位置"的深度，`count` 是"按区间"的
  重叠数。
* 基准：1M 区间 + 100k target，`pgr rg count` 224–235 ms vs `rgr count`
  753–758 ms（快 ~3.4×），输出 sort 后 diff 为空；RSS 42 vs 26 MB。
  详见 notes/benchmarks/bench-rg-count.md。

**`pgr rg prop`**——交集比例（✅ 已实现，2026-08-04，二分版）

* 输入：runlist JSON + `.rg` 文件；输出：每行追加 `prop`（与 runlist 交集
  占比），`--full` 追加 `length`、`size`，即 `range<TAB>prop[<TAB>length
  <TAB>size]`。
* 实现：`IntSpan::covered`（对有序 span 两次二分定位重叠段，O(log n +
  k)）+ `range_prop`；`.rg` 行解析复用 `usable_range`
  守卫（与 `rg_to_set`/`count` 同款）。
* 输出与 `rgr prop` 在 S288c fixture 上逐字节一致（`--full` 同样）。
* 基准：100k target + 154k spans，48.7 ms vs `rgr prop` 5.82 s（快
  ~120×，优化前 6.2 s 持平；火焰图定位到 `intersect` 的 O(n) VecDeque
  搬移后换算法）。详见 notes/benchmarks/bench-rg-prop.md。
* 场景：评估一组区间（如基因、链）被重复库/比对区间覆盖的比例。

**`pgr rg runlist`**——按 runlist 过滤（✅ 已实现，2026-08-04；沿用
rgr 命令名）

* 按 runlist 过滤 `.rg` 行：`overlap`（与 runlist 相交）/ `non-overlap` /
  `superset`（被 runlist 包含）。与 rgr 的 `runlist` 子命令同名同义，
  降低从 rgr 迁移用户的理解成本。
* 实现：与 prop 同一套 intersect 逻辑；行级输出保留原行结构。

### 3.3 中价值

**`pgr rg merge`（已实现，2026-08-04；反转此前"不做"决定）**

* 背景：`rgr merge` 是老的"片段重复（SD）发现"实现——把互相重叠 ≥ 阈值
  （0.95/0.98）的区间（links 文件中的 part）归并成代表区间，靠
  `replace` 应用映射。它是 SD 检测流程的一个中间环节，本身不是 SD 检测。
* 现状：pgr 已有完整的 `sd` 家族（`sd search/align/cluster/decompose/
  cover/cross/run`，见 notes/design/sd.md 与 notes/references/biser.md），
  覆盖了"找片段重复"的目标，且输出/语义比 rgr merge 更适合 pgr。
* 用户 2026-08-04 指示补做 .rg 版：`rg merge` 对单列 `.rg` 范围做互覆盖
  （reciprocal ≥ `--coverage`，默认 0.95）聚类，输出 `range<TAB>merged`
  映射（代表 = 并集 cover `chr(+):min-max`，单例省略）。实现用 coitrees
  找候选邻居 + union-find，O(n log n + k)，无需 rgr 的 O(n²)+petgraph 与
  `pl-2rmp` 分片。
* rgr 原测试基于多 part TSV（II.links.tsv），按决策 A 不迁移；新增 .rg
  版测试。

**`pgr rg span`（行级 span，✅ 已实现，2026-08-04）**

* 对 `.rg` 每行做 trim/pad/shift/flank/excise（5p/3p 方向），输出变换后的
  `.rg` 行（`--append` 可追加新 range 字段）。pgr 现有 `runlist span` 是
  JSON 集合级，二者语义不同。
* 场景：处理带链向的行（如 PSL/链文件），pgr 的 `Range` 已有
  `trim_5p/3p`、`shift_5p/3p`、`flank_5p/3p`（rgr span 直接复用这些）。

### 3.4 低价值 / 不建议

* `rg sort`（✅ 已实现，2026-08-04）：**外部 `sort` 不可完整替代**——排序键是解析出的
  `(chr, start, strand)` 三元组（`-k` 无法处理 `name.chr(strand)` 变长
  组合），且非法行置尾、`-H` 头保留、`-f` 指定字段、`--group` 分组均无
  对应（实测：带链向时 `chr1(-)` 被外部 sort 整体排到 `chr1(+)` 前、
  物种前缀 `S288c.I` 被当成主键、非法行被排到最前而非末尾）。对 pgr
  内部无加速价值（`coverage` 的事件排序自带），但若"输出有序 .rg 供
  用户/下游工具查看"的场景成立，pgr 原生实现成本很低（复用 `Range`
  解析 + `sort_by_key`，stable 排序顺带修正 rgr 的整行去重副作用），
  属中低优先级、可选做。
* `rg field`：`gff rg` 已覆盖"字段 → range"；通用 TSV 字段化属于 tva 域。
* `rg merge`：原判"找片段重复的旧实现、`pgr sd` 家族已替代（§3.2）"，后按
  用户指示补做 .rg 版（§3.3，2026-08-04 已实现）；`pl-2rmp` 仍不做
  （merge 的工程包装，.rg 版无需分片）。
* 通用 TSV 六件套：不进 pgr（tva 是正主，重复实现无意义）。

## 4. 与命名体系 / 既有讨论的衔接

* 按此前归纳的 pgr 命名体系，**一级命令按输入格式命名**：新家族
  `rg`（输入 `.rg` 单列）与 `fa`/`fas`/`psl`/`gff` 同规则；
  `count`/`prop`/`runlist`/`span`/`sort`/`merge` 是家族内**操作名**
  （不改变格式，行→行；`runlist` 沿用 rgr 的命令名，表示"按 runlist
  过滤行"）；若迁移 cover/coverage，用 `rg cover`/`rg coverage`（操作名，
  不用 `rg runlist`，避免与过滤命令撞名）。
* rgr 的教训：把"通用 TSV"和"range 专用"混在一个二进制里，靠
  `extract_rg` 启发式统一，命令语义变得模糊。pgr 新建 `rg` 家族后，
  输入域明确（`.rg` 单列），通用 TSV 交给 tva（TSV+range 交给 rgr），
  runlist JSON 交给 `runlist`，三者各司其职，避免重蹈"杂烩"覆辙。
* 与 mosdepth 讨论衔接：`coverage --per-base`（逐位点，待实现）是"按
  位置"维度；`rg count`（逐区间）是"按区间"维度；`count` 的 per-region
  输出即 mosdepth `thresholds.bed.gz` 的对应物。

## 5. 建议实施顺序

1. ✅ 搭建 `rg` 家族骨架（`make_subcommand`/`execute` + 注册 + `docs/rg.md`）。
2. ✅ `rg count`——coitrees 现成、需求直接（重叠计数），含基准。
3. ✅ `runlist cover/coverage` 迁移到 `rg cover`/`rg coverage`。
4. ✅ `rg prop`——IntSpan 现成，含 rgr fixture 回归。
5. ✅ `rg runlist`（overlap 过滤）——与 prop 同基础。
6. ✅ `rg sort`——按 (chr, start, strand) 稳定排序、非法行置尾。
7. ✅ 行级 `rg span`——trim/pad/shift/flank/excise（5p/3p）。
8. ✅ `rg merge`——互覆盖聚类 + union-find（coitrees 找邻居）。
