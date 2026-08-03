# fas-multiz 设计与实现

> **状态：已实现** — `libs::fas_multiz` 核心库与 `pgr fas multiz` CLI 均已落地，详见 §2.10；
> 阶段 1（全体物种 SP 打分 + 确定性合并次序）已于 2026-08-03 落地，见 §2.11；
> 阶段 1.5（yama 引擎直译：C/D/I 三状态 + 准自然 gap + 参考锚定 LB/RB）已于 2026-08-03 落地，见 §2.12。

本文档是 `libs::fas_multiz` 的设计稿，基于对 multiz 源码的分析（见 [multiz.md](../references/multiz.md)）。
涵盖 pgr 的两条路径与 multiz 的关系（§1）以及 fas-multiz 的设计与实现（§2）。

## 1. pgr 的两条路径与 multiz 的关系

`pgr` 里与 multiz 相关的工作分成两条路径，定位不同，但**不存在"策略层面"的对立**：

1. **快速交集捷径**：`pgr pl p2m`（内部用 `fas cover` / `spanr intersect` / `fas slice` / `fas join` / `fas subset`）只做区间集合运算和机械堆叠，快速得到"所有输入共同覆盖"的交集比对，不做 gap 求解。
2. **完整比对引擎**：`pgr fas multiz`（`libs::fas_multiz`）把 multiz 的 yama 引擎直译为 `.fas` 层 profile 合并，与 multiz 是同一套比对策略，输出 union/mesh 风格的完整比对。

`p2m` 之所以"快"，不是因为它采用了另一种比对策略，而是因为它**主动放弃了 multiz 会做的两件事**——gap 冲突求解和边缘覆盖保留——只回答一个更窄的问题："大家共同覆盖的那段区域长什么样"。

### 1.1 快速交集捷径：p2m + join

*   **工具**: `pgr pl p2m`，内部依次调用 `fas cover` → `spanr intersect` → `fas slice` → `fas join` → `fas name` / `fas subset`（后者按 name 列表统一各 block 的物种列顺序）。
*   **实现细节**:
    1.  **锚定**: 取第一个输入的第一个物种作为参考（Reference Target）。
    2.  **交集计算**: `fas cover` 提取每个输入的参考区间（`--trim 10` 切边），`spanr intersect` 计算所有输入共有的基因组区域，`spanr merge` 把交集范围合并回各输入。
    3.  **切片**: `fas slice` 按交集范围从原始 Pairwise 文件中提取序列。
    4.  **堆叠**: `fas join` 以参考序列的坐标范围为 Key，将来自不同文件的 Block 机械地堆叠在一起。
*   **特点**:
    *   **极速**: 仅涉及 I/O 和坐标计算。
    *   **局限**: 不处理 Gap。它假设所有输入文件在同一坐标下的参考序列片段是完全一致的。若不同 Pairwise Alignment 中参考序列的 Gap 状态不一致，直接 Join 会导致非参考序列错位。
*   **定位**: 这是 multiz 类比对的**低成本近似**：当只需要严格交集（core）且参考骨架一致时，可以跳过 DP 直接拿结果。它不是与 multiz 对立的"另一种策略"。`fas refine` 是独立的可选后处理命令，不在 `p2m` 流程内。

### 1.2 完整比对引擎：fas-multiz

*   `pgr fas multiz` 将 multiz 的 yama 引擎直译为 `libs::fas_multiz`（C/D/I 三状态 + 准自然 GAP + LB/RB 参考锚定 + 全体物种对打分，见 §2.12）。
*   它与 multiz 的策略**相同**：以参考坐标为主轴做一体化 DP，在合并时实时解决参考 gap 冲突，输出 union/mesh 风格比对。
*   与 `p2m` 的关系：`p2m` 解决不了的两类情况（参考 gap pattern 不一致、边缘覆盖差异），正是 fas-multiz 处理的场景；代价是需要跑 DP，比 `p2m` 慢。

### 1.3 结论

*   只想要"所有输入共同覆盖的快速交集" → `pgr pl p2m`（无 DP）。
*   想要 multiz 同款 union/mesh 完整比对、处理 gap 冲突 → `pgr fas multiz`（yama 直译）。
*   两者是"快而粗"与"慢而全"的取舍，不是两种对齐哲学。

### 1.4 在 pgr 中是否需要 Yama DP

*   **结论：需要，且已实现**。`libs::fas_multiz` 就是 yama 的直接移植（§2.12），`pgr fas multiz` 提供 CLI。
*   `p2m` 不需要 DP，不是因为它代表了另一种策略，而是因为它明确放弃了 gap 求解：只求交集、只信参考骨架一致的部分。
*   Yama DP 解决的两类问题：
    1.  **参考序列 gap 冲突**：不同输入中参考序列的 gap pattern 不一致时，通过 DP 在合并过程中实时插入/调整 gap，使所有序列在同一参考坐标系下保持一致。
    2.  **合并取舍**：通过 sum-of-pairs 打分决定列如何合并、哪些区域保留，而不是机械堆叠。
    在 `p2m` 里，这两类问题被"交集窗口的选取"直接绕开（只选大家一致的区域）；在 fas-multiz 里则交给 yama DP 解决。
*   因此正确的表述是：**p2m = 快速交集捷径（无 DP）；fas-multiz = multiz 同策略的完整比对（有 DP）**。二者互补，不存在策略层面的核心差异。

### 1.5 pgr 中 multiz 前置链路：LASTZ 与链化

在 UCSC 的典型 WGA 流程中，`multiz` 位于"pairwise 比对 + 链化 + net + mafFromNet"之后，只消费已经整理好的 MAF。`pgr` 目前在这一前置链路上，也已经有相当完整的 Rust 封装，主要对应到：

*   `pgr align lastz`：LASTZ 前端
    *   位置：`src/cmd_pgr/align/lastz.rs`。
    *   作用：包装 `lastz`，生成 LAV 格式输出，参数设计对齐 Cactus/UCSC 风格。
    *   特点：
        *   内置 UCSC 风格 preset（`set01`..`set07`），包括常见物种组合（Hg vs Pan/Mm/Bos/DanRer 等），每个 preset 绑定一套参数串和一个 4x4 替换矩阵（通过临时文件写给 `Q=` 选项）。
        *   自动加上 `--format=lav`、`--markend`、`--ambiguous=iupac`、`--querydepth=keep,nowarn:N` 等选项，行为贴近 Cactus repeat-masking 里的 LASTZ 调用约定。
        *   支持单文件和目录递归：对 target/query 目录做递归扫描（`.fa`/`.fa.gz`），生成笛卡尔积 job 列表。
        *   使用 `rayon` 并行跑多个 lastz 进程，并为每个 target–query 组合生成类似 `[t]vs[q].lav` 的输出文件名（带冲突规避逻辑）。
    *   对 multiz 的意义：
        *   对应于"blastz/lastz pairwise 比对"这一步，为后续链化、net 和 multiz/fas-multiz 提供高质量的成对比对基础。

*   `pgr psl chain`：PSL 链化（axtChain 风格）
    *   位置：`src/cmd_pgr/psl/chain.rs`，调用 `libs::chain` 中的 DP 引擎。
    *   作用：把 PSL 对齐 block 链成较长的 syntenic chain，逻辑类似 UCSC 的 `axtChain`/`chainNet` 里的链化步骤。
    *   打分与 gap 模型：
        *   使用 `SubMatrix` 作为替换矩阵，默认 Identity（匹配 +100 / 不匹配 -100），也可通过 `--score-scheme` 选择 HoxD55 或读取 LASTZ 格式打分文件。
        *   gap 成本由 `GapCalc` 驱动：
            *   线性模式：`--gap-model loose|medium`，对应 Kent 源码中针对远缘/近缘物种的 quasi-natural gap 曲线。
            *   仿射模式：`--gap-open` + `--gap-extend` 显式指定 open/extend，内部通过 `GapCalc::affine` 生成 gap 曲线。
        *   链化 DP 中的评分公式与 UCSC axtChain 一致：`Score = BlockScore + max(PrevScore - GapCost)`。
    *   结构与实现：
        *   依据 `(t_name, q_name, strand)` 分组，内部使用 KD-tree 等结构（见 `libs::chain`）加速前驱 block 搜索。
        *   允许在有 2bit 序列的情况下，用 `ScoreContext` 和 `calc_block_score` 精确重算每个 block 的序列得分，而不是只依赖 PSL 自带分数。
    *   对 multiz 的意义：
        *   在 `pgr` 里，这一步提供了与 UCSC 链化阶段等价的"整理过的 syntenic 对齐骨架"，可以作为（通过 AXT/MAF/FA 转换后）fas-multiz 的上游输入。

综上，`pgr align lastz` + `pgr psl chain` 组合，大致覆盖了 UCSC 链路中 "blastz/lastz 比对 + axtChain 链化" 这两步。它们提供了 multiz/fas-multiz 所需的 pairwise 对齐基础，而 `libs::fas_multiz` 则承担了更上游的 profile 合并角色：在已经有 syntenic 对齐骨架的前提下，对多个 `.fas` profile 做带状 DP 合并，构建 union/mesh 风格的多序列比对。

## 2. fas-multiz 设计与实现

在 `pgr` 中，`fas` 是"块级多序列比对"的核心抽象：每个 block 表示一段参考坐标下的多物种比对。`libs::fas_multiz` 就是在这一抽象上实现的 multiz 类功能：围绕 `.fas` profile 做合并，而不是在 MAF 文本层面复刻 `multiz`。

### 2.1 目标与输入/输出

*   **目标**：
    *   给定多个 block FA 文件（例如多个 pairwise-derived `.fas` 或不同 pipeline 生成的 `.fas`），在共享参考物种的坐标系下，将它们合并为一个"union/mesh 风格"的 block FA。
    *   和现有 `p2m + join`（可选加 `fas refine`）所产出的 **Core/Intersection** 结果互为补充：一个偏交集（core），一个偏并集（union/mesh）。
*   **输入**：
    *   `k` 个 `.fas` 文件（`k >= 2`），它们的 block 中均包含同名的参考序列（例如 `ref`）。
    *   可选的：一个核心交集区域（来自 `fas cover` + `spanr`），用于限制计算范围。
*   **输出**：
    *   一个新的 `.fas` 文件，包含合并后的多序列比对 block：
        *   在交集区域内，行为应与当前 `p2m + join` 相兼容。
        *   在边缘/非完全交集区域内，会尽量保留来自不同输入的对齐（union 行为）。

### 2.2 相对于 multiz 的主要差异

*   **工作层级不同**：
    *   multiz 直接操作 MAF，对齐的是两个 MAF profile。
    *   `pgr` 中的 fas-multiz 将直接操作 block FA，对齐的是若干 `.fas` profile。
*   **上游数据准备不同**：
    *   在 `pgr` 中，pairwise MAF/AXT 等通常已经通过 `pgr maf/axt to-fas` 等步骤规整为更统一的 block FA 表达。
    *   这意味着 fas-multiz 可以假设输入已经经过一次"标准化"，不需要自己处理复杂的 MAF 语法和扩展行。
*   **与现有命令的关系**：
    *   fas-multiz 更像是 `pgr fas join` 的"智能版/DP 版"：
        *   `fas join`：根据参考坐标"机械堆叠"，不处理 gap 冲突。
        *   fas-multiz：在堆叠时引入 profile–profile 的 DP/启发式，解决参考 gap 冲突和 block 选择问题。
    *   输出仍然是 `.fas`，可以直接接上 `fas refine`, `fas stat`, `fas to-vcf` 等命令。

### 2.3 数据流设计

1.  **标准化输入**：
    *   所有 upstream 比对结果（MAF/AXT 等）先通过现有命令统一转为 `.fas`。
    *   如有需要，可加一步 `fas normalize`（对序列名、物种名、参考 ID 做统一）。
2.  **block 级别的配对与聚类**：
    *   按参考物种与坐标对 block 做分组，将"位置相近"的 block 视为候选合并单元。
    *   这一层可以重用 `fas cover` / `spanr` 得到的区间信息。
3.  **profile 合并（multiz-like）**：
    *   对每个候选区间内的多个 block profile，执行 profile–profile 带状 DP（yama 直译，见 §2.12）：
        *   在参考坐标附近采用带状 DP（Radius R），解决不同 `.fas` 之间参考 gap 的不一致。
        *   根据 sum-of-pairs 打分决定保留哪些列/序列，以及如何插入额外 gap。
    *   输出合并后的单个 block（或少数几个 block）。
4.  **后处理与 refine**：
    *   输出的 `.fas` 可以再交给 `pgr fas refine` 做局部 MSA，以获得更"平滑"的 alignment（尤其是在非参考序列上）。

### 2.4 与现有 core 流程的互补关系

*   `p2m + join`：
    *   假设参考骨架在各数据源中基本一致。
    *   倾向于"只相信大家都同意的部分"（严格交集），适合构建 core genome。
*   fas-multiz：
    *   允许不同数据源在边缘和 gap pattern 上存在一定差异，通过 profile 合并策略尽量"合在一起"。
    *   输出更偏 union/mesh，适合探索 union pan-genome 或 WGA 风格的结果。

在实现层面，fas-multiz 作为一个独立子命令（`pgr fas multiz`），与 `p2m + join` 的适用场景不同：前者追求覆盖度（union），后者继续服务于一致性（intersection）。

### 2.5 命令行接口

*   子命令名称（示例）：
    *   `pgr fas multiz`
*   核心参数（示例）：
*   `-r, --ref <NAME>`：参考物种/序列名称，必须在所有输入 `.fas` 中存在。
*   `<infiles>...`：位置参数，输入的 block FA 文件，数量 `>= 2`，行为与 `pl p2m` 一致。
    *   `--radius <INT>`：带状 DP 半径 `R`，类似 multiz 中的 `R`，控制参考坐标附近的搜索宽度。
    *   `--min-width <INT>`：最小输出 block 宽度，对标 multiz 的 `M`。
    *   `-o, --out <FILE>`：输出 `.fas` 文件名。
    *   `--score-matrix <FILE>`：可选，指定替换矩阵（默认可复用 `libs/chain/sub_matrix.rs` 中已有配置）。
    *   `--mode <core|union>`：模式切换：
        *   `core`：在交集区域内行为尽量贴近 `p2m + join`，只对 gap 冲突做最小修复。
        *   `union`：尽量保留所有输入的对齐信息，生成 mesh 风格结果。

### 2.6 约束与实现注意事项

*   **参考骨架一致性**：
    *   要求所有输入 `.fas` 的参考序列来自同一基因组版本，且建议事先经过相同的 masking/裁剪流程。
*   **窗口化处理**：
    *   实现时应采用窗口化策略（例如按固定长度或按 block 切分），避免在超长区间上运行大规模 profile DP。
    *   每个窗口内的 profile 合并结果可以再交给 `fas refine` 做一次本地 MSA。
*   **打分与带状 DP**：
    *   可以重用 `pgr` 中现有的打分矩阵和 gap 参数（如 `libs/chain` 相关代码），避免在 `fas` 层重新定义一套 scoring。
    *   带状 DP 的半径 `R` 和最小宽度 `M` 建议和 multiz 保持同一数量级，以便结果直观可控。
*   **失败与降级策略**：
    *   当某个窗口内 profile DP 无法找到合理路径（打分过低或冲突过多）时，可以退回到简单的 `fas join` 行为，或干脆将该窗口标记为"未合并"，交给上游/下游流程决定如何处理。

### 2.7 与现有模块的集成点

*   **输入准备**：依赖现有的 `pgr axt/maf to-fas` 和 `fas` 系列命令，将所有上游结果规整为块级 `.fas`。
*   **区间计算**：复用 `fas cover` 和 `spanr` 的区间逻辑，定义候选合并窗口。
*   **比对与 refine**：在新实现的 fas-multiz 中完成 profile 合并后，调用现有 `pgr fas refine` 作为可选的精修步骤。
*   **下游分析**：输出 `.fas` 可以继续被 `fas stat`, `fas to-vcf`, `fas split` 等命令消费，与当前 `p2m + join` 的结果处于同一生态。

### 2.8 libs 实现概览

> 2026-02 更新：本节描述的是 `libs::fas_multiz` 在库层面的整体设计，与当前实现基本一致（包括 `FasMultizMode`/`FasMultizConfig`/`Window` 以及 `merge_window`、`merge_fas_files`、自动窗口推导等）。更细节的实现行为与局限见 2.10 小节。

本节给出 fas-multiz 在 Rust 中的 libs 级别设计。

*   **模块位置**：
    *   新增 `src/libs/fas_multiz/`，在 `src/libs/mod.rs` 中通过 `pub mod fas_multiz;` 暴露。
*   **依赖复用**：
    *   解析 `.fas`：复用 `libs::fmt::fas` 中的 `FasEntry`、`FasBlock`、`next_fas_block` 等。
    *   区间坐标：继续使用 `intspan::Range`。
    *   打分与碱基类型：复用 `libs::nt::NT_VAL` 以及 `libs::chain` 中已有的替换矩阵和 gap 参数。
    *   简单统计/评估：如有需要，可调用 `libs::alignment::alignment_stat` 做 sanity check。
*   **核心类型**：
    *   合并模式：
        ```rust
        pub enum FasMultizMode {
            Core,
            Union,
        }
        ```
    *   配置结构：
        ```rust
        pub struct FasMultizConfig {
            pub ref_name: String,
            pub radius: usize,
            pub min_width: usize,
            pub mode: FasMultizMode,
        }
        ```
    *   窗口定义：
        ```rust
        pub struct Window {
            pub chr: String,
            pub start: u64,
            pub end: u64,
        }
        ```
*   **对外 API 草图**：
    *   文件级合并（供 CLI 使用）：
        ```rust
        pub fn merge_fas_files(
            ref_name: &str,
            infiles: &[impl AsRef<Path>],
            windows: &[Window],
            cfg: &FasMultizConfig,
        ) -> anyhow::Result<Vec<FasBlock>>;
        ```
        *   读入多个 `.fas` 文件，根据给定窗口把 block 分组，对每个窗口调用 `merge_window`，最终返回按参考坐标排序的一组 `FasBlock`。
    *   单窗口合并（算法核心）：
        ```rust
        pub fn merge_window(
            ref_name: &str,
            window: &Window,
            blocks_per_input: &[Vec<FasBlock>],
            cfg: &FasMultizConfig,
        ) -> Option<FasBlock>;
        ```
        *   输入是某个窗口内来自多个文件的 block 集合，输出是一个合并后的 block（或在无法合理合并时返回 `None`）。
*   **merge_window 内部步骤概述**：
*   将每个输入中参考物种的 `FasEntry` 映射到统一的参考坐标网格上，得到多条略有差异的参考轨迹。
*   在参考轨迹之间执行带状 profile 对齐：DP 网格仍然只在参考坐标上展开，但每个对角单元的得分由两个 profile 的物种交集上的 sum-of-pairs 决定（共享物种的 base–base 使用 `libs::chain::SubMatrix::hoxd55` 的替换分数并做适当缩放，base–gap 使用统一的 gap 罚分，gap–gap 不计分），对两个及以上输入采用 progressive 带状 DP。
*   按照合并后的参考轨迹，对每个输入的非参考序列进行重采样：在缺失列处插入 gap，在 Union 模式下允许在参考 gap 位置引入新列，在 Core 模式下则尽量丢弃不一致列；非参考物种在 DP 打分中参与 sum-of-pairs，但坐标仍然沿参考轨迹重采样。
*   将重采样后的各物种序列按列拼接，构造新的 `FasBlock`，并为参考 entry 生成合适的 `Range`（可以取窗口的 Range 或交集 Range）。
*   如果在某个窗口内 profile 对齐得分过低或冲突过多，则返回 `None`，由调用者决定使用简单 `fas join` 还是跳过该窗口。

### 2.9 与 multiz-multiz 源码的异同

这里的 fas-multiz 是从 `multiz-multiz` 源码直译出来的 "pgr 版本"：DP 引擎逐条对应移植（§2.12），工作层级上做了调整（fas 层 vs MAF 层）。

*   **共同点（继承 multiz 的部分）**：
    *   都是以参考物种坐标为主轴，在参考坐标上定义窗口/段落，再在每个窗口内做 profile 合并。
    *   都采用带状 DP（或类似思想）限制搜索空间，在参考附近做局部优化，而不是在全空间做 MSA。
    *   在 union/mesh 场景下，都试图尽可能保留不同输入中的真实比对关系，只在必要时删除或压缩冲突列。
    *   支持"核心交集 + 扩展区域"的思路：核心部分倾向于各输入一致，边缘部分允许有差异并通过 DP 协调。
*   **差异（pgr 有意做的调整）**：
    *   **工作层级不同**：
        *   multiz-multiz 在 MAF 层操作，直接对齐两个 MAF profile。
        *   fas-multiz 在 block FA 层操作，输入是多个 `.fas` 文件，链路由 `pgr axt/maf to-fas` 标准化过，因此语法和元数据更简单。
    *   **DP 引擎复杂度不同**：
        *   multiz-multiz 的 `yama` 部分实现了一套完整的 profile–profile DP 引擎（sum-of-pairs + 替换矩阵 + 准自然 gap 模型 + LB/RB 参考锚定）。
        *   fas-multiz 的 `banded_align.rs` 是这套引擎的**直接移植**（C/D/I 三状态、GAP 构型查表、端部 gap 免费、LB/RB 锚定，见 §2.12），不在引擎复杂度上做简化；与 multiz 的差异只在工作层级（fas vs MAF）与整体合并次序（内容驱动贪心 vs guide tree）。
    *   **实现位置与职责边界不同**：
        *   multiz-multiz 是一个专门服务于 MAF/多序列比对构建的独立 C 项目。
        *   fas-multiz 被设计为 `pgr` 的一个 libs 模块（`libs::fas_multiz`），与现有 `fas cover/slice/join/refine` 等命令协作，而不是独立的 pipeline。
    *   **输入准备和预处理链路不同**：
        *   multiz-multiz 直接消费上游链路输出的 MAF（如 blastz/last 等）。
        *   在 pgr 中，上游的 pairwise 结果通常已经通过若干步骤转换、规范成 `.fas`，fas-multiz 可以假设这些输入已经做过一次清洗/规整。
    *   **目标偏好与使用场景不同**：
        *   multiz-multiz 更偏"通用 WGA 引擎"，追求在大范围基因组上做 mesh 式对齐。
        *   fas-multiz 明确被设计成 pgr 的一个"union/mesh complement"：在 core/intersection 流程之外，提供一个额外的 union 视角，并保持与现有 `p2m + join` 在交集区域内尽量兼容。

### 2.10 当前 fas-multiz 实现状态（2026-02）

> 本节描述的是当前 `pgr` 仓库中已经落地的 `libs::fas_multiz` 实现，可与前文对 multiz 及 fas-multiz 的设计描述对照阅读。单步合并的 DP 已升级为 yama 引擎直译（§2.11/§2.12），本节保留最初的工程描述作为实现记录。

**实现位置与对外 API**

*   模块位置：`src/libs/fas_multiz/`，通过 `pub mod fas_multiz;` 暴露为 `pgr::libs::fas_multiz`。
*   核心类型：与前文给出的设计保持一致，并在配置中加入了 DP 打分参数：
    *   `FasMultizMode { Core, Union }`
    *   `FasMultizConfig { ref_name, radius, min_width, mode, match_score, mismatch_score, gap_score }`
    *   `Window { chr, start, end }`
*   对外函数：
    *   `merge_window(ref_name, window, blocks_per_input, cfg) -> Option<FasBlock>`
    *   `merge_fas_files(ref_name, infiles, windows, cfg) -> Result<Vec<FasBlock>>`
    *   `merge_fas_files_auto_windows(ref_name, infiles, cfg) -> Result<Vec<FasBlock>>`

**窗口推导与 Core/Union 语义**

*   `merge_fas_files` 需要调用方显式给出 `windows`，行为与前文 2.3–2.8 小节给出的设计一致。
*   `merge_fas_files_auto_windows` 会：
    *   从所有输入 `.fas` 中提取参考物种 `ref_name` 的 `Range`，按 `radius` 向两侧扩展。
    *   按染色体合并重叠区间，再按 `min_width` 过滤过短窗口。
    *   按 `cfg.mode` 过滤窗口：
        *   `Core`：只保留"在所有输入中都有参考覆盖"的窗口（严格交集）。
        *   `Union`：只要有任意一个输入在该窗口有参考覆盖即可保留（并集风格）。

**窗口内合并逻辑（带状 DP 合并）**

*   一般情况（任意输入个数）：
*   对于给定窗口，先从每个输入文件中选出在窗口内与参考重叠的 block，组成 `blocks`。
*   若 `blocks` 为空，或（在 Core 模式下）某些输入找不到参考 block，则直接返回 `None`。
*   若 `blocks.len() >= 2`，先尝试 progressive 带状 DP 合并：
*   使用内部函数 `merge_blocks_with_dp`，按内容驱动的确定性顺序（§2.11）对 `blocks` 做两两 DP 合并。
*   每一步都要求参与合并的参考 entry 在去掉 `'-'` 后的序列完全相同（ungapped equal），否则这一轮 DP 失败。
*   在参考坐标网格上调用 `banded_align_refs`：
*   只在 diagonal ± `radius` 的带内做 DP。
*   对每个对角单元，按 multiz yama 语义打分（阶段 1.5 起，见 §2.12）：全体物种对（K×L 笛卡尔积）参与 sum-of-pairs，base–base 用替换矩阵（/50 缩放），base–gap 收一次 gap_extend，gap–gap 为 0；I/D 状态按"插入列非 dash 数 × 对方行数 × gap_extend"收费，横向/纵向移动按准自然 GAP 构型查表收 gap-open。
*   将 DP 生成的参考轨迹映射到所有物种：
*   对每一列，优先从前一个累积结果（或第一个输入）的对应位置取碱基，不存在时再从当前输入取；两边都缺失则填 `'-'`。
*   `Core` 模式下只合并在当前累积结果和新输入中都存在的物种；`Union` 模式下允许物种只存在于其中一边。
*   在 Core 模式下，任一步 DP 失败都会导致整个 progressive 合并失败，随后回退到"保守合并"逻辑。
*   在 Union 模式下，如果某一步 DP 失败，则跳过该输入，继续尝试将后续输入与当前累积结果进行 DP 合并；成功的部分会被保留，无法对齐的输入则在该窗口中被忽略。
*   progressive DP 完成后（无论是否跳过了一些输入），若至少完成了一次成功的 DP 合并，则直接返回最终累积的 block。
*   如果 progressive DP 入口阶段就失败（例如前两条参考轨迹 ungapped 不同，或带宽内找不到合理路径），则自动回退到"保守合并"逻辑：
*   要求所有候选 block 的参考 entry 完全相同（包含 gap），否则返回 `None`。
*   `Core` 模式下只保留在所有输入中都存在的物种；`Union` 模式下保留物种并集。
*   参考物种的 `Range` 继承自模板 block；其他物种继承其来自的原始 block。

**当前实现的局限与后续扩展方向**

*   多输入通过 progressive 两两合并完成，每步只看到"当前累积块 + 一个新输入"的物种集合，尚未实现一次决策看到全体输入的 SP-DP（见下方演进方向）；合并次序已由内容驱动的贪心策略确定化（§2.11），不再依赖输入文件顺序。
*   DP 网格是参考×参考的二维带状网格，与 multiz yama 的成对 profile 合并一致；非参考物种通过 sum-of-pairs 打分参与决策，但不拥有独立的坐标轴。
*   替换分数已经复用 `libs::chain::SubMatrix` 做 base–base 的 sum-of-pairs 打分：默认使用 `hoxd55`，也支持通过 `--score-matrix` 读取 LASTZ 格式文件或预设名称（例如 `hoxd55`），并通过简单缩放与当前 `match_score` 的量级对齐。gap 支持三类模型：`constant`、`medium`/`loose`、以及显式仿射：
    *   `constant`：直接使用 `gap_score` 作为统一线性 gap 罚分。
    *   `medium`/`loose`：从 `GapCalc::medium`/`GapCalc::loose` 的 quasi-natural 曲线中取 `len=1,2` 两点，反推出一组近似的仿射参数 `(open, extend)`，再按 HoxD55 的打分尺度和 `match_score` 做线性缩放，在带状 DP 中用"open + extend × length"的形式累积 gap 罚分，从而实现长度依赖的 quasi-natural 近似。
    *   显式仿射：当通过 `--gap-open`/`--gap-extend` 提供 open/extend 时，fas-multiz 在 DP 中直接使用这一组仿射参数（同样按 `match_score` 缩放）进行三状态的仿射 gap 计分。
*   已提供 CLI 子命令 `pgr fas multiz`，支持 `--mode core|union`、`--radius`、`--min-width`、`--gap-model`、`--gap-open`、`--gap-extend` 以及 `--score-matrix` 等参数；gap 配置风格与 `pgr psl chain` 保持一致，而替换矩阵也不再局限于内置的 HoxD55，可与链化阶段共享同一套 matrix 配置；`libs::fas_multiz` 仍作为底层引擎，便于在 pipeline 或其他子命令中复用。
*   在 gap 行为上，端部 gap 与 multiz 语义一致视为免费（不收 gap-open，extend 照收），由 DP 状态机直接处理（§2.12）；回溯后只裁剪两端"全物种都是 gap"的列，不裁剪单侧 gap 列（旧规则会误删真实内容）。
*   在上述基础上，仍可以在后续逐步接近 multiz 的完整行为，例如：
    *   将 progressive DP 升级为真正的多输入 profile–profile sum-of-pairs 动态规划。这里的"真正"并不是指在 K 条序列上做天真的 K 维 DP（那样复杂度是 O(L^K)，在 K 稍大时不可用），而是指在工程上尽量在同一个 DP 决策里综合所有输入对的 sum-of-pairs 打分，减少合并顺序对结果的影响。一个可行的演进路径可以分为三个阶段：
        1. **全体物种 SP 打分 + 合理的合并次序**：~~在仍保持当前"以参考为一维"的带状 DP 框架下，把 scoring 从"只看当前两个 block 的共同物种"升级为"在一个窗口内对全体物种做 sum-of-pairs（缺失视为 gap）"，并配合更合理的 merge 次序（例如基于 guide tree 或其他拓扑），这样即便 DP 依然是 pairwise reference–reference，决策时看到的是"全局 profile"的得分，progressive 的顺序敏感会明显减弱。~~ **✅ 已实现（2026-08-03，§2.11）**：共享物种自对 + 参考交叉对评分 + 内容驱动的确定性合并次序；输入顺序无关已有回归测试。
        2. **局部的小 K 多输入 DP（K≤3）**：在窗口长度较短、冲突较集中的局部，引入一个只针对少数输入（例如 2–3 条参考轨迹）的 exact 多输入 DP 分支；在这一分支里，状态是 (i,j,k...) 这样的多维索引，每一个 DP 列都按"全体物种的 SP 打分"计分，用来精修最困难的区域，而大多数区域仍然走带状 2D DP 路径，从而在不爆炸复杂度的前提下局部地"接近理想解"。
        3. **多轨迹但仍是 2D 网格的近似多输入 DP**：在充分掌握前两阶段行为的基础上，可以尝试构建"参考坐标 × 合法状态"的 2D DP 网格：横轴仍然是参考（或参考对）的坐标，纵轴是有限集合的"多物种轨迹状态"（例如用 bitmask 或离散状态表描述某一列中哪些物种前进一步、哪些物种打 gap），通过严格限制合法状态和转移（如必须顺着预先给定的链/轨迹前进，禁止任意插入/删除）来控制状态空间大小，在每个 DP 列上仍按全体物种的 sum-of-pairs 打分。这一层相当于在现有 fas-multiz 参考框架上实现一个工程化的"多输入 SP-DP 近似版"，在不引入完整 K 维 DP 的前提下向 multiz 的行为靠拢。
*   在 DP 失败时更智能地选择降级策略（退回 `fas join`、标记窗口未合并等）。

### 2.11 阶段 1 落地记录（2026-08-03）：全体物种 SP 打分 + 确定性合并次序

实现演进路径第 1 阶段，两处改动：

**1. 对角评分的物种集合（`banded_align.rs`）**

从"只对两个 block 的共同物种打分"升级为"共享物种自对 + 参考交叉对"：
每个物种的列都与参考锚点物种成对计分（星型拓扑），加上共享物种的自对。

实现过程中的两个实证发现（S288c 三输入真实数据，参考序列去 gap 后 3826 bp）：

*   **朴素全笛卡尔积（所有跨物种对）会拖偏参考锚点**：直接对全体物种两两
    计分后，DP 在真实数据上产生 155 列额外错位（合并参考去 gap 后 3981 bp
    而非 3826）。原因是物种携带真实 indel 时，对角线上"碱基-vs-gap"的
    逐对 gap 罚分（~-80）远大于逐对错配罚分（~-2），DP 倾向移位把碱基配
    对起来而避开 gap 列。共享物种自对 + 参考交叉对（不含物种间交叉对）
    消除了该问题。
*   **对角线上的 gap 贡献改为 0**：gap 成本完全由仿射 gap 转移
    （`gap_i`/`gap_j`）承担。逐对再收一次 gap 罚分会与替换分数失衡，
    同样诱发移位。

**2. 确定性合并次序（`merge.rs`）**

progressive 合并不再按输入文件顺序，而是按块内容做贪心聚合：先选物种数
最多的块，之后反复选择与已累计物种集重叠最大的块；并列时按参考区间 +
物种名排序的内容键打破。该顺序只依赖块集合本身，与输入顺序无关。

配套修正：合并后的参考序列只保留第一块的参考序列（不再用第二块碱基填充
第一块的 gap）。旧行为会把去 gap 后的参考膨胀 12+ bp，破坏下游合并依赖的
`ungapped_equal` 不变量；修正后合并参考恒等于输入参考（3826 bp）。

回归测试：`merge_window_output_independent_of_input_order`（同一组 3 个块
的 6 种输入排列产出完全一致的合并块：物种名 + 序列）。

**范围边界**：物种间交叉对（YJM-vs-Spar 这类）暂不进 DP 打分；若未来要
引入，需要先解决"逐对 gap 罚分与替换分数的量级校准"问题（阶段 1.5），
否则会重新拖偏参考锚点。

### 2.12 阶段 1.5 落地记录（2026-08-03）：yama 引擎直译

对照 multiz 源码（`multiz-multiz/mz_yama.c` / `mz_preyama.c` /
`mz_scores.c`）把 `banded_align.rs` 的 DP 升级为 yama 的完整语义：

**1. C/D/I 三状态 + 准自然 gap（`mz_yama` 直译）**

*   C（替换）/ D（删除，A 列插入 B 全 dash）/ I（插入，B 列插入 A 全 dash）
    三状态，每格记录 `flag_c | flag_d<<2 | flag_i<<4` 供回溯；
*   gap 成本用准自然 GAP(s,t,u,v) 查表：16 种"最后两条边"构型中 6 种收
    gap_open，其余 0（`mz_scores.c` 的构型定义逐条直译）；
*   对角列分 SS 为**全体物种对（K×L 笛卡尔积）**：base-base 用替换矩阵
    （/50 缩放），base-gap 收一次 gap_extend（multiz `SS('-',x)=-extend`），
    gap-gap 为 0——解决了 §2.11 记录的量级校准问题（multiz 的 per-pair
    gap 成本是 extend 而非 open+extend，不会拖偏锚点）；
*   I/D 状态按"插入列非 dash 数 × 对方行数 × gap_extend"收费；
*   **端部 gap 免费**：I 在末行、C 在起点列、D 在起点/终点列不收 gap-open
    （extend 照收），与 multiz "End-gaps are not charged a gap-open penalty"
    一致。

**2. 参考锚定 LB/RB（`mz_preyama` 直译）**

*   把两个参考序列去 gap 后逐位配对（第 k 个碱基 ↔ 第 k 个碱基），为每个
    ref_a 碱基列建立 (lb, rb) 点约束；
*   `smooth`：先单调化（lb 前向、rb 后向），再做半径"香肠"扩展；
*   DP 每行只在 [lb[i], rb[i]] 内计算，格点按行存储（变宽带，非固定对称
    band）；
*   **必要性实证**：无 LB/RB 时自由端 gap 会让 DP 把列数差（如 336 列）整段
    堆到块端，随后被边缘裁剪删掉真实内容（Spar 4057→3724）；LB/RB 把路径
    锚在参考对角线上，端部 overhang 被限制在 radius 内。

**3. 边缘裁剪改为"全 gap 列"判定**

旧裁剪删除"单侧列"（某侧 map 为 None），在自由端 gap 语义下会误删真实内容
（ref 末尾碱基、边缘物种碱基）。现改为只删**两端全物种都是 gap 的列**：
用两块的物种列剖面（`col_a`/`col_b`）判定该列是否真的空。

**验证**：S288c 三输入 Union 合并后 4193 列，参考去 gap 恒 3826，
RM11_1a 3834 / Spar 4057 / YJM789 3822 —— 与输入逐碱基一致，零丢失；
`merge_window_preserves_species_content`（微型内容守恒测试）与全部既有
测试通过。DP 的 gap 模型、LB/RB、端部处理均为 multiz 语义的直译，与
`multiz-multiz/` 源码逐条对应。

**P2 状态（2026-08-03）**：
*   **hox70 矩阵预设 ✅**：`SubMatrix::from_name` 增加 `hox70` 别名——pgr 的
    `hoxd55` 与 multiz 的 HOX70 数值完全相同（91/-114/-31/-123，gap
    open 400 / extend 30），只是命名不同；`--score-scheme hox70` 现可显式选择。
*   **v=0 模式（未做）**：两个参考都可调需要第二次 yama 对齐参考行；pgr 的
    渐进合并以累计块参考为锚（等价 multiz v=1），v=0 对窗口模型的价值不明，
    不做推测性实现。
*   **multiz 对齐回归（部分受阻）**：两工具输入格式（MAF vs fas）与分块语义
    （multiz 逐重叠区分块流 vs pgr 单窗口块）不同，字节级对比不可行。实测
    MAF→fas 转换后喂 pgr 时发现一个**预先存在的输入模型限制**：多块输入中
    部分重叠的块若只共享参考且参考去 gap 不等（`merge_conflicting_refs` 需要
    共享非参考物种打分 crossover），合并被拒绝。这是 pairwise 切片输入模型
    的边界，主流程（`tests/fas/*.slice.fas` 切片输入）不受影响；修复需要
    重新设计"仅参考共享"的 crossover 打分，留待真实需求出现。
