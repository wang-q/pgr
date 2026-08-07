# HV（Hypervector）设计笔记

> 2026-08-07 建立，持续更新。本文顺序：算法概览（§1）→ 影响算法设计的
> 因素（§2，核心决策记录）→ pgr 落地现状（§3）→ HyperGen 论文与代码
> 参考（§4，外部参考：HyperGen + hdlib，背景材料）→ 后续方向（§5）。
> 来源：外部参考详细分析见 [[../references/hv.md]]（HyperGen 论文与代码、
> hdlib；该文件将随后续文献持续扩充）。
> * 基准与验证：[[../benchmarks/bench-simd-hv-jaccard.md]]、
>   [[../benchmarks/dist-cohort-validation.md]]。

## 1. 算法概览

### 1.1 目标

把每个基因组的 k-mer 集合（通常几十万到几百万个）投影成一个**固定维度 D
的高维整数向量（HV）**，比较两个基因组 = 比较两个 HV（O(D)），用于大规模
集合的近邻/粗筛/聚类。集合大小与相似度直接从向量估计，不需要逐 k-mer
比对。

### 1.2 流水线

1. **k-mer 采样**：minimizer / closed syncmer（FASTA 侧）或 `.pgi` 的
   unique k-mers（索引侧），得到 k-mer hash 集合（大小 N）。
2. **随机投影（编码）**：每个 k-mer 的 hash 值作为伪随机数生成器（RNG）
   的种子，生成 D 维随机向量并聚合成 HV。两种稠密编码路径：
   * **bit（±1）路径**：每维 ±1；
   * **i8 路径**：每维一个随机字节（−128..127）；
   每个 k-mer 影响**所有** D 维（另有稀疏路径只碰 s 个维度，见 §3.1）。
3. **距离估计**（利用准正交性）：
   ```
   card  = ‖H‖² / D                     （集合大小）
   inter = H_A · H_B / D                （交集）
   J     = inter / (card₁ + card₂ − inter)  （Jaccard）
   mash  = −(1/k)·ln(2J/(1+J))          （Mash 距离；ANI = 1 − mash）
   ```

### 1.3 设计因素一览

算法骨架（采样 → 投影 → 距离）是稳定的，真正的实现自由度在"聚合那一步
怎么算得快、算得准"，即 §2 要展开的因素：

* **编码路径**（bit vs i8）：决定统计性质与每多少维消耗一次 RNG（§2.1）；
* **随机数生成速度**：随机位从哪里来、生成成本是否成为瓶颈（§2.2）；
* **接收容器宽度**（累加器类型）：聚合结果用什么类型接收（§2.3）；
* **位拆分/展开**：每 bit 的展开指令数，当前限速点（§2.4）；
* **幅度**：不同编码下每维数值的分布范围，及其对区分度的（无）影响
  （§2.5）。
* **采样方法**（FracMinHash vs minimizer vs syncmer）：决定被保留的
  k-mer 集合是否是"近似均匀随机子采样"，即 Jaccard/ANI 估计是否存在
  系统性偏差——这是比编码路径更早的一层偏差来源（§2.6）。

### 1.4 k-mer 为什么不能直接当种子

k-mer 本身不是哈希，而是 **2-bit 打包的序列整数**，结构性强：相邻 k-mer
共享 k−1 个碱基（数值几乎相同），低复杂度区（poly-A、串联重复）有大量
近重复 k-mer。流水线里的"种子"是**采样器哈希之后的值**（minimizer 用
rapid/fx/murmur，syncmer DNA 用乘性滚动哈希），再喂给
`RapidRng::seed_from_u64(hash)`——所以"拿 k-mer 当种子"与"拿 k-mer 哈希
当种子"在当前实现中是同一件事。

若跳过哈希、直接用原始 2-bit 编码当种子：

* 相似 k-mer → 相关随机向量 → 破坏准正交性，聚合时重复区把某些维度
  系统性拉偏，`‖H‖² ≈ N·D`、`dot ≈ shared·D` 的统计假设失效；
* canonical 语义要求正反链映射到同一值，原始编码做不到；
* 哈希把 k-mer 压到近似均匀的 u64，种子流独立性更好。

特例：`pgi to-hv` 的 `key_to_seed` 直接 XOR 折叠 2-bit key 当种子，没有
单独哈希——但 `hash_hv_sparse` 内部每个 (seed, dim) 都过 splitmix64
（强混合）兜底。**关键不是哈希放哪一步，而是种子里不能残留可预测结构**：
FASTA 路径把混合放在采样器，pgi 路径藏在 splitmix64，HyperGen 用 t1ha2。

**为什么"相似 k-mer → 相似向量"不是好事**：本方案的相似度语义是**精确
k-mer 集合重叠**（Jaccard → ANI 可换算）。相似但不相同的 k-mer 若也贡献
相关信号，"近但不等于"会被计入点积——低复杂度/重复区会让高估量取决于
重复与突变结构（同一基因组对谁都显得更近），Jaccard 无法无偏换算回 ANI。
想让近匹配被捕捉，正确的位置是**采样层**（k 大小、syncmer 参数、spaced
seeds），而不是让投影变模糊；准正交投影保持"精确计数"，敏感度由采样决定。
模糊相似度（SimHash 式排序）是另一类任务，与本方案的 ANI 估计语义不兼容。

**模糊匹配的适用场景（补充）**：对较远的近缘（如 95% ANI），精确 k-mer
重叠快速衰减（单碱基突变率 d 下 21-mer 精确存活率 (1−d)^21 ≈ 0.34），
"相似 k-mer 也计一点"确实能救回灵敏度——这正是 spaced seeds / 更小 k /
错配容忍种子的动机。但要实现错配容忍，应做**确定性分桶**（spaced seeds、
k-mer 及其 1-错配变体归入同一桶），保持"桶集合计数"语义并重新标定距离
公式；若坚持在投影层做 LSH 式模糊，等于明确放弃 ANI 数值语义、接受重复
敏感，仅适合"粗筛排序"这类不要求无偏数值的任务。

## 2. 影响算法设计的因素

这些是我们在实现中实际测量/推导出的因素，每一项都给出了对实现的决策。
所有性能数字默认 n=10k、D=4096、AVX2 主路径（Ryzen 9 7945HX），除非注明。

> **数字标注说明**：§2.2 的 2.68/4.36 ms 与 §2.4 的 1.41/1.37/0.71 ms 是
> **历史/中间测量**（旧 i8 实现、延迟 −N 优化前的逐组件分解）；当前主路径
> 数字汇总见 §3.3（2026-08-08 复测一致，见
> [[../benchmarks/bench-simd-hv-jaccard.md]]）。

### 2.1 编码路径：比特（±1）vs i8

**统计语义**

* **比特路径**：每个 k-mer 种子对每维贡献一个随机位 → ±1。每维
  H = 2Σb − N，值域 [−N, N]，围绕 0 **上下平衡**（std ~√N）；准正交性
  好，是"干净"的编码。
* **i8 路径**：每个种子对每维贡献一个随机字节（0..255 → −128..127）。
  字节均值 −0.5 → 每维 H 的期望 ≈ −N/2，**不上下平衡**（直流偏置），
  std ~73√N。这个偏置让点积出现 ~N²/4 的二次噪声底，距离语义被破坏
  （cohort 实测与身份率 Spearman ≈ 0，见 [[../benchmarks/dist-cohort-validation.md]] §2）。

**速度**（AVX2，n=10k、D=4096）

* RNG 调用：bit 每 64 维 1 次（一次 `rnd_at` 用满 64 位：低 32 位 → 前
  32 维、高 32 位 → 后 32 维）；i8 每 8 维 1 次 → bit 的 RNG 成本约为
  i8 的 1/8。
* 实测：bit **1.11 ms**、i8 2.11 ms（bit 快 ~1.9×）。

**三组实现对比**（n=10k、D=4096；标量 = 每 seed 串行 RNG 流 + 逐位累加，
数据见 [[../benchmarks/bench-simd-hv-jaccard.md]]）：

| 编码 | AVX2 手写 | wide（可移植路径） | 标量（串行流） | AVX2/wide | AVX2/标量 |
|---|---:|---:|---:|---:|---:|
| bit（±1） | **1.11 ms** | 6.64 ms | 9.01 ms | ~6.0× | ~8.1× |
| i8 | **2.11 ms** | 4.38 ms | 11.69 ms | ~2.1× | ~5.5× |

差异根源（**三者 lane 数相同——都是 8（256-bit i32），差距不在 lane
宽度**）：

* **wide vs 标量只有 1.4–2.7×**（bit 6.64/9.01、i8 4.38/11.69）：wide 的
  8-lane SIMD 只加速"位提取/字节转换 + 累加"部分，RNG 仍是串行
  `next_u64`（i8 每 8 维 1 次，依赖链主导，见 §2.2 的 rng-only 2.68 ms
  vs 总 4.36 ms）——这正是"8 lane 只有约两倍"观察的来源。
* **AVX2 vs wide**（bit 6.0×、i8 2.1×）：AVX2 额外带来**跳步 RNG +
  块主序**（RNG 调用大幅降低、跨 seed ILP）+ 手写展开（无 wide 的冗余
  装载/转换/每次循环重建 shift 数组）；i8 因 RNG 频率高（每 8 维 1 次），
  跳步收益被稀释。
* **AVX2 vs 标量**：bit ~8.1×、i8 ~5.5×（两块收益叠加）。

**决策**：bit（±1）为主实现；i8 仅作为"保语义"变体保留（FASTA 路径当前
仍走 i8，其量纲问题见 §3.4、修复计划见 §5.3）。

### 2.2 随机数生成速度

**瓶颈定位**：旧实现（i8 seed-major）是 RNG 主导——rng-only 2.68 ms vs
总 4.36 ms；串行 `RapidRng::next_u64` 有依赖链，且每 8 维就要一次
128-bit mix。

**关键性质**：`RapidRng` 的状态是常数步长计数器
（输出 j = mix(seed + j·SECRET0, …)），因此可以**跳步**（`rnd_at`）。
这让循环可以改成**块主序**：HV 分块常驻寄存器、遍历全部种子，每个种子
每 64 维只做 1 次 RNG（一次 u64 输出用满），且不同种子的 mix 相互独立
→ 指令级并行。

**现状**：bit 主路径的 RNG 独立成本只有 ~0.3–0.7 ms，且与 SIMD 展开在
**不同执行端口上重叠**（标量 128-bit mix 走标量 ALU，展开走向量端口），
不是限速点。

**决策**：主路径用 RapidRng。标量对照（RapidRng / SmallRng(xoshiro256++) /
StdRng(ChaCha12)，i16 累加、n=10k）三者差距很小（9.15 / 9.38 / 10.45 ms，
见 [[../benchmarks/bench-simd-hv-jaccard.md]] §2）——选 RapidRng **不是因为
标量最快**，而是它的常数步长计数器可 O(1) 跳步（`rnd_at`），这是块主序
循环的前提；SmallRng 虽快但走 `next_u64` 串行链、无跳步接口，拿不到块主序
收益。换更快的 RNG（xorshift/LCG 等）收益有限，不投入。

**RNG 候选实测（2026-08-08，AVX2 bit 主路径 64 位框架，n=10k/D=4096）**：
RNG 输出统一加 black_box（阻止 LICM 广播折叠）后，wyrand / 宏版 rapid /
常量 RNG 三者持平（1.09–1.13 ms vs 主实现 1.1145 ms，差异 <3%），
splitmix64 慢 ~4%；常量（RNG 免费）只快 ~2%——**广播 + 依赖链占 ~98%
成本，mix 计算本身只占 ~2%**。经典 RNG 全部落选：MT19937 标量慢 ~2.9×
（每 seed 需初始化 624 项状态数组 + temper），LCG/PCG 的 O(log j) 跳步
在块主序内层每 chunk 重算、慢 2.7–9.6×——"O(1) 跳步（counter+mix）"
是块主序框架的硬前提。真实候选均无法超越 RapidRng，换 RNG 无收益，
维持原决策；进一步优化方向是削减每 seed 广播/依赖链而非换 RNG（基准
变体 `bit_avx2_rng_*` 保留在 `benches/hv_benchmark.rs`，数据详见
[[../benchmarks/bench-simd-hv-jaccard.md]] RNG 候选对比）。

> 采样哈希侧（HyperGen 的 t1ha2）同样不引入：21-mer 吞吐比 pgr 现役
> rapidhash 慢 ~1.76×（rapid 16.9 µs vs t1ha2 29.7 µs / 10k）。

**RNG 性质与 HV 编码的适配**（2026-08-08，外部交叉验证见 rapidrand 基准
与 [[../benchmarks/bench-simd-hv-jaccard.md]]）：

* **原生输出宽度**：64 位原生（RapidRng / wyrand / splitmix64 /
  xoshiro256++ / PCG64）一次 mix 给完整 64 位；32 位原生（MT19937 /
  PCG32 / ChaCha12 字 / MINSTD 31 位）的 u64 是两个 u32 拼的——成本 ×2
  （rapidrand 表：Pcg32 u64≈2×u32、ChaCha12 u64≈2×u32）。"生成 u32 是
  普遍行为"只对经典款成立，现代快速款原生 64 位。
* **跳步能力三档**：O(1) counter+mix（RapidRng / wyrand / splitmix64）
  ——块主序的硬前提；O(log j)（LCG / PCG，每 chunk 重算，实测慢
  2.7–9.6×）；不可跳步（MT19937，标量慢 ~2.9×）。
* **单次 draw 速度档**（rapidrand，M1 Max）：快速款 ~0.5 ns（RapidRng /
  WyRand）、1–2 ns（xoshiro / PCG）、4–13 ns（ChaCha12/8/20）。
* **适配结论**：HV 块主序需要"64 位原生 + O(1) 跳步 + 快"三者齐备，
  缺一即慢或不可用；RapidRng 与 wyrand 类是最佳候选，两者实测持平。

> **边界澄清**：`hash_hv_sparse` 的 splitmix64 是投影维度的确定性派生
> （非 RNG 实验）；HyperGen 参考实现的 WyRng/t1ha2 仅见于
> [[../references/hv.md]]，不参与 pgr 的选择。

### 2.3 接收容器宽度（累加器类型）

* **i32**：无条件安全（任意 N、任意编码），当前所有路径采用。
* **i16 + 分段**（bit 路径候选）：每段 S 个种子，段内每维值域确定地落在
  [−S, S]，S ≤ 32767 就**确定不溢出**（不依赖概率）。实测反而慢
  **~2.5×**（3.50 vs 当时 i32 基线 1.40 ms）：16-bit 变量移位 `vpsrlvw` 延迟更高，
  且 32 维块从 4 条独立累加链降到 2 条——这循环是链式延迟主导，
  指令数少不顶用。**2026-08-08 64 位框架复测**（`hash_hv_bit_i16`，每 64
  维 4 个 i16 16-lane 寄存器，n=10k 不分段，值域 [−10k, 10k] 安全）：
  3.64 ms vs i32 主实现 1.11 ms，**慢 ~3.3×**——链 8→4 条 + `vpsrlvw`
  延迟 + 每块末尾 i16→i32 扩展转存（i32 路径无此步）；结论不变，且
  随链数增加差距略扩大。
* **延迟 −N 偏移**（bit 路径采用）：既然 ±1 数值围绕 0 平衡，就不必
  每种子做 `2b−1`；改为每种子只累加 `2b`、每块末尾一次性减 N——每组
  少 1 个向量 op，4 条累加链不变 → **1.40 → 1.14 ms**（n=100k 14.0 →
  11.4 ms；2026-08-08 64 位用满后再降至 1.11 ms / 11.2 ms，见 §3.3）。
* **i8 必须 i32**：直流偏置 + 值域 ±128N；260 万种子实测各维 ±1.3e6，
  i16 必溢出。

**决策**：bit 与 i8 都用 i32；bit 用延迟 −N 偏移吃下"数值平衡"的红利；
i16 lane 分段累加不采用。

### 2.4 位拆分/展开（当前限速点）

逐组件分解（deferred 前）：完整 1.41 ms，纯展开 **1.37 ms**（~97%），
纯 RNG 0.71 ms（与展开重叠，被掩盖）。采用延迟 −N 后完整降到 1.14 ms
（2026-08-08 64 位用满后再降至 1.11 ms），展开仍是主体。

每 bit 至少需要 ~5 个向量 op（variable shift → and → shift-left → add，
加上 broadcast 摊销）；i16 lane 不能减少每 bit 的 op 数（实测更慢，
见 §2.3）。

**决策**：进一步提速应从"减少每 bit 指令数"入手，而不是换 RNG 或换
lane 宽度。AVX-512 参考实现（旧式每步 sub，已同步 64 位用满）Zen 4 实测
1.44 ms，反而比当前 AVX2 主路径（1.11 ms）慢，作者决策不采用 AVX-512
（见 [[../benchmarks/bench-simd-hv-jaccard.md]] §2）。

### 2.5 幅度与区分度

**观察**：i8 每维幅度远大于 bit（std ~73√N vs √N，n=10k 时 ±~7k vs ±~100），
初看似乎"同样 4096 维下不同基因组差距更大"。

**实测结论（数值模拟，N=10k、D=4096，共享 2000/500/50）**：

| 编码 | 每维 std | Jaccard 估计误差 |
|---|---:|---:|
| bit ±1 | ~100 | 0.013 / 0.013 / 0.005 |
| i8 真实现（均值 −0.5） | ~7300 | 0.16 / 0.17 / 0.19（基本报废） |
| i8 零均值（±127） | ~12800 | 0.013 / 0.013 / 0.005（与 bit 相同） |

原因：Jaccard 是**比值度量**，幅度是分子分母的公共因子，信号与噪声同时
缩放、在比值中抵消（零均值 i8 与 bit 的误差逐位相同）；真正的 i8 幅度
大是直流偏置的副作用，偏置反而引入二次噪声底，**更糟**。

**决策**：区分度由共享种子占比与维度 D 决定，与编码幅度无关。想要
4096 维下拉大基因组间差距 → 提高 D 或改进采样（FracMinHash / syncmer
参数），而不是换回 i8。

### 2.6 采样方法：FracMinHash vs closed syncmer 对估计的影响

**前提**：§1.2 的距离公式（card / inter / J → ANI）成立的前提是"被保留
的 k-mer 是近似均匀随机子采样"——每个 k-mer 独立地以同一概率被保留，
保留比例 ≈ 交集比例，J 才可无偏地换算回 ANI。采样方法决定这个前提是否
成立，是比编码路径（§2.1）更早的一层偏差来源；换 bit/i8 救不回采样偏差。

**三种采样器在文献中的定位**：

* **FracMinHash**（HyperGen 使用）：保留 `h(x) < 阈值` 的 canonical
  k-mer（Irber 等 2022 提出，bioRxiv 2022.01.11.475838）。每个 k-mer 的
  保留概率相同且独立 → 保留集合近似均匀随机子集。Hera 等 2023（Genome
  Res 33(7):1061–1068）证明其 Jaccard / containment 估计的偏差很小且
  **可精确校正**（给出校正公式），并据此为突变率（即 ANI 语义）推导
  **点估计 + 置信区间**（论文演示参数 k=21、scale 0.1）。
* **minimizer**（pgr FASTA 路径的另一选项）：窗口内保留最小 hash。
  保留概率与 k-mer 在窗口中的哈希排名 / 局部序列结构相关，不是均匀子
  采样。Belbasi 等 2022（ISMB, *Bioinformatics* 38:i169–i176）证明
  minimizer 的 Jaccard 估计**有偏且不一致**——偏差不随序列长度增长
  消失，估计收敛到与真实 J 不同的值。MashMap 旧版（minimizer 版）的
  ANI 估计即受此影响，Kille 等 2023 的 minmer 方案（*Bioinformatics*
  39(9):btad512）专门为修正它而提出。
* **closed syncmer**（pgr 主力）：窗口内最小 s-mer 落在首/末位
  （Edgar 2021, *PeerJ* 9:e10805；syng 移植，见
  [[../references/syng.md]]）。文献支持分两层：
  * **open syncmer 与 FracMinHash 等价**（Liu & Koslicki 2023,
    bioRxiv 2023.11.09.566463）：在 k-mer 集合相似度（Jaccard /
    containment）意义上等价，且距离分布与保守性更好。关键假设：
    syncmer 的采样对**重叠 k-mer 有依赖**（共享子串影响被选概率），
    FracMinHash 则完全独立，等价性正由此论证；
  * **closed syncmer 的无偏性有同行评审直接证据**：Shibuya 等 2022
    （WABI, LIPIcs 242:14，doi:10.4230/LIPIcs.WABI.2022.14）用与 pgr
    同规则的 closed syncmer 做集合比对，称其为 minimizer 的
    "上下文无关"替代、给出无偏 Jaccard 估计；作者博士论文实验里
    minimizer 明显有偏，syncmer 与随机采样在无偏性上重合。注意：
    closed syncmer 没有 FracMinHash 那种可写出的偏差量级、校正与
    置信区间公式，"无偏"停留在论证 + 实验层面。

**项目实测佐证**（[[../benchmarks/dist-cohort-validation.md]]）：

* `dist seq`（k=8 closed syncmer 草图层）是各草图层里与身份率最贴近的
  （ρ=0.616–0.816）——closed syncmer 的排序语义可用；
* `dist pgi`（k=40 syncmer 集合）Spearman 0.59，**确定但有偏**（采样
  位置漂移）——closed syncmer 的偏差真实存在，量级可控但不为零。

**结论与建议**：

1. 排序 / 粗筛任务：closed syncmer 够用（cohort 验证），保留现状；
2. 需要**可标定、可带置信区间**的 ANI：优先 FracMinHash（或 open
   syncmer），Hera 2023 提供校正与 CI 公式，HyperGen 已验证端到端
   （见 [[../references/hv.md]] §1）；
3. minimizer 仅作为历史兼容 / 对比基线，不应用于需要数值精度的新路径；
4. 与 §1.4 / §2.5 呼应：敏感度与偏差由**采样层**决定，投影层只负责
   精确计数——换编码路径救不回采样偏差。

## 3. 实现现状（pgr 落地）

> **当前推荐配置速查**：编码用 bit（±1）+ i32 累加 + AVX2 跳步 RNG
> （块主序）、D=4096；FASTA 侧采样用 closed syncmer（`dist seq` 默认
> k=8/w=5）；需要可标定 ANI 时优先 FracMinHash（§2.6）。⚠️ FASTA
> `dist hv` 路径（`load_hv_from_fasta` / `load_hv_from_fasta_syncmer`）
> 仍走 i8，量纲问题未修（§3.4）。

### 3.1 `src/libs/hv.rs` 当前实现与分派

| 函数 | 编码 | 说明 |
|---|---|---|
| `hash_hv_bit` | 稠密位编码 ±1，i32 | AVX2 主路径：跳步 RNG + 块主序 + 延迟 −N；每 64 维 1 次 RNG（用满 64 位，低 32 位→前 32 维、高 32 位→后 32 维） |
| `hash_hv_i8` | 稠密 i8 累加，i32 | AVX2 路径：每 8 维 1 次 RNG；保留原语义（含直流偏置） |
| `hash_hv_sparse` | 稀疏 ±1，i32 | splitmix64 派生，每 k-mer 只更新 `s`（默认 3）个随机维度；`.hv` v2 采用 |
| `hv_norm_l2_sq` / `hv_cardinality` / `hv_dot` | — | wide SIMD 范数；cardinality=‖H‖²/D；dot 按 √D 归一 |
| `calc_distances` | — | jaccard / containment / mash 多口径输出 |
| `load_hv_from_fasta(_syncmer)` | i8 | FASTA → minimizer / closed syncmer 集合 → `hash_hv_i8` |

分派链（x86-64）：`avx2` 运行时检测 → AVX2 intrinsics；否则 wide 可移植
路径（aarch64 自动用 NEON，其余标量）。所有路径输出**逐位一致**（含
非 32 倍数维度的尾部处理），由 `test_hash_hv_bit_serial_vs_simd` /
`test_hash_hv_i8_serial_vs_simd` 及 avx2 变体保证。
AVX-512 实现只保留在 `benches/hv_benchmark.rs` 作参考对照，不参与分派。
`cargo check --target aarch64-apple-darwin --lib` 交叉编译通过。

### 3.2 消费链与 `.hv` 格式

* `pgr dist hv`：FASTA 路径（minimizer/syncmer → `hash_hv_i8` →
  `calc_distances`）；`.hv` 路径（`pgi to-hv` 产物直接比较，稀疏余弦）。
* `pgr pgi to-hv`：把 `.pgi` 的 unique k-mer keys 投影成**稀疏** HV；
  `.hv` v2 格式：`PGV1` magic + version 2 + `k/dim/sparse/n_kmer/name` +
  i32 数组。稀疏投影 + 存储真实 n_kmer 是 v2 的关键修复。
* `pgr dist hv a.hv b.hv`：余弦相似度 → `inter = cos·√(n1·n2)`，
  集合大小用文件头存储的 n_kmer。

### 3.3 基准与验证记录

* [[../benchmarks/bench-simd-hv-jaccard.md]]（2026-08-06/07/08）：AVX2 bit
  主路径 1.11 ms（n=10k、D=4096；2026-08-08 升级为 64 位用满后实测），
  相对旧 wide bit ~6.0×、相对旧 wide i8 ~3.9×；i8 保语义 2.11 ms
  （~2.1×）；限速分解、累加器宽度实验、幅度模拟、RNG 候选对比结论均
  记录于此。
* [[../benchmarks/dist-cohort-validation.md]]（10 株 E. coli × 45 对）：
  稠密 i8 饱和（±1.3e6，Spearman ≈ 0）→ **稀疏 v2 修复**（与 `dist pgi`
  mash 排序 Spearman 0.969、45 对 0.12 s、共享计数平均误差 2.39%）；
  `dist seq`（k=8 syncmer）仍是与身份率最贴近的草图层（ρ=0.616）。

### 3.4 已知问题（待修）

FASTA 路径 `hash_hv_i8` + `calc_distances` 存在**量纲不匹配**（§2.1 的
直流偏置所致）：`hv_cardinality = ‖H‖²/D` ≈ N·(E[b²]+(N−1)/4) 量级
（N=3000 时 ≈ 6140·N，不是 N），点积被 ~N²/4 项主导，Jaccard 随 N 增大
趋向只依赖集合大小的常数（等大时 ~0.5）。
数值模拟（N=3000、shared=500）：reported 0.154 vs 真值 0.091（已固化为
`test_hash_hv_i8_jaccard_dc_bias`）。稀疏 `.hv` 路径不受影响；该 FASTA
路径应改用 bit 编码或稀疏投影，先用两株 E. coli 对照 `dist seq` /
`dist pgi` 实测确认（§5.3）。

## 4. 外部参考（HyperGen + hdlib）

> 外部参考分析已迁至 [[../references/hv.md]]（该文件将随后续文献持续
> 扩充）。参考用途：了解 HDC 草图的主流做法与参数选择。实现决策以
> §1/§2 为准。

## 5. 后续方向与待办

### 5.1 pbit 决策 B 的参考方案

[[pbit.md]] 决策 B（HV sketch 内嵌）触发条件是"无源 FASTA、仅归档、需距离
粗筛"——正是 HyperGen 的完整场景。触发时可直接参考
[[../references/hv.md]] §1–3 的方案：FracMinHash 采样 + 稠密 ±1 编码
（i32）+ 无损量化 bitpack + 预存范数 + 点积 ANI，注意其中 i16 值域限制
与 batch 全对输出。

### 5.2 可借鉴工程点

* 无损量化 + bitpack：pgr `.hv` 现在存原始 i32，可无损压到 ~6–16 bit/维；
* 预计算 L2 范数 / 存储 n_kmer（pgr 稀疏版已做后者）；
* 对称模式只算上三角、输出按 ANI 排序 + 阈值过滤；
* 批量搜索的 GEMM 化（HyperGen 论文 future 方向，pgr `linalg` 已是 wide SIMD）；
* 结构/权重编码（hdlib 的 bind/permute、weighted bundle）：未立项，
  触发条件见 [[../references/hv.md]] §4。

### 5.3 待验证 / 待办

1. **FASTA `dist hv` 路径量纲问题**（§3.4）：先实测确认，再决定改
   `hash_hv_bit` / 稀疏投影；
2. HV 编码 SIMD 深挖（todo.md §4）：**已解决（2026-08-07，见 §2.2/§2.4）**；
3. 若决策 B 立项，按 §5.1 做设计。

## 参考

* [[pbit.md]]（决策 B 与 `.hv` 消费链）
* [[../benchmarks/bench-simd-hv-jaccard.md]]（HV 编码基准与限速分解）
* [[../benchmarks/dist-cohort-validation.md]]（饱和问题与稀疏 v2 验证）
* [[../todo.md]]（§4 HV SIMD 疑虑）
* [[../references/hv.md]]（外部参考：HyperGen 论文与代码、hdlib）
* 采样方法文献（§2.6）：Edgar 2021, *PeerJ* 9:e10805（syncmer 定义）；
  Belbasi et al. 2022, *Bioinformatics* 38(Suppl 1):i169–i176,
  doi:10.1093/bioinformatics/btac244（minimizer Jaccard 有偏且不一致）；
  Irber et al. 2022, bioRxiv 2022.01.11.475838（FracMinHash）；
  Hera et al. 2023, *Genome Res* 33(7):1061–1068,
  doi:10.1101/gr.277651.123（FracMinHash 校正 + 置信区间 / ANI）；
  Liu & Koslicki 2023, bioRxiv 2023.11.09.566463（open syncmer ≡
  FracMinHash，且重叠 k-mer 采样有依赖）；Shibuya et al. 2022, WABI,
  LIPIcs 242:14, doi:10.4230/LIPIcs.WABI.2022.14（closed syncmer
  无偏 Jaccard）；Kille et al. 2023, *Bioinformatics* 39(9):btad512
  （minmer：minimizer 偏差修正）。
