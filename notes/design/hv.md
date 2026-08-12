# HV（Hypervector）设计笔记

> 2026-08-07 建立，持续更新。本文顺序：算法概览（§1）→ 影响算法设计的
> 因素（§2，核心决策记录）→ pgr 落地现状（§3）→ HyperGen 论文与代码
> 参考（§4，外部参考：HyperGen + hdlib，背景材料）→ 后续方向（§5）。
> 来源：外部参考详细分析见 [[../references/hv.md]]（HyperGen 论文与代码、
> hdlib；该文件将随后续文献持续扩充）。
> * 基准与验证：[[../benchmarks/bench-simd-hv-jaccard.md]]、
>   [[../design/hv.md]]。

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
  （cohort 实测与身份率 Spearman ≈ 0，见 [[../design/hv.md]] §2）。

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

**决策**：bit（±1）为主实现；i8 仅作为"保语义"变体保留（2026-08-08 起
FASTA 路径已改用 bit，量纲问题修复见 §3.4；i8 不再有生产消费者）。

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
加上 broadcast 摊销）。两条"减少每 bit 指令数"的替代路线均被实测否决：
i16 lane（见 §2.3，慢 2.5–3.3×）与 **pshufb 4-bit 查表**（2026-08-08
实测 5.21 ms vs 主实现 1.11 ms，慢 ~4.7×——字节级 LUT 展开需要标量
nibble 提取 + `vpmovzxbd` 字节→i32 扩展，每 bit op 数反而更多）。
srlv 家族已接近该方向的最优。

**决策**：进一步提速应从"减少每 bit 指令数"入手，而不是换 RNG 或换
lane 宽度——该方向已被 i16/pshufb 实测否定，剩余杠杆在广播/依赖链
（const 实验显示其占 ~98% 成本）。AVX-512 参考实现（旧式每步 sub，已
同步 64 位用满）Zen 4 实测 1.44 ms，反而比当前 AVX2 主路径（1.11 ms）
慢，作者决策不采用 AVX-512（见 [[../benchmarks/bench-simd-hv-jaccard.md]]
§2）。

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
  k-mer（bioRxiv 2022 提出，Irber）。每个 k-mer 的
  保留概率相同且独立 → 保留集合近似均匀随机子集。Genome Res 2023 证明
  其 Jaccard / containment 估计的偏差很小且
  **可精确校正**（给出校正公式），并据此为突变率（即 ANI 语义）推导
  **点估计 + 置信区间**（论文演示参数 k=21、scale 0.1）。
* **minimizer**（pgr FASTA 路径的另一选项）：窗口内保留最小 hash。
  保留概率与 k-mer 在窗口中的哈希排名 / 局部序列结构相关，不是均匀子
  采样。Bioinformatics 2022 证明
  minimizer 的 Jaccard 估计**有偏且不一致**——偏差不随序列长度增长
  消失，估计收敛到与真实 J 不同的值。MashMap 旧版（minimizer 版）的
  ANI 估计即受此影响，minmer 方案（Bioinformatics 2023）专门为修正它
  而提出。
* **closed syncmer**（pgr 主力）：窗口内最小 s-mer 落在首/末位
  （PeerJ 2021, Edgar；syng 移植，见
  [[../references/syng.md]]）。文献支持分两层：
  * **open syncmer 与 FracMinHash 等价**（bioRxiv 2023, Koslicki）：
    在 k-mer 集合相似度（Jaccard /
    containment）意义上等价，且距离分布与保守性更好。关键假设：
    syncmer 的采样对**重叠 k-mer 有依赖**（共享子串影响被选概率），
    FracMinHash 则完全独立，等价性正由此论证；
  * **closed syncmer 的无偏性有同行评审直接证据**：Shibuya 等 2022
    （WABI 2022）用与 pgr
    同规则的 closed syncmer 做集合比对，称其为 minimizer 的
    "上下文无关"替代、给出无偏 Jaccard 估计；作者博士论文实验里
    minimizer 明显有偏，syncmer 与随机采样在无偏性上重合。注意：
    closed syncmer 没有 FracMinHash 那种可写出的偏差量级、校正与
    置信区间公式，"无偏"停留在论证 + 实验层面。

**项目实测佐证**（[[../design/hv.md]]）：

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

> **实证（2026-08-08，`dist seq --sampler frachash` 落地后）**：用全
> k=40 集合（不经采样）对照，e2348 × cft073 的真值 Jaccard = 0.451，
> FracMinHash 估计 0.417（**无偏**），而 `dist pgi`（syncmer 8/5）仅
> 0.095——syncmer/minimizer 的位置采样在近缘（含共享重复元件）基因组间
> 有系统性偏差，FracMinHash 的独立随机采样不受影响。**数值 ANI 必须用
> FracMinHash**（详见 [[../design/hv.md]]）。

### 2.7 稀疏投影（`.hv` v2 路径）

> **定位说明（2026-08-08 澄清）**：稀疏编码**不是有意设计的产品方案**，
> 而是历史对话中在"尽全力提高速度"的压力下形成的实现（为绕开稠密 i8
> 的饱和问题与速度瓶颈）。它当前是 `.hv` v2 的实际路径（`pgi to-hv` →
> `dist hv`），但**应视为待重新审视的候选**，而非确定的实施方案——作为
> 实施方案依据时需与稠密 bit 路线重新权衡（见 §5.4）。

第三种编码路径，与稠密 bit/i8 并列：`hash_hv_sparse` 用 splitmix64 派生
随机维度，**每个 seed 只更新 s 个随机维度**（±1），`.hv` v2 采用
（`pgi to-hv` → `dist hv`）。

**理论依据（稀疏随机投影）**：每个元素经哈希映射到 s 个随机维度、±1
计数，属于稀疏随机投影 / 特征哈希家族——Achlioptas (2003, PODS/JCSS)
证明随机 ±1 矩阵保持内积结构（JL 型保距）；Li 等 (2006, KDD)
证明**稀疏化**随机投影仍保持该性质；feature hashing (ICML 2009) 是同一
模式的工业实践。无偏性可严格推导：设两集合共享 `shared` 个元素，每个
元素选某维概率 s/D、符号 ±1 均匀，则 `E[dot] = shared·s`、
`E[‖H‖²] = n·s`，余弦期望 `cos ≈ shared/√(n₁·n₂)`；且固定 D 时
误差需区分两个方差视角（完整推导见
[[../references/hv.md]] §6.5，仿照 DotHash Theorem 2）：
* **投影随机性视角**（固定集合、独立投影）：归一化后 `E = shared` 无偏，
  `Var = shared/s + (|A||B| − 2·shared)/D`——**s 影响方差**（shared/s
  项），但典型参数下占比 ≤18%（s=1）/ <7%（s≥3），且随 s 增大快速
  消失；
* **集合随机性视角**（不同集合对、单次投影——pgr 实际评估方式）：
  相对方差 ≈ n²/(shared²·D)，**与 s 无关**（s 消去）；pgr 投影是确定性
  的（seed 来自 k-mer hash），跨集合评估时集合随机性主导，s 的影响被
  淹没——50 组扫描 MAE 平坦即此。
* 结论：**s 对期望完全无影响、对方差有小的有限影响**（<20%）；"s 只
  决定速度"在集合随机性视角下成立，投影随机性视角下 s 还有一个小
  方差项。注：Li 2006 的"稀疏度放大方差"针对固定数据维度的 JL 降维
  语境，与上述公式不同但方向一致（s 小方差大）。

> **产品先例（诚实标注）**：该思想（稀疏随机投影 / feature hashing）
> 在工业界有先例（Google feature hashing、SimHash、count sketch），但
> "基因组 k-mer 集合 → 每元素 s 个 ±1 桶 → 余弦 ≈ 集合重叠"这个具体
> 组合**无直接产品对标**：HyperGen 用稠密（[[../references/hv.md]]）、
> Mash/sourmash 用 minhash、Dashing 2 用 SetSketch。pgr 稀疏路径的
> 合理性目前主要靠 cohort 实测（见下）+ 无偏性推导，无成熟产品背书。

**性能**（2026-08-08，n=10k、D=4096，AVX2 无关的标量路径）：

| s | medium 耗时 | vs 稠密 bit（1.11 ms） |
|---:|---:|---:|
| 1 | 0.022 ms | ~50× 快 |
| **3（默认）** | **0.055 ms** | **~20× 快**（vs i8 ~38×） |
| 8 | 0.155 ms | ~7× |
| 16 | 0.365 ms | ~3× |
| 64 | 1.622 ms | ~0.7× |

耗时近似随 s 线性（每 seed s 次 splitmix + 内存更新）。

**距离语义**：稀疏 HV 的每维是碰撞计数，`‖H‖²/D ≈ n²s²/D²`，稠密
cardinality 公式不成立；`.hv` v2 用**余弦 + 文件头存储的 n_kmer**
（`inter = cos·√(n1·n2)`），与稠密 Jaccard 公式不同。

**s 参数与估计质量（2026-08-08 系统扫描）**：50 组独立随机集合对
（N=3000、shared=500、D=4096），s=1..4096 全范围扫描，**MAE 完全平坦**
（0.010–0.013，无单调趋势，s=4096 与 s=1 误差相同）——**s 不是精度
杠杆**：误差由 D 决定（~1/√D，随机投影的固有噪声），s 只决定速度。
固化为 `test_hash_hv_sparse_s_error_scan` / `_jaccard_s_scan`。
cohort 验证（s=3）与 `dist pgi` mash 排序 Spearman 0.969、共享计数
平均误差 2.39%（见 [[../design/hv.md]] §2）。

**D 与编码成本的解耦（2026-08-08 关键实验）**：稀疏编码成本
**O(n·s)，与 D 无关**（s=1 时 D=4096/16384/65536 编码 0.022/0.023/
0.024 ms，几乎不变）；精度 `MAE ∝ 1/√D`（D×4 → MAE÷2.5，D×16 →
MAE÷5.2，实测）。对照稠密 bit 编码 O(n·D)（D 翻 16 倍编码也 ~16×）。
固化为 `test_hash_hv_sparse_d_error_scan`。

**决策（2026-08-08 重估）**：

* **s 取 1 而非 3（已落地 2026-08-08）**：s 不影响精度，s=3 是历史默认
  （commit d967d16, 2026-08-02）无依据；`pgr pgi to-hv` 默认已改 1
  （0.022 ms，s=3 的 2.5× 快）。cohort 复测（5 株 × 10 对，2026-08-08）：
  s=1 与精确 `dist pgi` mash 排序 Spearman 0.988、最大差异 0.0025
  （见 [[../design/hv.md]]）。**完整 45 对复测（2026-08-12，10 株本地
  数据）：mash 排序 Spearman 0.9814、Pearson 0.9969、max |Δ| 0.0035**
  ——全 cohort 上排序保真度略低于 10 对子集但仍很高，s=1 默认成立。
* **大 s 是错误用法**：s=64 时 1.62 ms，比稠密 bit（1.11 ms）还慢——
  稀疏实现是随机内存访问（s 线性、缓存不友好），稠密是连续 SIMD；
  s 大既不提精度又更慢，纯亏损。稀疏的有效区间是 s 小（1–3）。
* **精度靠 D 而非 s**：MAE ≈ 常数/√D，D 可按精度需求标定；稀疏用大 D
  免费换精度（编码不变），稠密不行。
* **是否用稀疏**（§5.4 决策）：编码是瓶颈的大规模场景（4 万 cohort）→
  稀疏 s=1 + 大 D 有实质优势；小规模 / 需 Jaccard 数值语义 → 稠密 bit。
* 稀疏路径不受 §2.3 的累加器宽度问题与 §3.4 的直流偏置问题影响。

## 3. 实现现状（pgr 落地）

> **当前推荐配置速查**：编码用 bit（±1）+ i32 累加 + AVX2 跳步 RNG
> （块主序）、D=4096；FASTA 侧采样用 closed syncmer（`dist seq` 默认
> k=8/w=5）；需要可标定 ANI 时优先 FracMinHash（§2.6）。⚠️ FASTA
> `dist hv` 路径（`load_hv_from_fasta` / `load_hv_from_fasta_syncmer`）
> 已改用 bit 编码（2026-08-08，量纲问题修复，§3.4）；`.hv` 索引路径走
> 稀疏 s=1（2026-08-08 默认 3→1，s 不影响精度），其选择为历史产物、
> 待 §5.4 重新审视。

### 3.1 `src/libs/hv.rs` 当前实现与分派

| 函数 | 编码 | 说明 |
|---|---|---|
| `hash_hv_bit` | 稠密位编码 ±1，i32 | AVX2 主路径：跳步 RNG + 块主序 + 延迟 −N；每 64 维 1 次 RNG（用满 64 位，低 32 位→前 32 维、高 32 位→后 32 维） |
| `hash_hv_i8` | 稠密 i8 累加，i32 | AVX2 路径：每 8 维 1 次 RNG；保留原语义（含直流偏置） |
| `hash_hv_sparse` | 稀疏 ±1，i32 | splitmix64 派生，每 k-mer 只更新 `s` 个随机维度；`.hv` v2 采用（历史产物，待 §5.4 重新审视；2026-08-08 默认 s 3→1） |
| `hv_norm_l2_sq` / `hv_cardinality` / `hv_dot` | — | wide SIMD 范数；cardinality=‖H‖²/D；dot 按 √D 归一 |
| `calc_distances` | — | jaccard / containment / mash 多口径输出 |
| `load_hv_from_fasta(_syncmer)` | bit | FASTA → minimizer / closed syncmer 集合 → `hash_hv_bit`（2026-08-08 由 i8 改） |

分派链（x86-64）：`avx2` 运行时检测 → AVX2 intrinsics；否则 wide 可移植
路径（aarch64 自动用 NEON，其余标量）。所有路径输出**逐位一致**（含
非 32 倍数维度的尾部处理），由 `test_hash_hv_bit_serial_vs_simd` /
`test_hash_hv_i8_serial_vs_simd` 及 avx2 变体保证。
AVX-512 实现只保留在 `benches/hv_benchmark.rs` 作参考对照，不参与分派。
`cargo check --target aarch64-apple-darwin --lib` 交叉编译通过。

### 3.2 消费链与 `.hv` 格式

* `pgr dist hv`：FASTA 路径（minimizer/syncmer → `hash_hv_bit` →
  `calc_distances`）；`.hv` 路径（`pgi to-hv` 产物直接比较，稀疏余弦）。
* `pgr pgi to-hv`：把 `.pgi` 的 unique k-mer keys 投影成**稀疏** HV；
  `.hv` v2 格式：`PGV1` magic + version 2 + `k/dim/sparse/n_kmer/name` +
  i32 数组。稀疏投影 + 存储真实 n_kmer 是 v2 的关键修复（注：稀疏是
  历史性能优化压力下的产物，非有意设计，见 §2.7 定位说明 / §5.4）。
* `pgr dist hv a.hv b.hv`：余弦相似度 → `inter = cos·√(n1·n2)`，
  集合大小用文件头存储的 n_kmer。

### 3.3 基准与验证记录

* [[../benchmarks/bench-simd-hv-jaccard.md]]（2026-08-06/07/08）：AVX2 bit
  主路径 1.11 ms（n=10k、D=4096；2026-08-08 升级为 64 位用满后实测），
  相对旧 wide bit ~6.0×、相对旧 wide i8 ~3.9×；i8 保语义 2.11 ms
  （~2.1×）；限速分解、累加器宽度实验、幅度模拟、RNG 候选对比结论均
  记录于此。
* [[../design/hv.md]]（10 株 E. coli × 45 对）：
  稠密 i8 饱和（±1.3e6，Spearman ≈ 0）→ **稀疏 v2 修复**（与 `dist pgi`
  mash 排序 Spearman 0.969、45 对 0.12 s、共享计数平均误差 2.39%）；
  `dist seq`（k=8 syncmer）仍是与身份率最贴近的草图层（ρ=0.616）。

### 3.4 已知问题（FASTA 量纲问题已修复 2026-08-08）

FASTA 路径原走 `hash_hv_i8` + `calc_distances`，存在**量纲不匹配**
（§2.1 的直流偏置所致）：`hv_cardinality = ‖H‖²/D` ≈ N·(E[b²]+(N−1)/4)
量级（N=3000 时 ≈ 6140·N，不是 N），点积被 ~N²/4 项主导，Jaccard 随
N 增大趋向只依赖集合大小的常数（等大时 ~0.5）。数值模拟（N=3000、
shared=500）：reported 0.154 vs 真值 0.091（`test_hash_hv_i8_jaccard_dc_bias`
记录该缺陷）。

**已修复（2026-08-08）**：`load_hv_from_fasta` / `load_hv_from_fasta_syncmer`
改用 `hash_hv_bit`（±1 平衡），量纲问题消失——bit 的 `‖H‖²/D ≈ N`，
模拟实测 Jaccard 0.102 vs 真值 0.091（`test_hash_hv_bit_jaccard_accurate`）；
两株 E. coli（MG1655 × Sakai）`dist hv --sampler syncmer` 实测输出合理
（jaccard 0.927、containment 0.989、不再饱和）。顺带提速（bit 1.11 ms
vs i8 2.11 ms）。稀疏 `.hv` 路径不受影响。

## 4. 外部参考（HyperGen / hdlib / 测距聚类文献）

> 外部参考分析已迁至 [[../references/hv.md]]（该文件将随后续文献持续
> 扩充）：§1–4 为 HyperGen 与 hdlib（HDC 主流做法），§5 为测距/聚类
> 文献。参考用途：了解 HDC 草图的主流做法与参数选择。实现决策以
> §1/§2 与 §6 审计为准。

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

1. **FASTA `dist hv` 路径量纲问题**（§3.4）：**已修复（2026-08-08）**——
   `load_hv_from_fasta` / `_syncmer` 改用 `hash_hv_bit`，模拟 + 两株
   E. coli 实测确认（见 §3.4）；
2. HV 编码 SIMD 深挖（todo.md §4）：**已解决（2026-08-07，见 §2.2/§2.4）**；
3. 若决策 B 立项，按 §5.1 做设计。

### 5.4 重新审视 `.hv` 路径的稀疏选择（2026-08-08 立项）

稀疏投影是历史性能优化压力下的产物（§2.7 定位说明），当前是 `.hv` v2
实际路径但无产品先例、验证以 cohort 实测为主。作为实施方案前应重新
权衡：稠密 bit（§2.1–2.5 全套调优 + 64 位用满，1.11 ms）vs 稀疏 s=1
（0.022 ms，~50× 快，且编码成本与 D 无关——可用大 D 免费换精度，见
§2.7 关键实验；代价是距离语义不同、无先例）。决策点：

* **编码是否瓶颈**：大规模 cohort（4 万基因组，编码次数多）→ 稀疏 s=1
  的 ~50× 编码优势 + 大 D 免费精度是实质性的；小规模 → 稠密更简单。
* **距离语义**：稀疏的余弦 + n_kmer 是否可接受，或改回稠密 Jaccard
  （数值 ANI / 可标定场景更倾向稠密或 FracMinHash）。
* **若保留稀疏**：s=1（s 不影响精度，见 §2.7）；D 按 MAE ≈ 常数/√D
  标定；不用中间 s（大 s 比稠密慢且无精度收益）。

### 5.5 文献驱动的未来方向（2026-08-08）

基于 [[../references/hv.md]] §5 的文献审计（评估详见 §6），按优先级：

1. **DotHash 误差界 → D/s 选择理论**（Nunes 2023, §5.1）：把 DotHash 的
   误差概率界引入 pgr，推导"给定误差容忍度所需的 D"（替换 D=4096 的
   经验默认），并评估"每元素 s 个随机 ±1 桶"与 s=1 / 稠密的误差对比。
2. **FracMinHash 落地**（Irber 2022 / Hera 2023）：**已实现（2026-08-08）**——
   `dist seq --sampler frachash`（FracMinHash 采样器：canonical k-mer 保留
   hash < u64::MAX/scale，k 默认 21/7、`--scale` 默认 1000）；MG1655×Sakai
   ANI 估计 97.7% vs 真值 97.3%，scale=1000/100 一致；排序与全 k=40
   集合真值 Spearman 1.0（见 [[../design/hv.md]]）；
   **`--ci` 输出 ANI 95% 置信区间**（正态近似；Hera 校正公式留作后续）。
3. **minmer 替代 minimizer**（Kille 2023）：无偏 Jaccard + MashMap 10×
   快；`seq_mins` 的 minimizer 采样（dist seq / fa）可升级，消除 §2.6
   的 minimizer 偏差。
4. **Yu 2022 conservation 理论**：采样器选择的定量框架（syncmer 闭式解），
   补 §2.6 的"排序/粗筛够用"论断。
5. **ProbMinHash 加权 Jaccard**（Ertl）：若支持 k-mer 多重度（拷贝数），
   对应 §4.4 的 weighted bundle 方向。
6. **大规模聚类集成**（§5.2）：RabbitTClust（4 万 cohort 对标）、
   GSearch / HNSW（ANN 搜索）。
7. **采样层新方法**（§5.3）：spaced seeds（ntHash2）、strobemers、
   mod-minimizer（长 k-mer），待测距应用需求驱动。

## 6. 决策理论依据审计（2026-08-08）

对照 [[../references/hv.md]] §5 文献，逐条审计 §2 的 7 个算法决策。
评估等级：强（定理/证明直接支持）、中（思想支持但参数无依据）、弱
（工程决策，文献无指导）。汇总：

| # | 决策 | 评估 | 关键依据 / 缺口 |
|---|---|---|---|
| 1 | 编码：bit ±1 为主 | 中偏强 | Kanerva 2009 准正交性、SimHash/DotHash ±1 超向量；i8 无文献先例 |
| 2 | RapidRng | 弱（工程合理） | HDC 文献只要求"随机"，不约束具体 RNG |
| 3 | i32 累加器 | 强 | 值域数学证明 + HyperGen i16 溢出教训 |
| 4 | srlv 展开 | 无理论（工程） | i16/pshufb 实测否决，实验闭环 |
| 5 | 幅度无关 | 强 | Jaccard 比值度量 + DotHash 点积框架 |
| 6 | syncmer + FracMinHash | 最强 | 8 篇文献；缺口：Yu 2022 未引用、minmer 仅一句 |
| 7 | 稀疏（默认 s=1） | **强（理论已补齐）** | **仿 DotHash 推导完整：无偏 + Var = shared/s + (|A||B|−2·shared)/D**（references/hv.md §6.5）；s/D 选择可标定，默认已改 1 |

### 6.1 逐条详情

1. **编码路径 bit（±1）**：Kanerva 2009（[[../references/hv.md]] §5.4）
   的"高维随机 ±1 向量准正交"是 bit 的理论源头；SimHash 与 DotHash
   （§5.1，Theorem 2：点积无偏估计交集）同属该框架。i8 字节编码在 HDC
   文献中无先例，直流偏置是 pgr 独有历史实现（§3.4）——文献支持
   "围绕 0 平衡的 bipolar 表示"，即 bit 路径。
2. **RapidRng**：HDC 文献只要求"均匀 + 独立"的随机向量；8 个 RNG 的
   实测对比（[[../benchmarks/bench-simd-hv-jaccard.md]]）是工程证据。
3. **i32 累加器**：值域确定性分析（bit 路径 [−N,N] 不溢出）+ HyperGen
   i16 溢出教训（[[../references/hv.md]] §2）。
4. **srlv 展开**：纯 SIMD 工程；i16 / pshufb 替代实测否决，无理论缺口。
5. **幅度与区分度无关**：Jaccard 比值度量的数学性质；DotHash 的点积/
   范数框架支持"点积估计交集、范数估计基数"。
6. **采样方法**：最强文献支撑——Bioinformatics 2022（minimizer 有偏）、
   PeerJ 2021（syncmer, Edgar）、bioRxiv 2022（FracMinHash, Irber）、
   bioRxiv 2023（syncmer≡FracMinHash, Koslicki）、Bioinformatics 2023
   （minmer 无偏）、WABI 2022（closed syncmer 无偏）、Genome Res 2023
   （FracMinHash 校正）。缺口：**Yu 2022（Bioinformatics, local
   k-mer selection conservation 理论，minimap2 实证
   8.2%）未引用**；minmer 条目可强化（无偏 + 10× 快）。
7. **稀疏投影（s 默认已改 1）**：DotHash Theorem 2 直接支持"随机超向量叠加点积无偏
   估计交集"。**2026-08-08 已仿照 Theorem 2 为 pgr 的 s 桶 ±1 构造补全
   理论**（references/hv.md §6.5）：无偏性 E = shared（√s 归一化后）、
   方差 Var = shared/s + (|A||B|−2·shared)/D、Chebyshev 误差界可直接套用
   （`Pr(|est−shared| ≥ ε) ≤ Var/ε²`）——**稀疏投影现在有 DotHash 级
   推导**，非仅有实测。**文献家族确认（references/hv.md §6.6）**：Li
   2006 稀疏分布（±√s、非零概率 1/s）与 pgr 期望等价、Weinberger 2009
   特征哈希指数尾界、Count-Min Sketch 桶结构与 pgr 同构、Achlioptas
   2001 稀疏保距奠基——s 桶构造不是孤立方案。s 的影响是小方差项
   （典型参数 <20%）；D 是精度主导（1/√D）。原缺口（无偏/误差界）
   已闭合；s 的历史默认已按 §2.7 决策改为 1（2026-08-08），按需标定 D。

## 参考

* [[pbit.md]]（决策 B 与 `.hv` 消费链）
* [[../benchmarks/bench-simd-hv-jaccard.md]]（HV 编码基准与限速分解）
* [[../design/hv.md]]（饱和问题与稀疏 v2 验证）
* [[../todo.md]]（§4 HV SIMD 疑虑）
* [[../references/hv.md]]（外部参考：HyperGen / hdlib（§1–4）+ 测距聚类
  文献（§5），§2.6 采样方法与 §6 审计的文献链）
* 采样方法文献（§2.6）：PeerJ 2021（syncmer 定义, Edgar）；
  Bioinformatics 2022（minimizer Jaccard 有偏且不一致）；
  bioRxiv 2022（FracMinHash, Irber）；
  Genome Res 2023（FracMinHash 校正 + 置信区间 / ANI）；
  bioRxiv 2023（open syncmer ≡ FracMinHash, Koslicki；且重叠 k-mer 采样
  有依赖）；WABI 2022（closed syncmer 无偏 Jaccard）；Bioinformatics
  2023（minmer：minimizer 偏差修正）。

---

## 证据附录：HV 距离 vs ANI 金标准标定（真实 Enterobacterales 基因组）

> 目的：给 `design/genome-nn-query.md` §7 的 P1 补第一块硬证据——HV
> 距离与真实 ANI 的相关性、分辨率区间与召回。日期：2026-08-08。
> 结论先行：HV 距离只在 ANI 90–98% 区间中等可靠（Spearman 0.5–0.6），
> **≥99% 与 <85% 失效**；同条件下 Mash 与 ANI 几乎完全相关（ρ≈0.97–0.99）。
> 物种内（≥98% ANI）选参考 / 聚类不应以 HV 距离为主。

## 数据与 cohort（严格挂靠 NWR 指导文件）

- 数据源：`~/data/Escherichia/`（Enterobacterales + Pasteurellales，
  150k 组装，132,572 个通过 QC，见 `~/Scripts/genomes/groups/Escherichia.md`）。
- cohort：135 个基因组，全部来自各物种 **NR.lst**（非冗余代表）且都在
  `summary/pass.lst`（QC 通过）内；物种标签取自 `summary/genome.taxon.tsv`。
  亲缘分层：E. coli NR 40（近缘，ANI≥98%）+ 其他 Escherichia 种 60
  （E. albertii/fergusonii/marmotae/ruysiae/whittamii/sp、Pseudescherichia，
  中等，ANI 88–97%）+ Yersinia 36（远缘，ANI<88%）。全两两 = 9,045 对。
- 样本清单与映射：`/tmp/hv_calib/cohort.meta.tsv`（name/species/path，
  临时目录，可重建）。

## 方法

- **HV 距离**：`pgr dist hv --list-files --parallel 8`，默认 DNA minimizer
  k=21/w=5；D=4096 与 D=16384 各跑一遍（输出第 7 列 = Mash 式距离
  d = −(1/k)ln(2J/(1+J))，1−d ≈ ANI 估计值）。
- **Mash 距离**：直接用 NWR 已算好的每基因组 `.msh` sketch
  （`MinHash/<species>/msh/<name>.msh`），`mash triangle -E -p 8`（与
  NWR `dist.sh` 同款参数），9045 对全覆盖。
- **ANI 金标准**：skani 0.1.0 `dist --ql/--rl -t 8 --min-af 0`（取满覆盖；
  极远缘无比对命中者视为未知，5,937/9,045 对有 ANI）。skani 是全基因组
  ANI（近似 BLAST-ANI），本文以它为真值；GSearch 用 BLAST-ANI/FastANI，
  二者同级别。
- **指标**：Spearman（排序）/Pearson/RMSE（HV 估计 ANI vs skani ANI，
  RMSE 按 0–1 标度）；recall@10 = 以 ANI 为真值取 top-10，与 HV/Mash
  距离 top-10 的交集比例（自比对剔除，按查询新颖度分层）。分析脚本：
  `/tmp/hv_calib/analyze_ani.py`、`analyze_dim.py`（临时）。

## 结果

### HV(1−d) vs skani ANI（D=4096）

| 分层 | n | Spearman | Pearson | RMSE(0–1) |
|---|---|---|---|---|
| 全部 | 5,937 | 0.882 | 0.481 | 0.191 |
| 同种内 | 1,115 | 0.610 | 0.124 | 0.168 |
| 种间 | 4,822 | 0.861 | 0.504 | 0.195 |
| ANI ≥99% | 122 | **0.383** | 0.238 | 0.045 |
| ANI 95–99% | 1,339 | 0.608 | 0.115 | 0.154 |
| ANI 90–95% | 2,427 | 0.496 | 0.184 | 0.072 |
| ANI 85–90% | 721 | 0.462 | 0.378 | 0.137 |
| ANI <85% | 1,328 | **0.054** | 0.068 | 0.344 |

### Mash vs skani ANI（参考，同数据）

| 分层 | n | Spearman | Pearson |
|---|---|---|---|
| 全部 | 5,937 | −0.990 | −0.982 |
| 同种内 | 1,115 | **−0.974** | −0.987 |
| 种间 | 4,822 | −0.983 | −0.977 |

### HV 距离分位数（分辨率直观检查，D=4096）

| ANI 区间 | n | hv_dist q05 / q50 / q95 |
|---|---|---|
| ≥99% | 122 | 0.016 / 0.043 / 0.065 |
| 95–99% | 1,339 | 0.034 / 0.058 / 0.099 |
| 90–95% | 2,427 | 0.087 / 0.109 / 0.139 |
| 85–90% | 721 | 0.108 / 0.126 / 0.206 |
| <85% | 1,328 | 0.152 / 0.202 / 1.000 |

相邻 ANI 区间的 hv_dist 分布大量重叠（95–99% 与 90–95% 的 q05–q95
区间几乎连续）——排序分辨率差的直接体现。

### recall@10（真值 = skani ANI top-10，135 个查询）

| 方法 | 总体 | 新颖度 ≥98% | 新颖度 90–98% |
|---|---|---|---|
| HV（D=4096） | 0.622 | 0.612 (n=124) | 0.727 (n=11) |
| Mash | 0.762 | 0.762 (n=124) | 0.764 (n=11) |

### D=4096 vs D=16384（分辨率是否随维度改善）

| ANI 区间 | Spearman D=4096 | Spearman D=16384 |
|---|---|---|
| ≥99% | 0.383 | 0.387（**无改善**） |
| 95–99% | 0.608 | 0.608（无改善） |
| 90–95% | 0.496 | 0.577 |
| 85–90% | 0.462 | 0.558 |
| <85% | 0.054 | 0.322 |
| recall@10 总体 | 0.622 | 0.629 |

## 结论

1. **HV 距离的可靠区间是 ANI 90–98%（及部分 85–90%），且仅中等
   （Spearman 0.5–0.6）**；≥99% 的近缘株排序几乎失效（ρ≈0.38，且
   D=16384 无改善——不是维度饱和度，是方法固有噪声）；<85% 远缘完全
   失效（ρ≈0.05–0.32）。
2. **Mash 是 ANI 的可靠代理**（同种内 ρ=−0.97），recall@10 比 HV 高
   14 pp——种内近缘排序任务 Mash 明显更优。
3. **对设计的直接影响**：物种内（≥98% ANI）聚类 / 选参考应以
   `dist mash` / `dist frac` 为主；HV 适合做**嵌入 / 粗筛 / 查询路由**
   （85–98% 带），不适合做 ANI 级精排，更不能替代 skani/fastANI 标定。
4. 增大 D（4096→16384）只改善中远缘，代价是 4× 内存/时间；对
   近缘分辨率无帮助，参数上不必为近缘场景升级 D。
5. 待补：完整度鲁棒性、sampler/k 扫描、HNSW 检索在真实 HV 上的
   ANI-truth 召回（§7.2 ② 的图检索部分）。

## 复现

```bash
mkdir -p /tmp/hv_calib
# 1. 从 NR.lst/pass.lst 抽样 135 基因组 -> cohort.meta.tsv / cohort.fa.lst
# 2. HV:  pgr dist hv cohort.fa.lst --list-files --parallel 8 -o hv.tsv
# 3. Mash: mash triangle -E -p 8 -l cohort.msh.lst > mash.tsv
# 4. ANI:  skani dist --ql cohort.fa.lst --rl cohort.fa.lst -t 8 --min-af 0 -o ani.full.tsv
# 5. 分析: python3 analyze_ani.py / analyze_dim.py
```

## 补充：四种距离统一对标 ANI（#1，2026-08-08）

同 cohort、同 ANI 真值，追加 `pgr dist frac --merge`（k=21, scale=1000）
与 `pgr dist mini --merge`（k=21/w=5, rapid）。Spearman 绝对值：

| 方法 | 全部 | 同种内 | ≥99% | 95–99% | 90–95% | 85–90% | <85% |
|---|---|---|---|---|---|---|---|
| HV | 0.882 | 0.610 | 0.383 | 0.608 | 0.496 | 0.462 | 0.054 |
| Mash | 0.990 | 0.974 | 0.805 | 0.972 | 0.961 | 0.924 | 0.676 |
| **frac** | **0.991** | **0.973** | 0.796 | **0.970** | **0.961** | 0.919 | 0.673 |
| mini | 0.917 | 0.612 | 0.401 | 0.607 | 0.632 | 0.583 | 0.600 |

recall@10（真值 = skani ANI top-10）：

| 方法 | 总体 | ≥98% | 90–98% |
|---|---|---|---|
| HV | 0.621 | 0.612 | 0.727 |
| Mash | 0.762 | 0.762 | 0.764 |
| frac | 0.757 | 0.759 | 0.736 |
| mini | 0.629 | 0.613 | 0.809 |

**结论补充**：① `dist frac` 与 Mash 同为 ANI 的可靠代理（ρ 0.97–0.99，
recall 与 Mash 相当），"frac 用于 ANI 估计"的既有建议得到实证支持；
② **minimizer 采样（mini）与 HV 有相同的近缘分辨率缺陷**（同种内与
≥99% 区间 ρ≈0.4–0.6）——问题是采样层（minimizer）而非 HV 编码本身；
③ mini 在 <85% 区间（ρ 0.60）明显好于 HV（0.05），远缘下 minimizer
草图仍保有信息，HV 的 4096 维编码是远缘失效的主因。

---

## 证据附录：HV 最近邻：HNSW vs 精确扫描的召回率与延迟（4096 维）

> 目的：回答 `design/genome-nn-query.md` §6.4 的待定问题——4096 维 HV 下
> ANN（HNSW）的召回率是否值得换取查询加速。日期：2026-08-08。
> 硬件：AMD Ryzen 9 7945HX（单线程，release 构建，LTO）。

## 方法

- **合成数据**：每个"基因组" = 2048 个共享 k-mer 核心 + 512–4096 个私有
  随机 u64 哈希（按基因组独立种子），用真实 HV 管线 `hash_hv_bit` 编码为
  4096 维 `Vec<i32>`，再 L2 归一化为 f32。归一化后欧氏距离排序 ≡ cosine
  排序（HV 近邻排序的实际语义，见 `design/genome-nn-query.md` §6.1）。
  共享核心模拟物种内近缘基因组的结构。
- **精确基线**：`linalg::dot_product`（f32x8 SIMD）+ 全量排序取 top-10。
  这是保守估计——生产精确扫描可用 top-k 堆（O(N log k)），会更快。
- **HNSW**：rust-cv `hnsw` 0.11.0（纯 Rust，dev-dependency），
  M=12 / M0=24，PRNG=Pcg64，构建 `ef_construction`=64（另跑 200 对照）；
  查询 ef ∈ {10, 20, 50, 100, 200, 400}。
- **指标**：recall@10 = |ANN top-10 ∩ 精确 top-10| / 10，50 个留出查询
  的平均；查询延迟为 criterion 每查询时间（50 查询/迭代）。
  复现：`cargo bench --bench hv_ann_recall`；环境变量
  `PGR_HV_ANN_SIZES`（逗号分隔规模）与 `PGR_HV_ANN_EFC`（构建 ef）。

## 结果（ef_c = 64）

| N | 精确扫描 µs/查询 | ef=10 | ef=20 | ef=50 | ef=100 | ef=200 | ef=400 |
|---|---|---|---|---|---|---|---|
| 1k | 445 | recall 0.980 / 138 µs | 0.982 / 141 | 0.982 / 148 | 0.982 / 148 | 0.982 / 149 | 0.982 / 147 |
| 10k | 9,960 | recall 0.854 / 403 µs | 0.904 / 458 | 0.924 / 523 | 0.946 / 569 | 0.956 / 663 | 0.956 / 656 |
| 30k | 30,300 | recall 0.734 / 654 µs | 0.788 / 797 | 0.852 / 983 | 0.886 / 1,121 | 0.906 / 1,326 | 0.916 / 1,513 |

构建时间（ef_c=64）：N=1k 1.6 s；N=10k 18.6 s；N=30k 63.3 s
（约 O(N log N)：外推 N=100k ≈ 4–5 min，N=1M ≈ 45 min，单线程）。

### ef_c=200 对照（区分"图质量"与"高维上限"）

| N | ef=10 | ef=20 | ef=50 | ef=100 | ef=200 | ef=400 |
|---|---|---|---|---|---|---|
| 10k | 0.852 | 0.904 | 0.926 | 0.950 | 0.956 | 0.956（构建 18.8 s） |
| 30k | 0.734 | 0.790 | 0.864 | 0.894 | 0.908 | 0.922（构建 67.6 s） |

## 结论

1. **4096 维下 HNSW 召回上限明显衰减**：N=30k 时即使 ef=400 + ef_c=200，
   recall@10 也只有 0.92；ef=10 时仅 0.73。N=10k 时 ef=20 即达 0.90，
   之后提升缓慢。召回上限主要由**查询侧 ef**决定，构建更用心（ef_c
   64→200）只带来 0.6–1.2 个百分点的改善——高维诅咒是主因。
2. **速度增益真实但代价是召回**：N=30k 时 HNSW 比精确扫描快 20–46×
   （0.65–1.5 ms vs 30 ms/查询）；N=10k 时快 15–25×；N=1k 时仅 ~3×。
   速度门槛在 ~10k 之后才明显。
3. **ef 收益递减**：ef 10→400，N=30k 召回 0.73→0.92，查询延迟
   0.65→1.5 ms（近似线性增长）。用 ef≈100–200 是召回/延迟折中区。
4. **对 pgr 的建议**（更新 `design/genome-nn-query.md` §6.4）：
   - ≤10k：精确扫描（SIMD 10 ms/查询）简单且无召回损失，直接够用；
   - 10k–30k 且接受 ~0.9 召回：HNSW 可考虑（0.5–1.3 ms/查询），
     否则继续精确；
   - >30k 或百万级：4096 维 HNSW 召回天花板 <0.92，构建分钟–小时级，
     应先降维（PCA 256–512 维）再评估 ANN，或回到 SQLite 精确扫描 +
     分桶/倒排路线。
5. 三个 dev-dependencies（`hnsw` / `space` / `rand_pcg`）只服务于本实验
   与后续 ANN 再评估，不进入主程序；基准代码保留以便复现与换数据重测。

---

## 证据附录：真实 cohort 上 HV 最近邻检索：精确 / HNSW / 物种路由（#6/#7）

> 日期：2026-08-08。135 个真实基因组（同 `hv.md`），
> HV 向量来自 `pgr pgi to-hv`（k=40 syncmer、D=4096、sparse=1），
> L2 归一化后以欧氏排序 ≡ cosine。真值：skani ANI top-10（生物学）与
> HV 精确 top-10（图检索误差）。对应 `design/genome-nn-query.md`
> §7.4 #6/#7。复现：`cargo bench --bench hv_ann_real`（env 指向
> /tmp/hv_calib 的数据）。

## 结果

| 变体 | recall_ANI@10 | recall_HV@10 | 平均查询 µs |
|---|---|---|---|
| 精确扫描 | 0.663 | 1.000 | 331 |
| 全局 HNSW ef=10 | 0.664 | 0.993 | 177 |
| 全局 HNSW ef=50 | 0.664 | 0.996 | 312 |
| 全局 HNSW ef=100 | 0.664 | 0.996 | 363 |
| 路由 R=1（物种，ef 无关） | 0.542 | 0.699 | 88 |
| 路由 R=2（物种） | 0.641 | 0.899 | 132 |

## 结论

1. **图检索层误差可忽略**：全局 HNSW 的 recall_HV@10 ≥ 0.993（ef≥10），
   recall_ANI 与精确扫描完全相同（0.664 vs 0.663）——HV→ANI 的 0.66
   差距全部来自**距离层**（HV 排序 vs ANI 排序），与
   `hv.md` 一致。P1 ② 收尾：检索层证据齐了。
2. **物种硬路由在本 cohort 上反而有害**：R=1 时 recall_HV 跌到 0.699，
   R=2 回升到 0.899——原因是很多物种 clade 成员 <10（E. whittamii 3、
   多数 Yersinia 种 1–4），真 top-10 必然跨 clade，硬路由漏掉真近邻。
   这补上了 §6.5 合成实验没覆盖的失败模式：**路由键的 clade 必须足够
   大（≥K 甚至数倍 K），否则路由变成截断**；小 cohort 或小 clade 下
   全局检索（或精确扫描）更稳。
3. **规模提示**：135 个基因组时精确扫描 331 µs/查询（SIMD 4096 维），
   HNSW 约 2× 更快（177 µs）——规模账维持"≤10k 精确够用"结论。
4. 局限：ANI-truth 仅统计有 ≥10 个已知 ANI 邻居的查询（远缘查询被
   排除）；路由实验的物种键来自标签，真实部署若用聚类键需保证
   clade 规模。

## 后续

- 用 E. coli NR 子集（≥1,000 成员的单物种）重测路由：此时物种内
  clade 足够大，路由应恢复 §6.5 合成实验的正收益——验证"clade 规模
  是路由生效前提"。

## 补充：路由正向规模案例（2,088 E. coli，2026-08-08）

方法：2,088 个 E. coli NR 全量 Mash 矩阵（`mash triangle`）→ necom
ward 聚类 → cut 成 C=8/16 个 clade；查询按 **HV（pgi-to-hv）点积路由**
到 top-R 个 clade 代表，路由 clade 内 **Mash 精确** top-10，recall vs
全量 Mash top-10（自比对剔除）。

| C | clade 大小 min/中位/max | R | HV 路由准确率 | 路由后 recall@10 vs 全量 Mash |
|---|---|---|---|---|
| 8 | 10/179/894 | 1 | 0.754 | 0.750 |
| 8 | 10/179/894 | 2 | 0.942 | 0.940 |
| 16 | 1/41/493 | 1 | 0.460 | 0.455 |
| 16 | 1/41/493 | 2 | 0.701 | 0.702 |

结论：
1. **路由后 recall ≈ 路由准确率**（4 组数据全部吻合到 ±0.004）——
   §6.5 的"期望召回 ≈ (1−m)·R₁"线性公式得到规模级定量确认：路由对
   则 recall≈1，路由错则 recall≈0。
2. **HV 可做种内路由先验**：C=8 时 R=1 准确率 0.754、R=2 达 0.942；
   用 2/8（25%）的库换取 94% 的 Mash recall——"HV 路由 + clade 内
   精确检索"的 §6.5 推荐形态在成员充足的 clade 上成立。
3. **clade 规模再次是前提**：C=16（中位 41、最小 1 成员）时准确率跌到
   0.46（R=1）——与 #7 小 clade 结论一致。
4. 局限：clade 划分用全量 Mash（含查询自身距离，轻微泄漏）；路由本身
   用 HV（留出），机制结论不受影响。

---

## 证据附录：HV 最近邻：hnsw_rs 多层 HNSW vs 单层 HubNSW（4096 维）

> 目的：验证 GSearch 官方推荐的高维对策——HubNSW 单层化
> （`hnsw_rs` 的 `modify_level_scale`，arXiv 2412.01940）在 4096 维 HV 上
> 是否改善召回；同时与 rust-cv `hnsw` 0.11 的召回/延迟做跨实现对比
> （后者见 `hv.md`）。日期：2026-08-08。
> 硬件：AMD Ryzen 9 7945HX（单线程，release 构建，LTO）。

## 方法

- **合成数据与精确基线**：与 `hv.md` 完全一致（同一
  种子、同一 `hash_hv_bit` 管线、同一 SIMD 点积 top-10 基线），保证两
  份结果可直接对比。
- **HNSW 实现**：`hnsw_rs` 0.3.4（GSearch 同源 Rust 重写，dev-dependency），
  M=24（所有层）、max_layer=16、`ef_construction`=64、`DistL2`
  （L2 归一化向量上 ≡ cosine 排序）；**scale 因子**：1.0（多层 HNSW）
  vs 0.2（HubNSW 单层，此时 P(level≥1) ≈ e^-16 ≈ 0）。
- **指标**：recall@10 = |ANN top-10 ∩ 精确 top-10| / 10，50 个留出查询
  的平均；查询延迟为 criterion 每查询时间。
  复现：`cargo bench --bench hv_ann_hubnsw`；环境变量
  `PGR_HV_ANN_SIZES` / `PGR_HV_ANN_EFC` / `PGR_HV_ANN_EFS`。

## 结果（ef_c = 64）

| N | 精确 µs/查询 | scale | ef=10 | ef=20 | ef=50 | ef=100 | ef=200 | ef=400 |
|---|---|---|---|---|---|---|---|---|
| 1k | 343 | 1.0 | 0.954 / 706 | 0.988 / 1,019 | 1.000 / 1,195 | 1.000 / 1,299 | 1.000 / 1,281 | 1.000 / 1,303 |
| 1k | 343 | 0.2 | 0.954 / 659 | 0.988 / 939 | 1.000 / 1,212 | 1.000 / 1,220 | 1.000 / 1,249 | 1.000 / 1,297 |
| 10k | 8,409 | 1.0 | 0.806 / 1,002 | 0.896 / 1,588 | 0.960 / 2,928 | 0.986 / 5,022 | 0.998 / 7,847 | 0.998 / 9,977 |
| 10k | 8,409 | 0.2 | 0.810 / 976 | 0.898 / 1,534 | 0.968 / 2,900 | 0.992 / 4,957 | 0.998 / 7,832 | 0.998 / 9,811 |
| 30k | 24,300 | 1.0 | 0.652 / 1,190 | 0.778 / 1,946 | 0.902 / 3,727 | 0.952 / 6,239 | 0.974 / 10,759 | 0.986 / 17,947 |
| 30k | 24,300 | 0.2 | 0.642 / 1,156 | 0.788 / 1,895 | 0.918 / 3,694 | 0.970 / 6,157 | 0.984 / 10,405 | 0.990 / 16,973 |

表中单元格为 "recall@10 / 平均查询延迟 µs"。
构建时间（ef_c=64，单线程，两种 scale 相同）：N=1k ≈ 0.8 s；
N=10k ≈ 28 s；N=30k ≈ 115–119 s（约 rust-cv 的 2×）。

## 结论

1. **HubNSW 单层化只有微弱正收益**：4096 维下两种 scale 召回几乎持平，
   中高 ef 单层略好（30k 时 +1.6 pp @ ef=50、+1.8 pp @ ef=100、
   +0.4 pp @ ef=400），低 ef（10）反而略低 1 pp；延迟完全一致。
   GSearch 官方"高维更准更省"的说法在 HV 上有微弱证据，但不是关键变量。
2. **实现差异比层数更关键**：同样数据下 `hnsw_rs` 召回远高于 rust-cv
   `hnsw` 0.11——N=30k 时 rust-cv 上限 recall@10=0.92（ef=400 +
   ef_c=200），`hnsw_rs` 达 0.974–0.990（ef=200–400）。代价是查询慢
   7–10×（30k：ef=200 时 10.4–10.8 ms vs 1.33 ms；ef=400 时 17–18 ms
   vs 1.51 ms）且构建慢 ~2×。因此 `hv.md` 里
   "4096 维 HNSW 召回天花板 <0.92"是 rust-cv 特定实现的结论，不是
   HNSW 算法本身的结论。
3. **速度增益大幅缩水**：N=30k 时 hnsw_rs 相对精确（24.3 ms）只快
   2.3–6.5×（ef=200 时 10.4 ms，2.3×；ef=50、recall 0.90–0.92 时
   3.7 ms，6.5×）；N=10k 时 ef=50（recall 0.96–0.97）快 2.9×。
   不再有 rust-cv 的 20–46×——hnsw_rs 用更重的候选计算换回了召回。
4. **对 pgr 的建议更新**（`design/genome-nn-query.md` §6.4）：
   - ≤10k：精确扫描（10k 时 8.4 ms/查询）仍然够用，结论不变；
   - 10k–30k：若接受 0.90–0.97 召回，hnsw_rs 在 ef=20–50 给出
     1.5–3.7 ms/查询（比精确快 3–5×）；rust-cv 更快的 0.5–1.3 ms
     对应 0.85–0.92 召回——两者是"召回优先 vs 速度优先"的实现选择；
   - >30k：hnsw_rs 召回上限虽高（0.99）但查询已接近精确扫描量级，
     ANN 收益所剩无几；**降维（PCA 256–512）仍是首选**，可在低维下
     同时获得高召回与高速度。
5. dev-dependency `hnsw_rs`（连带 probminhash，GSearch 参考实现）只
   服务于本实验与后续 ANN 再评估，不进入主程序。

---

## 证据附录：HV 最近邻：知识路由（clade 分片）vs 全局 HNSW（4096 维）

> 目的：验证 `design/genome-nn-query.md` §6.5 的设想——用外部知识
> （系统发育 / 生物学性状）选代表节点做查询路由，能否改善 4096 维 HV
> 的 ANN 查询。日期：2026-08-08。
> 硬件：AMD Ryzen 9 7945HX（单线程，release 构建，LTO）。

## 方法

- **合成 cohort 带显式 clade 结构**：16 个 clade，每 clade 有独立
  2,048 个 k-mer 核心；每个"基因组" = 256 全局核心 + 2,048 clade 核心
  + 512–4,096 私有 k-mer → `hash_hv_bit` → L2 归一化 f32。clade 内
  相似度远高于 clade 间（共享 2,304 vs 256 个 k-mer）。
- **对比**（总构建成本相近）：
  * 全局 HNSW：`hnsw_rs` 0.3.4，M=24、ef_c=64、DistL2，N=10k/30k；
  * 路由 HNSW：16 棵 clade 内 HNSW（同一参数）；查询先按"到各 clade
    代表（clade 首个基因组）的精确点积"路由，取 top-R（R=1/2/4）
    clade 搜索，按距离合并取 top-10。
- **指标**：recall@10 对**全局精确 top-10**；路由准确率 = 查询真 clade
  落在代表 top-R 内的比例；查询延迟含路由成本（16 次点积）。
  复现：`cargo bench --bench hv_ann_clade`；环境变量
  `PGR_HV_ANN_SIZES` / `PGR_HV_ANN_EFS` / `PGR_HV_ANN_EFC` /
  `PGR_HV_ANN_CLADES`。

## 结果

路由准确率：两种规模下查询真 clade 都在代表 top-1（50/50）——合成数据
的先验是完美的，本实验测的是"完美先验下路由的上限收益"。

| N | 精确 µs/查询 | 变体 | ef=10 | ef=20 | ef=50 | ef=100 | ef=200 | ef=400 |
|---|---|---|---|---|---|---|---|---|
| 10k | 8,265 | global | 0.974 / 781 | 0.994 / 968 | 0.994 / 1,197 | 0.994 / 1,259 | 0.994 / 1,409 | 0.994 / 2,932 |
| 10k | 8,265 | routed R=1 | 0.970 / 681 | 0.988 / 859 | 0.988 / 988 | 0.988 / 1,029 | 0.988 / 1,099 | 0.988 / 1,158 |
| 30k | 23,988 | global | 0.822 / 926 | 0.840 / 1,300 | 0.884 / 2,102 | 0.908 / 2,533 | 0.908 / 2,758 | 0.908 / 3,185 |
| 30k | 23,988 | routed R=1 | 0.940 / 920 | 0.982 / 1,310 | 0.984 / 2,121 | 0.988 / 2,453 | 0.988 / 2,545 | 0.988 / 2,659 |
| 30k | 23,988 | routed R=2 | 0.940 / 1,408 | 0.982 / 2,128 | 0.984 / 3,346 | 0.988 / 4,041 | 0.988 / 4,633 | 0.988 / 5,220 |
| 30k | 23,988 | routed R=4 | 0.940 / 2,536 | 0.982 / 3,521 | 0.984 / 5,505 | 0.988 / 7,065 | 0.988 / 8,658 | 0.988 / 10,165 |

表中单元格为 "recall@10 / 平均查询延迟 µs"（R=2/4 在 10k 与 R=1 相同，
省略）。构建时间：10k 全局 8.8 s vs 16×0.4 s ≈ 6.4 s；30k 全局
41.8 s vs 16×2.3 s ≈ 37 s——路由方案总构建成本相近甚至略省。

## 结论

1. **知识路由有效，机制 = 有效 N 变小**：30k 时同 ef 下 routed 召回比
   全局高 12–16 pp（ef=10：0.940 vs 0.822；ef=20：0.982 vs 0.840；
   ef=50：0.984 vs 0.884），延迟几乎相同（ef≤20 时 routed 还略快）。
   全局在 30k 卡在 0.908（ef≥100），routed 卡在 0.988。
2. **收益随 N 增长**：10k 时全局 HNSW 在 ef≥20 已达 0.994，路由反而
   略低（0.988，受 clade 内小图自身召回限制）——≤10k 无路由必要。
3. **R>1 在完美先验下无收益**：R=2/4 召回完全不变，只增加延迟
   （搜索更多图）；R>1 的价值仅在先验不可靠时兜底。
4. **路由正确率是决定因素**：若先验错误，误路由查询的 top-10 全部落在
   错误 clade（只共享 256 全局核心，距离远大于真近邻），recall ≈ 0；
   期望召回 ≈ (1−m)·R₁（m = 误路由率），容错曲线近似线性。真实数据
   clade 间有共享内容，误路由不会全丢，但方向不变：先验质量直接决定
   收益，因此先验应作软路由（R>1 / 重排兜底），不宜硬性只搜一个 clade。
5. **对 pgr**：Necom 聚类 / 构树产出的分组与代表天然可作路由表，
   "clade 代表精确路由 → clade 内 HNSW / 精确"是 10k–30k 场景里比全局
   HNSW 更优的形态；配合 HGT / 质粒等先验不一致的兜底（R>1 或 ANN
   候选重排）即可落地。

---

## 证据附录：bac120 标记基因路由准确率（#19）

> 日期：2026-08-08。对应 `design/genome-nn-query.md` §7.4 #19。
> 目的：§6.5 的"生物学先验"落地方案——用保守标记基因做廉价查询路由，
> 测它在真实数据上的准确率。

## 方法

- 数据：135 cohort 基因组；`Domain/bac120.fa.gz`（369,814 个标记蛋白）
  与 `Domain/seq_asm_f3.tsv`（蛋白→组装→标记映射，旧名体系，按 GCF
  accession 映射到 cohort）。
- 路由键：**8 个全覆盖 bac120 标记**（PF01025.14、PF00466.15、
  TIGR00019、PF02576.12、TIGR00061、TIGR00054、TIGR00043、TIGR00029；
  135/135 覆盖）拼接成每基因组蛋白集。
- 距离：`pgr dist frac --protein --merge --scale 100`（默认 scale 1000
  对蛋白太稀，sketch 仅 1 元素——本身是个发现）。
- 评估：留一法——查询的最近 aa 邻居所属物种 = 路由结果；准确率 =
  路由物种 == 查询自身物种。对照：HV（pgi-to-hv 点积）路由、ANI
  top-1 物种（标签噪声上限）。

## 结果

| 路由方式 | 留一法准确率 |
|---|---|
| 8×bac120 标记蛋白（aa frac） | **0.756** |
| HV（pgi-to-hv 点积） | 0.822 |
| ANI top-1 物种（金标准上限） | 0.800 |

- 全局精确（HV 点积排序）recall@10 vs ANI = 0.711（n=135）。
- 路由后"只搜路由物种 clade"的 recall 无法统计（n=0）：路由目标物种
  clade 大多 <10 成员——再次印证 #7 的"小 clade 路由失效"。

## 结论

1. **8 个保守标记基因的 aa 距离即可路由 75.6% 的查询到正确物种**，
   接近 ANI 金标准上限（80%），只比全基因组 HV 路由低 7 pp——作为
   §6.5 的"生物学先验"，成本极低、效果可用。
2. 物种标签与 ANI top-1 只有 80% 一致（近缘株标注噪声/种内近邻跨
   物种），路由准确率的天花板就在 ~80%，标记路由已接近。
3. 蛋白 frac 默认 scale=1000 对短蛋白（~300 aa）过稀（sketch 1 元素，
   距离全 0/1），蛋白场景应降 scale（如 100）——建议 `dist frac
   --protein` 文档提示。
4. 路由 recall 仍需大 clade 场景验证（E. coli 全量 NR 是天然试验场）。

## 复现

```bash
# 1) cohort_marker.tsv：seq_asm_f3.tsv 按 accession 映射到 cohort
# 2) pgr fa some bac120.fa.gz <prots> → 按基因组拆分
# 3) pgr dist frac genomes.lst --list-files --merge --protein --scale 100
# 4) 留一法分析：/tmp/hv_calib/analyze_marker.py
```

---

## 证据附录：采样/参数扫描与组装质量对距离误差的影响（#3/#4）

> 日期：2026-08-08。30 基因组子集（10 E. coli + 其他 Escherichia +
> Pseudescherichia + Yersinia），435 对 vs skani ANI。对应
> `design/genome-nn-query.md` §7.4 #3/#4。

## #3 sampler/k/D 扫描（Spearman vs ANI，负值=距离对相似度）

| 变体 | 全部 | ≥99% | 95–99% | 90–95% | 85–90% | <85% |
|---|---|---|---|---|---|---|
| frac s=100 | −0.991 | −0.758 | −0.949 | −0.959 | −0.870 | −0.096 |
| **frac s=1000（默认）** | **−0.991** | −0.726 | −0.948 | −0.961 | −0.873 | −0.048 |
| frac s=10000 | −0.986 | −0.644 | −0.941 | −0.913 | −0.818 | −0.078 |
| mini k21w5（默认） | −0.886 | −0.614 | −0.471 | −0.679 | −0.454 | +0.126 |
| mini k31w5 | −0.894 | −0.614 | −0.560 | −0.707 | −0.368 | +0.114 |
| mini k15w5 | −0.859 | −0.588 | −0.433 | −0.567 | −0.518 | +0.080 |
| mini k21w10 | −0.887 | −0.588 | −0.480 | −0.682 | −0.450 | +0.105 |
| mini hasher fx / murmur | ≈ rapid（无差别） | | | | | |
| hv syncmer（s8/w5） | −0.866 | −0.616 | −0.922 | −0.688 | +0.182 | −0.196 |

RMSE(1−d, ANI/100)：frac s=1000 = 0.0115（≈1.15 个 ANI 点）；
mini k21w5 = 0.0414；hv syncmer = 0.0983。

结论：
1. **frac 全区间最优且默认 scale=1000 已接近上限**（s=100 仅近缘段
   略好 0.03，s=10000 近缘段变差）——默认参数合理，无需改。
2. mini 的 k/w/hasher 对 ANI 相关几乎无影响（k=31 中段略好）；
   minimizer 类采样近缘分辨率缺陷是结构性的（ρ≈−0.6）。
3. hv syncmer 在 95–99% 段（ρ=−0.92）明显好于 minimizer HV，但
   85–90% 段异常（ρ=+0.18），整体不稳定；不构成推荐切换。

## #4 组装质量（N50/大小/contig）对 frac ANI 误差的影响

frac s=1000 估计 ANI（1−d）与 skani ANI 的绝对误差：

| 层 | 误差中位数（ANI 点） | 与误差的 Spearman |
|---|---|---|
| 全部 | 0.96 | S_min −0.47；N50_min +0.18；C_min −0.16 |
| 同种内 | 0.23 | N50_min **−0.41**；C_min **+0.30**；size_ratio −0.23；S_min +0.11 |

结论：种内误差主要由**碎片化**驱动（N50 低、contig 多 → 误差大），
配对成员大小差异也增加误差——与 #2 完整度结论一致。装配质量差的
基因组应优先过滤或标注（种内聚类/标定场景）。

## 复现

```bash
# 30 基因组子集（subset30.lst），各变体 --list-files --merge：
# pgr dist frac --scale {100,1000,10000}
# pgr dist mini -k/-w 变体与 --hasher 变体
# pgr dist hv --sampler syncmer
# 分析脚本：/tmp/hv_calib/analyze_params.py、analyze_sizebias.py
```

---

## 证据附录：`dist mash` 与 Mash 2.3 字节级兼容对照

> 目的：验证 `pgr dist mash`（bottom-k MinHash）与参考工具 Mash（Ondov et
> al. 2016，本地 `Mash-master` + 系统 mash 2.3）在相同 k / sketch size 下
> 的距离输出完全一致。日期：2026-08-08。

## 算法核对结论（Mash-master 源码）

- canonical k-mer：正链与反向互补**字节比较取小**（`memcmp`），不是哈希取小；
- 哈希：`MurmurHash3_x64_128(kmer_bytes, seed=42)` 取前 64 位；
- 过滤：k-mer 含非 ACGT 直接跳过（整 k-mer 不参与）；
- bottom-k：全局保留 sketchSize（默认 1000）个最小唯一哈希；
- 距离：合并排序 → shared / denom（denom clamp 到 sketchSize）→
  Jaccard → `-ln(2J/(1+J))/k` clamp [0,1]。

**关键坑 1（Jaccard）**：Mash 的 Jaccard = `compareSketches` 的匹配数 /
sketchSize（合并两个排序 sketch、最多 sketchSize 步的匹配数），不是
标准集合 Jaccard。Mash-master/test 的 genome1×genome2 完整交集 581 个
哈希，但 Mash 报告 456/1000。`dist mash` 严格按 Mash 语义实现，
与 `mash dist` 的 shared/denom 完全一致（20/20 对，见对照 2）。

**关键坑 2（Containment，2026-08-08 发现并修复）**：`dist mash` 原实现
的 containment 用了与 Jaccard 相同的 merge-common（`common / a.len()`），
当集合满 sketchSize 时数值恰好等于 Jaccard——这不是标准 containment。
系统对照（5 株 × 20 对）显示它比标准值（完整集合交集 / query 集合大小，
Mash `within` 语义）系统性低估 0.12–0.18（相对偏差约 25–35%）。已修复：
containment 改为完整集合交集 / 第一个集合大小；修复后 20/20 对与标准值
完全一致（Δ=0.0000），Jaccard/距离不受影响。

**关键坑 3（Undersized sketch 的 denom，2026-08-08 审计修复）**：Mash
`compareSketches` 的 Jaccard 分母是 merge 遍历的 `denom`（一方提前耗尽时
补上剩余未遍历哈希，clamp 到 sketchSize），不是固定 sketchSize。原实现
除以 sketchSize，导致 sketch 未满（短 contig/质粒、高 k、小序列）时距离
错误——两个相同 46-hash sketch（k=15/s=1000）Mash 报 46/46、距离 0，
pgr 报 46/1000、距离 0.1621。已修复：分母与 `union` 输出均改为 Mash 的
denom 语义；修复后小 sketch 与 Mash 一致（2/2 → 0、0/92 → 1），满 sketch
不受影响。

**内存（2026-08-08 审计修复）**：草图构建原为全量物化（每条序列的完整
哈希 Vec + 缓冲、`--merge` 累积整个文件），内存 O(基因组长度)（4.6 Mb
单条 record 峰值 RSS 约 90 MB，100 Mb 基因组约 1 GB）。已改为滚动窗口
流式（O(k) 窗口）+ 增量 bottom-k 累积器（O(sketch_size)）：4.6 Mb 基因组
merge 模式峰值 RSS 降到 15.8 MB，且不随序列长度增长。流式实现与旧全量
逻辑逐哈希对照一致（单元测试 `test_for_each_mash_hash_matches_reference`）。

## 对照 1：Mash-master/test 三对（k=21, s=1000）

| 对 | mash dist | pgr dist mash | shared |
|---|---:|---:|---:|
| genome1 × genome2 | 0.0222766 | 0.0223 | 456/1000 |
| genome1 × genome3 | 0.0000000 | 0.0000 | 1000/1000 |
| genome2 × genome3 | 0.0222766 | 0.0223 | 456/1000 |

## 对照 2：5 株 E. coli 真实基因组 20 对（k=21, s=1000，整文件 `--merge`）

数据：`/tmp/pgr_cohort/data/{mg,sa,se,e2,cf}.fa.gz`（MG1655 / Sakai / SE11 /
E2348 / CFT073）。pgr 侧用 `--merge`（与 `mash dist` 的整文件语义一致）；
不加 `--merge` 时 pgr 是逐 contig 比较，行数会不同，注意区分。

| 对 | mash dist | pgr dist mash | diff |
|---|---:|---:|---:|
| mg×sa | 0.0227107 | 0.0227 | <1e-4 |
| mg×se | 0.0162749 | 0.0163 | <1e-4 |
| mg×e2 | 0.0304804 | 0.0305 | <1e-4 |
| mg×cf | 0.0312752 | 0.0313 | <1e-4 |
| sa×se | 0.0240597 | 0.0241 | <1e-4 |
| sa×e2 | 0.0328290 | 0.0328 | <1e-4 |
| sa×cf | 0.0350415 | 0.0350 | <1e-4 |
| se×e2 | 0.0339141 | 0.0339 | <1e-4 |
| se×cf | 0.0356220 | 0.0356 | <1e-4 |
| e2×cf | 0.0174761 | 0.0175 | <1e-4 |

全部 20 对一致，差异仅为 pgr 输出 4 位小数的舍入。结论：`dist mash`
与 Mash 字节级兼容，可作 Mash 的直接替代（同参数）。

## 性能

`benches/dist_sketch_benchmark.rs`（criterion，本机）：
`dist mash` 草图构建约 54 MiB/s（k=21 单线程、流式 bottom-k + 预筛，
2026-08-08 优化后；优化前 24 MiB/s）；对比 `dist mini` 约 104 MiB/s
（窗口采样，哈希次数少 5 倍）、`dist frac` 约 113 MiB/s（语义不同，
仅量级参考）。Mash C++ 单线程约 66 MiB/s（实测估算）。

### 并行扩展（2026-08-08，5 株 E. coli × 4 query，3 次取最快）

sketch 加载改为文件级并行后（`load_entries` rayon，三个草图命令共用）：

| -p | pgr dist mash（前） | pgr dist mash（后） | mash |
|---:|---:|---:|---:|
| 1 | 1.34 s | 0.57 s | 0.45 s |
| 4 | 1.34 s | 0.23 s | 0.18 s |
| 8 | 1.34 s | 0.24 s | 0.19 s |

单对（2×4.6 Mb）：pgr 0.24 s vs mash 0.18 s（1.3×）。剩余差距来自
纯 Rust murmur3 与 FASTA 读取（Mash 为 C++ 优化），未换哈希实现以保持
字节级兼容。`dist frac` 同样受益（0.70→0.31 s @ -p4）、`dist mini`
（0.53→0.30 s）。

### `dist frac` containment 同 k 对照（2026-08-08）

澄清此前"~10% Jensen 偏差"（`dist-cohort-validation.md`）：那是 frachash
k=21 vs 全 k=40 真值的 **k 不匹配假象**。同 k=21 全 canonical k-mer 集合
真值下，5 株 × 10 对：

| 对 | 真值 containment | frac (scale=1000) | 偏差 |
|---|---:|---:|---:|
| mg×sa | 0.6690 | 0.6591 | −1.5% |
| mg×se | 0.7474 | 0.7347 | −1.7% |
| sa×se | 0.5884 | 0.5840 | −0.7% |
| e2×mg | 0.4976 | 0.5012 | +0.7% |
| e2×sa | 0.5168 | 0.5049 | −2.3% |
| e2×se | 0.5021 | 0.5004 | −0.3% |
| cf×mg | 0.4819 | 0.4861 | +0.9% |
| cf×sa | 0.4875 | 0.4823 | −1.1% |
| cf×se | 0.4844 | 0.4795 | −1.0% |
| cf×e2 | 0.6680 | 0.6772 | +1.4% |

偏差正负对称、幅度 ~1–2%（scale=1000 时采样 SE≈0.7%），属采样方差而非
系统性偏差。Hera 2023 校正因子 (1−(1−s)^|A|) 对 |A|≥10⁵ ≈ 1，大肠杆菌
场景无校正效果；仅极短序列（|A|<~100）+ 大 scale 才有意义。`dist frac`
containment/ANI 保持现状。

---

## 证据附录：距离消费者 cohort 验证（dist pgi / dist hv / dist seq）

> 目的：验证 `.pgi` 的两个距离消费者（`dist pgi` 精确归并、`pgi to-hv` +
> `dist hv` 近似向量）在 10 株 E. coli 上的距离语义，与 `dist seq` 草图距离
> 及比对身份率（`pgr align pgi` 实测，见 [[../design/pgi-align.md]] §2.4）
> 对照。日期：2026-08-02 初测、2026-08-05 复测。

## 方法与数据

- 10 株 × 45 对：`dist pgi`（k=40、syncmer 8/5、精确集合归并）、
  `dist hv`（`pgi to-hv` dim=1024/4096/8192，k=40）、
  `dist seq`（syncmer k=8/w=5 草图）；
- 参考真值：45 对比对身份率（扩展 PSL 块）。

## 排序一致性（Spearman ρ，与 1 − 身份率）

### 初测（2026-08-02）

| 方法 | ρ | 备注 |
|---|---:|---|
| `dist seq`（k=8 草图） | **0.816** | 全基因组组成采样，最贴近真值 |
| `dist pgi`（k=40 精确） | 0.539 | 确定性但有偏（见下） |
| `dist hv`（dim 1024，稠密饱和） | **−0.05** | 饱和退化，不可用（见下） |

### 复测（2026-08-05）

| 方法 | ρ | 备注 |
|---|---:|---|
| `dist seq`（k=8 syncmer 草图，首 contig 对，jaccard） | **0.616** | 最贴近身份率真值 |
| `dist pgi`（k=40 精确，mash） | 0.590 | 确定性但有偏（见发现 1） |
| `dist hv`（稀疏 v2，dim 4096，mash） | 0.547 | 与 dist pgi 的 mash 排序一致（ρ=0.969） |

身份率真值 = pooled PSL identity（`(matches+rep)/block_len`，口径见
[[../design/pgi-align.md]] §2.2），2026-08-05 用当前二进制重算 45 对：
分布 0.9544-0.9879（均值 0.9682）。排序结论不变：`dist seq` 最贴近
身份率、`dist pgi` 次之、`dist hv` 与 `dist pgi` 高度一致。

### s=1 复测（2026-08-08，5 株 × 10 对）

`pgr pgi to-hv` 默认稀疏度 s 由 3 改为 1 后（`design/hv.md` §2.7 决策：
s 不影响精度、s=1 编码快 2.5×），用 5 株（MG1655 / Sakai / SE11 /
E2348 / CFT073，GCF 基因组）× 10 对验证：

| 对 | `dist hv` s=1 mash | `dist pgi` mash | 差 |
|---|---:|---:|---:|
| mg1655–se11 | 0.0126 | 0.0124 | +0.0002 |
| mg1655–sakai | 0.0180 | 0.0176 | +0.0004 |
| sakai–se11 | 0.0186 | 0.0186 | 0.0000 |
| mg1655–cft073 | 0.0256 | 0.0260 | −0.0004 |
| se11–cft073 | 0.0261 | 0.0270 | −0.0009 |
| sakai–cft073 | 0.0269 | 0.0278 | −0.0009 |
| e2348–cft073 | 0.0448 | 0.0438 | +0.0010 |
| sakai–e2348 | 0.0568 | 0.0557 | +0.0011 |
| mg1655–e2348 | 0.0569 | 0.0553 | +0.0016 |
| se11–e2348 | 0.0594 | 0.0569 | +0.0025 |

**Spearman ρ = 0.9879**（mash 排序），最大差异 0.0025——s=1 与精确
`dist pgi` 高度一致，排序与 s=3 时代一致（原 s=3 测 ρ=0.969，45 对）。
再次确认 s 不影响距离排序（`design/hv.md` §2.7 理论 + 实测）。

### FracMinHash 无偏性验证（2026-08-08，全集合对照）

新增 `dist seq --sampler frachash`（§5.5 方向 2）后，用**全 k-mer 集合
（不经采样）**验证各采样器的 Jaccard 估计（e2348 × cft073，k=40）：

| 方法 | Jaccard | 备注 |
|---|---:|---|
| **全 k=40 集合（真值）** | **0.4508** | 所有 canonical k-mer 精确交集 |
| **FracMinHash（k=40, scale=100）** | **0.4171** | 随机采样，**无偏**（差 0.034） |
| `dist pgi`（k=40 syncmer 8/5） | 0.0948 | **syncmer 位置偏差，严重低估 4.7×** |

**重要发现**：`dist pgi` 的 syncmer 8/5 采样在两个（近缘、含共享重复元件的）
基因组上有**位置偏差**——即使全集合重叠 45%，syncmer 位置因微小序列差异
而错开，交集比例远低于真值。FracMinHash（独立随机采样）无此问题，给出
接近全集合的 Jaccard。这也解释了此前 `dist pgi` "确定性但有偏"（§2 发现 1）
的量级。**FracMinHash 是 pgr 距离体系中唯一无偏的采样器**（排序类任务
各采样器大致可用，数值 ANI 必须用 FracMinHash）。

**偏差范围细化（2026-08-08，5 株 × 10 对 vs 全 k=40 真值）**：

| 组 | 对数 | 平均 ANI 偏差 |
|---|---:|---:|
| 不含 e2348 的对 | 6 | **~0.3 pp**（`dist pgi` 可靠） |
| 涉及 e2348 的对（EPEC，噬菌体/IS/质粒丰富） | 4 | **~3.3 pp**（数值不可信） |

`dist pgi` 不是"本身不准"：对"干净"基因组对偏差小（~0.3% ANI），
偏差集中在**重复/移动元件丰富的株**（syncmer 位置在重复区更难对齐，
交集低估被放大）。排序整体仍可用（偏差方向一致、按相似度缩放）。

**Jaccard vs Containment（2026-08-08，10 对 vs 全集合真值）**：

| 指标 | 平均相对偏差 | 排序 Spearman |
|---|---:|---:|
| Containment | **34.3%** | **0.661** |
| Jaccard（→ mash） | 39.3% | 0.515 |

containment 略准：分母是单集合（syncmer 采样率 ~57%），偏差被稀释；
jaccard 的分母 union 受 inter 低估放大。**但两者相对偏差仍 ~35%**、
排序 ρ 仅 0.5–0.7——syncmer 偏差对两个指标都有显著影响，`dist pgi`
的数值判断（无论 jaccard 还是 containment）均不可靠，无偏数值 ANI
必须用 `dist seq --sampler frachash`。

**FracMinHash 的 containment（2026-08-08，10 对 vs 全 k=40 真值）**：

| 方法 | 相对偏差 | 排序 Spearman |
|---|---:|---:|
| **FracMinHash containment** | **9.9%** | **0.976** |
| pgi containment（对照） | 34.3% | 0.661 |

FracMinHash containment = inter/card1（`dist seq --sampler frachash` 的
containment 列），**期望无偏**（共享 k-mer 判定一致，分子分母采样率 p
抵消）；残差 ~10% 是比值估计的 **Jensen 偏差**——实测 scale=10/100/1000
偏差几乎相同（一阶理论：Cov 与 E² 都 ∝ p²，p 抵消），**增大采样不能消除，
需偏差校正公式才能进一步精确**（Hera et al. 2023 方向）。排序 ρ=0.976
对粗筛/质粒检测足够；数值精确的 containment 留作校正公式后续。

> **更正（2026-08-08，同 k 对照实验）**：上述 ~10% 残差是 **k 不匹配的
> 假象**——frachash 用 k=21、真值用全 k=40 集合，两套 k-mer 的共享比例
> 本来就不同。同 k=21 下（5 株 × 10 对，全 canonical k-mer 集合真值）：
> containment 偏差仅 0.3–2.3%（正负对称，即采样方差，SE≈0.7%），
> jaccard 偏差 0.1–1.6%。读 Hera et al. 2023 原文后确认：其校正因子
> (1−(1−s)^|A|) 对 |A|≥10⁵ 的基因组 ≈ 1（4.6 Mb 时 ≈ e^(−4600)≈0），
> 仅对 |A|<~100 且 scale≥0.1 的极短序列有意义；Theorem 8 的 CI 需对 p
> 数值求解，与现有正态近似差异小。**结论：大肠杆菌/默认 scale 场景
> 无需 Hera 校正**，`dist frac` 的 containment/ANI 保持现状（数据见
> `design/hv.md` 同 k 对照节）。

### `dist hv --sampler syncmer` 偏差对照（2026-08-08，5 株 × 10 对）

此前 syncmer 偏差数据只覆盖 `dist pgi`（k=40 smer=8/w=5）。补上 `dist hv
--sampler syncmer`（smer=8/w=55，HV 投影）的对照，让 syncmer 偏差在两个
实验入口都可体验（hv 为 FASTA 直算、pgi 需先建索引）：

| 对 | hv syncmer jac | frac jac | 真值 k=21 jac | hv con | frac con | 真值 con |
|---|---:|---:|---:|---:|---:|---:|
| cf×e2 | 0.9467 | 0.5219 | 0.6680 | 0.9673 | 0.6772 | 0.5150 |
| cf×mg | 0.9175 | 0.3398 | 0.4819 | 0.9328 | 0.4861 | 0.3393 |
| cf×sa | 0.9190 | 0.3108 | 0.4875 | 0.9594 | 0.4823 | 0.3132 |
| cf×se | 0.9292 | 0.3180 | 0.4844 | 0.9605 | 0.4795 | 0.3217 |
| e2×mg | 0.9204 | 0.3484 | 0.4976 | 0.9394 | 0.5012 | 0.3462 |
| e2×sa | 0.9198 | 0.3246 | 0.5168 | 0.9652 | 0.5049 | 0.3301 |
| e2×se | 0.9216 | 0.3309 | 0.5021 | 0.9617 | 0.5004 | 0.3294 |
| mg×sa | 0.9269 | 0.4486 | 0.6690 | 0.9894 | 0.6591 | 0.4494 |
| mg×se | 0.9452 | 0.5471 | 0.7474 | 0.9947 | 0.7347 | 0.5535 |
| sa×se | 0.9281 | 0.4259 | 0.5884 | 0.9583 | 0.5840 | 0.4327 |

**发现**：hv syncmer 的 jaccard（0.92–0.95）远高于 k=21 真值（0.48–0.75），
看似"虚高"，但这不是 HV 方法的缺陷——**k=8 全集合真值本身 jaccard≈0.999**
（8-mer 空间 65536 基本饱和，近缘基因组 8-mer 几乎全重合），hv syncmer
相对 k=8 基线**低估 ~5–8%**（syncmer 位置漂移使共享 8-mer 丢失），与
`dist pgi`（k=40 锚定）的偏差**机制相同、方向一致**（都是交集低估），
只是 k=8 基线极高掩盖了量级。

**排序一致性**（10 对距离排序 Spearman）：

| 对照 | ρ |
|---|---:|
| hv syncmer vs 真值 | 0.636 |
| frac vs 真值 | 0.782 |
| hv syncmer vs frac | 0.612 |

结论强化：syncmer 家族（`dist pgi` 与 `dist hv --sampler syncmer`）数值
均不可信（pgi 严重低估、hv 因 k=8 饱和无区分度），排序也弱于 `dist frac`
（ρ 0.64 vs 0.78）。体验入口：`dist hv --sampler syncmer`（FASTA 直算）与
`dist pgi`（索引）；数值 ANI 一律 `dist frac`。

**syncmer 家族（dist pgi 与 dist seq --sampler syncmer 同机制）**
**（2026-08-08 补充）**：两者都是 closed syncmer 采样，位置偏差同源。
k=21 时 `dist seq syncmer` 比 frachash 低 ~2%（e2-cf 0.5075 vs
0.5191、mg-se 0.5451 vs 0.5548）；pgi（smer=8 锚定 k=40）更严重
（4.7×），因"小 s-mer 窗口锚长 k-mer"放大位置漂移。**数值距离
（jaccard/containment/ANI）一律用 FracMinHash**；syncmer 家族只承担
排序/粗筛（FASTA 侧 `dist seq syncmer`）与比对锚点（pgi/align）。
**分层使用建议**：**初筛用 `dist hv`**（.hv 路径，O(D) 固定比较——
大规模近邻/粗筛/聚类，hv.md §1.1 的定位）；候选对的**中等精度距离**
用 `dist seq --sampler frachash`（无偏 + CI）或 `dist pgi`（注意重复
元件株的偏差）；最终验证用 `align pgi`。

> **注意（2026-08-08 用户澄清）**：上述"位置偏差"是**采样集合交集**的
> 固有行为，不是 pgi 的缺陷——pgi 的 syncmer 是 **align pgi（FastGA
> 兼容基因组比对）的锚点基石**：closed syncmer 由 k-mer 序列自身判定
> （同步性），共享 k-mer 在两个基因组中确定性双选 → 锚点完整保留；
> FracMinHash 共享 k-mer 的双选概率是 1/scale（hash 相同、判定完全相关，
> 非 1/scale²）——但 scale=1000 时锚点只保留 1/1000（4.6Mb 约 4600 个）
> 且无窗口保证（随机分布、长 gap 无锚），链化不可用。**不能改 pgi
> 采样器**。正确分工：pgi + syncmer = 比对/排序；无偏数值 ANI =
> `dist seq --sampler frachash`（§5.5 方向 2，含 CI）。

### FracMinHash 排序验证（2026-08-08，5 株 × 10 对 vs 全集合真值）

用**全 k=40 canonical k-mer 集合**（不经任何采样）作为真值，对比
`dist seq --sampler frachash`（k=40、scale=100）的 mash 排序：

* **Spearman ρ = 1.0000**（10 对相对顺序完全一致）；
* FracMinHash mash 平均比真值高 0.0027（Jaccard 略低估）——比值估计的
  Jensen 偏差（E[X/Y] ≠ E[X]/E[Y]），scale 越大越小，不影响排序；
* 最近对 mg1655–se11、最远对 sakai–cft073 均与真值一致；e2348–cft073
  实际是第二近（J=0.45，真值），此前 `dist pgi` 的"远"是 syncmer 偏差。

## 三个发现

### 1. `dist pgi` 是"采样集合的精确距离"，不是全组成距离

两索引归并对共享 k-mer 集合的计数是**确定性的**（零采样方差），但集合本身
只含 syncmer 位置采样的 40-mer：点突变会改变局部 s-mer 哈希、进而改变
syncmer 选择，两侧位置集漂移远快于真实组成差异。MG1655 vs Sakai（97.3%
身份率）的 k-mer containment 只有 0.54、jaccard 0.33，mash 0.0175 高估
分歧。Spearman 0.590 说明排序大致跟随但受株系重复内容干扰
（如 ec2011c_3493–se11 身份率 97.9% 却给出 mash 0.038）。

### 2. `dist hv` 对大规模集合饱和退化（dim 无关）

`hash_hv_i8` 下每个 k-mer 种子向**所有维度**累加随机 i8，260 万种子使各维
值 ~±1.3e6；共享种子主导点积，`dot ≥ card` 后被 `calc_distances` 截断，
近缘株 containment 全部饱和为 1.0、mash 压缩 ~10×（0.0002–0.0046），
与身份率 Spearman ≈ 0。dim 1024 → 4096 → 8192 不改善（期望与 dim 无关）。
**结论**：`pgi to-hv` + `dist hv` 的"4 万级粗筛"定位需要重做 HV 距离
（如阈值化双极编码 + 余弦/Hamming，或按种子数缩放），当前实现不可用于
近缘基因组排序。

> **已修复（2026-08-02，复测 2026-08-05）**：改为**稀疏投影 + 余弦**
> （`.hv` v2）：
> 每个 k-mer 只更新 `--sparse`（2026-08-08 起默认 1，原默认 3）个随机
> 维度（±1），文件头存
> `sparse` 与 `n_kmer`；比较时余弦相似度估计共享数
> `inter = cos × √(n1·n2)`。默认 `--dim 4096`。
> 结果：与 `dist pgi` 的 mash **Spearman 0.969**（当前），
> 共享 k-mer 计数平均误差 2.39%（最大 12.94%，复测与初测 2.4%/13% 一致）；
> 45 对比较 8.25s → **0.12s（复测，~71×）**（初测 14.4s → 0.29s，50×）。
> 稠密编码（`hash_hv_i8`，dim 1024–8192 均饱和）与余弦改算均无效
> （余弦范围 0.994–0.999，信息已淹没），稀疏是有效修复。

### 3. 粗筛分层（2026-08-08 更新）

* **初筛（超大规模近邻过滤）**：`dist hv`（.hv 稀疏路径，**O(D) 固定
  比较**，s=1 已验证排序 ρ=0.988 vs `dist pgi`，§1）或 FASTA 侧
  `dist seq`（k=8 syncmer，Spearman 0.616 vs 身份率，快且稳健）。
* **二次距离（候选对）**：`dist seq --sampler frachash`（无偏 + CI）或
  `dist pgi`（已建索引时，注意重复元件株的偏差，见上）。
* **最终验证**：`align pgi`（精确比对）。

> **更新（2026-08-02，复测 2026-08-05）**：`.hv`（稀疏）现在是 `dist pgi`
> 的 ~70× 快、0.97 排序一致的近似层；`dist seq` 仍是与身份率最贴近的
> 草图层（当前 0.62 vs 0.59），因为 k=40 syncmer 集合受采样位置漂移限制。

## 附带交付：`dist hv` 支持 .hv 文件

按设计文档 §7.2 补齐了消费者链：`pgr dist hv a.hv b.hv` 直接比较
`pgi to-hv` 的输出（校验 k/dim 一致，复用 `calc_distances` 核心），
含集成测试（自比较 mash=0/jaccard=1、维度不匹配报错）。

## 相关文档

- 索引格式与消费者规划：[[../design/pbit.md]]（.pgi 距离消费者层级）
- 比对身份率真值：[[../design/pgi-align.md]] §2.4
