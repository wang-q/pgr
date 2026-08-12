# pgr dist

`pgr dist` 模块提供序列和向量的**距离/相似度计算**功能。它是构建系统发育树（Phylogeny）和聚类分析（Clustering）的核心前置步骤。

## 核心定位

- **定位**：多模式距离计算器（序列 / 索引 / 超向量）。
- **输入**：FASTA/蛋白序列文件、`.pgi` 索引、`.hv` 超向量文件。
- **输出**：Pairwise TSV 格式（`Name1 Name2 Distance ...`），可用于下游分析或矩阵构建。
- **互补**：
  - 上游：`pgr fa` (序列格式处理), anchr `fq` (reads 处理),
    `pgr fa six-frame` (蛋白翻译)。
  - 下游：聚类/构树工具。

## 核心概念（Concepts）

距离计算的输入是序列，但直接比较全部内容太慢，所有子命令都先做一件事：
**采样（sketch）**——从序列中挑出一小部分"代表性 k-mer"组成一个紧凑的
指纹，再比较指纹。代价是结果变为近似（由采样方式决定偏差与速度的权衡）。

**k-mer**：序列中长度为 k 的连续子串，是序列比较的基本单位。`-k` 越大，
子串越独特（区分度越高），但对变异越敏感。

**采样器（各子命令的区别所在）**：

- **MinHash（`dist mash`）**：对所有 k-mer 哈希后，保留最小的若干个
  （bottom-k）。两个指纹共享的 MinHash 比例估算集合相似度。Mash 工具
  即此法，`dist mash` 与其输出逐字节一致。
- **Minimizer / 最小子（`dist mini`）**：每个滑动窗口里取哈希最小的
  k-mer。采样密集、速度最快，适合排序/粗筛；但估计有偏，数值 ANI 不可用。
- **FracMinHash（`dist frac`）**：每个 k-mer 独立按 1/scale 概率保留。
  是**无偏**随机子集，Jaccard/containment 可直接换算 ANI（推荐数值用）。
- **Syncmer / 同步子（`dist hv --sampler syncmer`、`dist pgi`）**：当窗口
  内**最短的子串（s-mer）**出现在窗口端点时才采样。保证覆盖但稀疏；
  缺点是位置随突变漂移，距离估计有偏（仅作实验/排序，见下）。

**度量（输出列的含义）**：

- **Jaccard index**：交集 ÷ 并集，衡量两个指纹重叠程度。
- **Containment**：交集 ÷ 第一个集合大小，衡量 B 覆盖 A 的比例（有方向）。
- **Mash distance / ANI**：由 Jaccard（或 containment）按模型换算的进化
  距离；`1 − distance` 近似平均核苷酸一致性（ANI，百分比）。

**Hypervector / HV（`dist hv`）**：把序列的 k-mer 投影到固定维度（默认
4096）的高维随机向量并叠加成一个向量。比较两个向量（点积/余弦）即可
近似距离，与序列长度无关、速度极快，适合超大规模初筛；精度低于直接
集合比较。

> 详细数学与文献（采样器偏差、FracMinHash 无偏性、syncmer 位置漂移等）
> 见 `notes/design/hv.md` 与 `notes/references/hv.md`。

## 子命令详解

草图距离命令族（`mini` / `mash` / `frac`）界面统一：输入（单/双文件、
`--merge`、`--list-files`）、并行（`-p`）、输出格式（`Name1 Name2
[Total1 Total2 Inter Union] Mash Jaccard Containment`）、`--sim`/`--zero`
语义一致；区别在采样算法与数值语义（注意：`mash` 的 `Inter`/`Union` 是
Mash `compareSketches` 的 shared/denom 语义，未满 sketch 时 `Union` 小于
`--size`；`mini`/`frac` 的 `Union` 是标准集合并集）。

### 1a. `pgr dist mini`: Minimizer 草图距离（排序/粗筛）
*窗口内最小哈希 k-mer 采样，速度最快，适合大规模排序。*

- **算法**：窗口 $w$ 内取哈希最小的 $k$-mer（minimizer）；DNA 默认
  `-k 21 -w 5`、蛋白默认 `-k 7 -w 1`。
- **参数**：`-k`/`--kmer`、`-w`/`--window`、`--hasher`
  （`rapid`/`fx`/`murmur`/`mod`）、`--protein`。
- **注意**：minimizer 的 Jaccard 估计**有偏且不一致**（Belbasi et al.
  2022）——本命令用于快速排序/粗筛，**数值 ANI 用 `dist frac`**。

### 1b. `pgr dist mash`: Mash 兼容 MinHash 草图距离
*bottom-k MinHash sketch，与 Mash（Ondov et al. 2016）字节级兼容。*

- **算法**：canonical k-mer（正/反链字节比较取小）→
  MurmurHash3_x64_128（seed 42）→ 保留 `--size` 个最小唯一哈希；
  Jaccard = Mash `compareSketches` 语义（合并两个排序 sketch、最多
  `--size` 步的匹配数 / merge denom——sketch 未满时 denom 补上剩余
  未遍历哈希并 clamp 到 `--size`，不是标准集合 Jaccard）；
  Containment = 完整集合交集 / 第一个集合大小（Mash `within` 语义，
  以第一个输入为 query）。
- **参数**：`-k`（默认 21）、`--size`（默认 1000，Mash 默认）、
  `--seed`（默认 42，Mash 默认）。
- **兼容性**：与 `mash dist` 相同 k/size 时输出一致（已验证：
  E. coli MG1655×Sakai k=21/s=1000 → 456/1000、距离 0.0222766；
  完整 20 对真实基因组对照见
  `notes/benchmarks/bench-dist-mash-compat.md`）。注意默认输出过滤差异：
  Mash 默认输出 shared=0 的行（距离 1），pgr 家族默认过滤、加 `--zero`
  才输出——行内容一致，默认行集合不同。
- **注意**：MinHash Jaccard 对等大小集合无偏（Broder 1997）；containment
  以第一个输入为分母、有方向性，对大小悬殊集合有偏。无 `--ci`
  （Jaccard 语义与 FracMinHash 不同）。

### 1c. `pgr dist frac`: FracMinHash 草图距离（无偏数值 ANI）
*每个 canonical k-mer 独立以 1/scale 概率保留——Jaccard/containment 无偏，
推荐用于数值 ANI。*

- **算法**：FracMinHash（Irber et al. 2022），`hash < u64::MAX/scale`。
- **参数**：`-k`（DNA 默认 21、蛋白 7）、`--scale`（默认 1000，越小越密、
  方差越低）、`--ci`（追加 ANI 的 95% 置信区间，Hera et al. 2023 提供
  更紧的校正界）。
- **验证**：与全集合真值排序 Spearman 1.0、Jaccard 无偏（详见
  `notes/benchmarks/dist-cohort-validation.md`）。

### 2. `pgr dist hv`: 基于 Hypervector 的序列距离
*利用超维计算（Hyperdimensional Computing, HDC）技术，将序列采样（minimizer 或 closed
syncmer）映射为固定维度的向量。*

- **核心概念**:
  - 将 k-mer 映射为高维空间（如 4096 维）中的随机向量。
  - 通过向量叠加（Superposition）表示整条序列。
  - 具有全息特性，对噪声鲁棒，且计算速度极快（位运算）。
- **优势**:
  - 维度固定，计算复杂度与序列长度无关。
  - 适合超大规模数据集的快速预筛选。
- **`.hv` 文件模式**: 输入为 `.hv` 文件（`pgr pgi to-hv` 产物）时直接比较，
  无需重新采样序列；要求两侧 `k`、维度与稀疏更新数（`--sparse`）一致。
  比较用余弦相似度恢复共享 k-mer 数（`inter = cos × √(n1·n2)`），是
  `pgr dist pgi` 的约 50× 快近似（排序相关性 ρ≈0.97）。
- **参数**（与草图命令族共享 sampler/hash 参数，含义与默认值相同）:
  - `--dim`: 向量维度 (默认 4096，建议为 32 的倍数以利对齐/性能，实现不强制)。
  - `--sampler`: `minimizer` (默认) 或 `syncmer`（syncmer 默认 DNA `-k 8 -w 55`、
    蛋白 `-k 7 -w 5`）。**syncmer 的 `-k` 是 s-mer 长度**（判定窗口端点最小的
    短子串），`-w` 是窗口内 s-mer 数，实际采样跨度 = k+w−1 碱基（如 DNA
    默认 k=8/w=55 → 跨度 62），与 mini/mash/frac 的 k-mer 语义不同；
    `--sampler minimizer` 时 `-k`/`-w` 与草图家族一致。`syncmer` 作为实验
    选项保留：closed syncmer 采样
    非均匀，Jaccard/containment 估计有偏（详见
    `notes/benchmarks/dist-cohort-validation.md`），用于与 minimizer /
    FracMinHash 结果对照体验偏差；数值 ANI 用 `dist frac`。
  - `--hasher`: 哈希算法（`rapid`/`fx`/`murmur`/`mod`，默认 `rapid`）。仅对
    `minimizer` 采样生效；`syncmer` 采样忽略 `--hasher`（DNA 用 2-bit 规范化
    滚动哈希、蛋白用 RapidHash）。
  - `-k`/`--kmer`: k-mer 长度 (默认 7)。
  - `-w`/`--window`: Minimizer 窗口大小 (默认 1)。
  - `--sim`: 将 Mash 距离转为相似度输出。
  - `--list-files`: 将输入视为文件列表。
  - `-p`/`--parallel`: 并行线程数。

### 3. `pgr dist pgi`: 基于 .pgi 索引归并的精确距离
*将两个 `.pgi` 索引的排序 k-mer 流线性归并，计算确定性的 Jaccard/Containment/Mash 距离。*

- **输入**：两个 `.pgi` 索引（`pgr pgi build` 生成）。
- **算法**：两排序流线性归并（O(|K1|+|K2|)），共享 k-mer 集合精确计数。
- **要求**：两侧索引采样参数（k/syncmer/window）必须一致，否则报错。
- **输出格式**：`Name1 Name2 Total1 Total2 Inter Union Mash Jaccard Containment`。
- **Containment**：与草图命令族一致，以**第一个**索引为分母、有方向性
  （`inter / total1`），交换两个参数会改变其值。
- **注意**：该距离是"采样集合的精确距离"（确定性、零采样方差），但 syncmer
  采样位置随变异漂移（重复元件株数值偏差可达 ~3% ANI；containment 略稳于
  Jaccard）；排序大致可用。数值 ANI 用 `dist frac`；初筛用 `dist hv`。
  `dist pgi` 是 align 生态的辅助消费者（索引约 FASTA 的 27×，为比对而建）。

## 典型用法

### 场景 A：基因组快速比较（Mash 兼容）
```bash
# 与 Mash 输出一致（k=21、size=1000）
pgr dist mash genome1.fa genome2.fa --merge -k 21 --size 1000

# 输出: File1 File2 Total1 Total2 Inter Union Mash Jaccard Containment
```

### 场景 B：无偏数值 ANI
```bash
# FracMinHash + 95% ANI 置信区间
pgr dist frac genome1.fa genome2.fa --merge -k 21 --scale 1000 --ci

# 快速排序/粗筛（minimizer）
pgr dist mini genes.fa -p 4 > dist.tsv
```

### 场景 C：基因组文件间 Hypervector 比较
```bash
# 将每个文件的所有序列合并为一个 hypervector，计算两个基因组之间的距离
pgr dist hv genome1.fa genome2.fa
```

## 如何选择距离方法（2026-08-08 分层建议）

| 需求 | 推荐 | 理由 |
|---|---|---|
| **超大规模初筛**（近邻过滤/聚类候选） | `dist hv`（`.hv` 路径） | O(D) 固定比较、`.hv` 仅 FASTA 的 1/87，最快的粗筛层 |
| **无偏数值 ANI** | `dist frac`（`--ci` 输出置信区间） | FracMinHash 独立随机采样，Jaccard 无偏（与全集合真值排序 Spearman 1.0） |
| **Mash 兼容距离** | `dist mash` | 与 `mash dist` 输出一致（k/size 相同），生态对照/迁移 |
| **快速排序/粗筛** | `dist mini` | minimizer 草图，最快；数值有偏（仅排序用） |
| **已建 `.pgi` 时的距离** | `dist pgi` | 零额外 I/O，但注意 syncmer 采样偏差（重复元件株 ~3% ANI；containment 略稳于 Jaccard）；排序大致可用 |
| **最终验证/排序** | `align pgi` | 精确比对（索引是 align 的必需品，`dist pgi` 只是其辅助消费者） |

**要点**：`.pgi` 索引约 FASTA 的 27×（为链化比对设计），**距离计算不应为它建索引**；初筛用 `.hv`、数值 ANI 用 FracMinHash。详见
`notes/benchmarks/dist-cohort-validation.md`。

## 未来规划 (Roadmap)

通用向量距离度量（euclid/cosine/jaccard）目前仅在 `libs::linalg` 中提供，
`pgr dist vector` CLI 已于 2026-07 移除，暂无恢复规划。
