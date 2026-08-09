# 替换 anchr trim 流水线中的 BBTools（分析）

> 2026-08 整理。用户最终需求：anchr 的 `templates/trim.tera.sh` 用 BBTools
> 做 read 清洗，结果满意但速度不满意，目标是用 pgr 替换 BBTools。
> 配套：[fq-trim-q.md](fq-trim-q.md)（已实现，覆盖 sickle 部分）、
> [seq-reader.md](seq-reader.md)（FAFQ/BGZF 基础设施）。

> **定位（2026-08 修正）**：`pgr fq trim-q` 只替换流水线中的 **sickle**（第 9 步，
> 多阈值质量/长度参数扫描），**不替换 BBTools 任何组件**。BBTools 8 步（第 1-8）
> 是另一个替换目标，逐项梳理见 §3。

## 1. 流水线拆解

`trim.tera.sh` 对输入（R1/R2 或单端）依次执行（每步有 `if [ ! -e ]` 缓存）：

| # | BBTools 工具 | 作用 | 参数要点 |
|---|---|---|---|
| 1 | `clumpify.sh` | read 排序（可选去重 `dedupe dupesubs=0`） | 为后续阶段加速 |
| 2 | `filterbytile.sh`（可选） | 按 flowcell tile 过滤低质量 read | `opt.tile==1` |
| 3 | `bbnorm.sh`（可选） | kmer cutoff：滤掉低深度 kmer 的 read | `passes=1 bits=16 min=cutoff` |
| 4 | `reformat.sh`（可选） | 降采样 | `samplebasestarget=sample` |
| 5 | **`bbduk.sh trim`（核心）** | 接头 kmer 修剪 + 质量/长度修剪 | `ktrim=r k=trimk mink=11 hdist=1 tbo tpe maxns=0 minlen=60 qtrim=r trimq=15 ftm=5` |
| 6 | `bbduk.sh filter` | kmer 匹配过滤 adapter/artifact 库 | `k=matchk cardinality` |
| 7 | `kmercountexact.sh` | kmer 直方图 + peaks | `k=cutk` |
| 8 | `repair.sh` | 交错 → R1/R2/singles | `out/out2/outs` |
| 9 | `sickle` | 多阈值质量/长度修剪（参数扫描） | `-q/-l` 遍历，`Q{qual}L{len}/` 目录 |

第 9 步（sickle）已由 `pgr fq trim-q` 覆盖（anchr 的 sickle 调用只有
`-q/-l/-t sanger`，未用 `-n` 截断，trim-q 均已覆盖），见
[fq-trim-q.md](fq-trim-q.md)。注意：bbduk trim（第 5 步）里的
`qtrim/minlen/maxns/ftm` 是 **bbduk 的参数，不属于 trim-q 的替换范围**。

## 2. 速度瓶颈在哪

1. **8 个独立 Java 进程**：每个启动 JVM（秒级）+ GC + `-Xmx` 内存；pgr 是
   Rust 原生，启动毫秒级。
2. **中间文件反复 gz**：每步输出 `*.fq.gz`，下一步再解压——压缩/解压放大
   IO 多次。pgr 命令可流式管道衔接，中间不落盘。
3. BBTools 算法本身不慢（clumpify/bbduk 都是高效实现），慢的是**工程外壳**；
   但替换它意味着把其中几个算法在 pgr 里重新实现。

## 3. BBTools 8 步逐项梳理

| # | 工具 | 职责与算法本质 | pgr 基础 | 替换评估 |
|---|---|---|---|---|
| 1 | `clumpify.sh` | 按 kmer 签名排序/聚类 reads（可选 `dedupe` 去重）；为后续 kmer 类步骤加速 | 无 | 工作量大（TB 级外部排序）；若 pgr 后续步骤本身够快，其价值下降 → 四期 |
| 2 | `filterbytile.sh`（可选） | 按 Illumina flowcell tile 坐标过滤低质量 reads | 无 | 小众、老数据专属 → 可不做 |
| 3 | `bbnorm.sh`（可选） | kmer cutoff：全数据集 kmer 深度表（`bits=16` 近似哈希，省内存），滤掉含低深度 kmer 的 reads（去测序错误） | `KmerTable` 精确计数（u128+u32） | 中：内存策略需权衡（近似哈希 vs 精确表）→ 三期 |
| 4 | `reformat.sh`（可选） | 降采样到目标碱基数 | 无 | 小：`fq sample` 流式抽样 → 一期 |
| 5 | **`bbduk.sh trim`（核心）** | 接头 kmer 修剪（`ktrim=r`/`mink=11`/`hdist=1`/`tbo`/`tpe`）+ `maxns=0` + `qtrim=r` + `minlen` + `ftm=5` | 无（trim-q 只替代 sickle，不是这里） | **大**：`tbo`/`tpe` 是 BBDuk 特有 → 二期 |
| 6 | `bbduk.sh filter` | 参考库（adapter/artifact）kmer 匹配，过滤命中 reads + `cardinality` 统计 | `KmerTable`/`canonical_keys` | 中：kmer 命中过滤 → 三期 |
| 7 | `kmercountexact.sh` | 精确 kmer 计数直方图 + peaks（估计基因组大小/深度） | `KmerTable.counts` 可直接出直方图 | 小-中 → 三期 |
| 8 | `repair.sh` | 交错文件 → R1/R2/singles | 无（`fq interleave` 的反操作） | 小：`fq split` → 一期 |

## 4. 建议推进顺序

**一期（快赢，小命令）**
- `fq split`（交错 → R1/R2/singles，对齐 repair.sh）；
- `fq sample`（按碱基目标/比例降采样，对齐 reformat.sh）。

**二期（核心：bbduk trim 替换）**
- 接头 kmer 修剪：参考序列建 kmer 表（`k=trimk`、`mink=11`、`hdist=1`），
  `ktrim=r` 从 read 3' 端匹配并修剪；`tbo`（双端 overlap 检测接头）与
  `tpe`（成对均衡修剪）是 BBDuk 特有算法，需逐行移植
  （BBTools-40.01 `jgi/BBDuk.java` + `bbmerge`/`tbo` 实现）。
- 同命令整合 `maxns=0`、`qtrim=r`（复用 trim-q 的质量算法）、`minlen`、`ftm=5`。
- pgr 基础：`libs/kmer`（KmerTable、canonical_keys）+ SIMD 经验。

**三期（kmer 深度类）**
- bbnorm cutoff：KmerTable 建表 → 每 read 过滤低深度 kmer（内存策略待权衡）；
- bbduk filter：kmer 命中过滤 + 统计；
- kmercountexact：kmer 直方图 + 峰值。

**四期（独立项）**
- clumpify（read 排序/去重）、filterbytile（可跳过）。

## 5. 验收标准（替换后的对比）

- 同一输入，pgr 流程与 BBTools 流程的**输出 read 集合一致**（trim 步骤逐条
  比对 name + 序列；接头修剪部分允许 hdist 边缘差异，需用户确认容忍度）；
- 端到端墙钟时间显著下降（8 个 JVM 进程 → Rust 单进程流式）；
- 中间文件可选（管道模式不落盘）。

## 6. 待用户确认

1. 一期 `fq split` / `fq sample` 是否先做？
2. 接头修剪的语义对齐：完全复刻 bbduk（tbo/tpe/hdist）还是先做
   cutadapt 式简化版（无 tbo/tpe）？复刻工作量明显更大。
3. 中间文件策略：流水线用管道串联（pgr 命令流式、不落 gz）是否可接受？

---

*参考来源: [trim.tera.sh](../../../anchr/templates/trim.tera.sh)（anchr 项目，只读） |
BBTools-40.01 源码 | [fq-trim-q.md](fq-trim-q.md) | [seq-reader.md](seq-reader.md)*
