# cutadapt: Mott/BWA 式质量修剪与接头去除

> 整理于 2026-08，源自对 `cutadapt-main/`（v5.2, 2025-10-23）源码的分析。目的：承接 [sickle.md](./sickle.md) 中"现代质量修剪算法对比"的结论——滑窗仍是默认，cutadapt 的 `-q` 提供了**唯一真正不同的算法方向**（Mott/BWA 累积质量法）。本文聚焦**按质量分数修剪**（非接头去除），为 pgr 的 `fq trim-q` 提供算法与 CLI 参考。接头/适配器部分仅作背景，pgr 当前只关心质量修剪。

## 1. 简介

`cutadapt`（Marcel Martin, 2011）是 FASTQ 预处理的事实标准工具，以**接头/引物去除**见长，同时提供质量修剪、长度过滤、poly-A 修剪、expected-errors 过滤等。本文只分析其质量修剪。

- **质量修剪入口**：CLI 的 `-q, --quality-cutoff [5'CUTOFF,]3'CUTOFF`（见 [cli.py](file:///home/wangq/Scripts/pgr/cutadapt-main/src/cutadapt/cli.py#L268)）。
- **核心实现**：`quality_trim_index`（[qualtrim.pyx](file:///home/wangq/Scripts/pgr/cutadapt-main/src/cutadapt/qualtrim.pyx#L22)），Cython 加速，单条读段 O(n)。
- **算法来源**：注释明确说明与 **BWA 的 `bwa_trim_read`** 相同（`qualtrim.pyx:29-33`）。这与经典 Mott 算法同源（BWA `-q` 即 Mott 算法），累计和取最小。
- **变体**：`nextseq_trim_index`（NextSeq polyG 暗循环）、`poly_a_trim_index`（poly-A/poly-T）、`expected_errors`（Edgar 2015 期望错误数）。

> **范围说明**：pgr 若实现 `fq trim-q`，应移植 `quality_trim_index`（Mott/BWA 累积质量法）作为 `--method mott`，与继承自 sickle 的滑窗（`--method sliding`）并存，见 §4。

## 2. 核心算法：`quality_trim_index`（Mott/BWA 累积质量法）

这是 cutadapt 质量修剪的核心，约 50 行 Cython。**5' 端与 3' 端分别独立计算切点**。

### 2.1 算法步骤（`qualtrim.pyx:22-73`）

```
参数：qualities（ASCII 质量串）、cutoff_front、cutoff_back、base（默认 33）

5' 端：
  s = 0; max_qual = 0; start = 0
  for i in 0..n:
    s += cutoff_front - (qual[i] - base)     # cutoff 减实际质量
    if s < 0: break                          # 累积和为负，停止
    if s > max_qual: max_qual = s; start = i+1

3' 端（反向）：
  s = 0; max_qual = 0; stop = n
  for i in reversed(0..n):
    s += cutoff_back - (qual[i] - base)
    if s < 0: break
    if s > max_qual: max_qual = s; stop = i

if start >= stop: start, stop = 0, 0         # 无有效区间
return (start, stop)
```

**关键数学直觉**：定义 `score(i) = sum_{j<i} (cutoff - q_j)`（对 5' 端）。当累积和**首次由正变负**（`s < 0`）时停止；切点选在累积和**达到最大**（`max_qual`）之后的位置。等价于 BWA/Mott：在"前面质量足够好"的连续段内，找到累积质量相对 cutoff 的**局部最大点**作为修剪边界。

- **5' 端**：从 5' 到 3' 扫描，`s` 是 `cutoff - q` 的累积。`s < 0` 表示从此处开始质量持续低于 cutoff，应停止 5' 修剪；`s > max_qual` 时更新 `start = i+1`，即质量相对高点的下一个碱基。
- **3' 端**：从 3' 到 5' 反向扫描，对称逻辑，`stop = i`。
- **区间判定**：`start >= stop` 时整段无效，返回 `(0,0)`（即不修剪，或全部丢弃由上层决定——注意这里**不是**丢弃整条，而是置空区间，`QualityTrimmer` 会返回空 read）。

### 2.2 与滑窗（sickle/Trimmomatic）的本质区别

| 维度 | 滑窗（sickle） | Mott/BWA（cutadapt `-q`） |
| :--- | :--- | :--- |
| 判定单位 | 固定窗口内**平均**质量 | 逐碱基累积，**每碱基**都参与 |
| 切点 | 窗口首/末首个越阈值碱基 | 累积和局部最大点 |
| 中间低质量区 | 无法处理（只看两端窗口） | **可定位并修剪中间低质量区** |
| 复杂度 | O(n)（差分滑动） | O(n)（单遍累积） |
| 阈值语义 | 平均质量阈值 | 逐碱基 cutoff（`cutoff - q` 的累积） |

**优势**：Mott 能处理"读段中部有局部低质量区"的情况，滑窗只能修两端。**劣势**：语义更抽象（累积和），不如"窗口平均质量"直观。

### 2.3 三个变体（`qualtrim.pyx`）

1. **`nextseq_trim_index`**（`-q` 的 NextSeq 变体，`--nextseq-trim`）：NextSeq 双色编码中"暗循环"（无颜色）通常被读成高质量 G，出现在读段 3' 端。算法与 `quality_trim_index` 的 3' 端相同，但把 **G 碱基的质量强制设为 `cutoff - 1`**（`qualtrim.pyx:108`），使其不贡献正累积，从而把 polyG 尾巴当低质量去掉。
2. **`poly_a_trim_index`**（`--poly-a`）：poly-A/poly-T 尾巴检测，'A'(或'T') 得 +1，其他碱基 −2，累计 score 最大处为切点；错误率上限 0.2（`errors * 5 <= len`），小于 3 的尾巴忽略。
3. **`expected_errors`**（`--max-ee`）：用 Edgar et al. (2015) 公式从 Phred 质量计算期望错误数 `sum(10^(-Q/10))`，用于按总错误数过滤（非修剪）。C 实现 `expected_errors_from_phreds`。

## 3. 质量修剪的修饰器封装与 CLI

### 3.1 `QualityTrimmer` 修饰器（`modifiers.py:840`）

cutadapt 采用"修饰器（modifier）"流水线架构——每个操作是一个 `SingleEndModifier`，对 `read` 返回 `read[slice]`：

```python
class QualityTrimmer(SingleEndModifier):
    def __call__(self, read, info):
        start, stop = quality_trim_index(read.qualities, self.cutoff_front, self.cutoff_back, self.base)
        self.trimmed_bases += len(read) - (stop - start)
        return read[start:stop]
```

- 统计修剪碱基数（`trimmed_bases`）供报告。
- 返回 `read[start:stop]`，即保留 `[start, stop)` 区间（**左闭右开**）。

### 3.2 CLI 参数（`cli.py:268`）

- `-q, --quality-cutoff [5'CUTOFF,]3'CUTOFF`：可指定单个（默认只修 3' 端）或 `5,3` 两个值。
- `--quality-base N`：Phred 偏移，默认 33（Sanger）。
- `-Q`：R2 的独立 cutoff（配对，默认继承 R1）。
- 双端时 `-q 5` 未给 R2 则 R2 复制 R1 的修剪器（`cli.py:1065`）。

**`parse_cutoffs`**（`cli.py:419`）解析 `"5"` → `(0,5)`（只修 3'），`"6,7"` → `(6,7)`。

### 3.3 执行顺序（pipeline）

cutadapt 的 read 修饰器按固定顺序执行，质量修剪在**接头修剪之后**。典型顺序：去接头 → 合并/质控 → 质量修剪 → poly-A 修剪 → 长度/错误过滤。pgr 移植时若只想做纯质量修剪，顺序不存在依赖问题。

## 4. 对 pgr 的启示：`fq trim-q`

### 4.1 移植核心：Mott/BWA 累积质量法

`quality_trim_index` 算法体量极小（约 50 行 Cython），Rust 移植预计 40–60 行，无复杂数据结构，是标准的"单遍累积 + 局部最大"：

```rust
/// Return `(start, stop)` of the good-quality segment (Mott/BWA bwa_trim_read).
/// `qual` is the ASCII quality string; `base` is the Phred offset (default 33).
fn quality_trim_index(qual: &[u8], cutoff_front: i32, cutoff_back: i32, base: u8) -> (usize, usize) {
    let n = qual.len();
    let score = |q: u8, cutoff: i32| cutoff - (q as i32 - base as i32);
    // 5' end
    let (mut start, mut s, mut max_qual) = (0, 0i32, 0i32);
    for i in 0..n {
        s += score(qual[i], cutoff_front);
        if s < 0 { break; }
        if s > max_qual { max_qual = s; start = i + 1; }
    }
    // 3' end (reverse)
    let (mut stop, mut s, mut max_qual) = (n, 0i32, 0i32);
    for i in (0..n).rev() {
        s += score(qual[i], cutoff_back);
        if s < 0 { break; }
        if s > max_qual { max_qual = s; stop = i; }
    }
    if start >= stop { (0, 0) } else { (start, stop) }
}
```

### 4.2 与滑窗共存（`--method`）

根据 [sickle.md](./sickle.md) §4.5 的结论，`fq trim-q` 建议提供两种质量修剪算法：

- `--method sliding`（默认）：忠实移植 sickle/Trimmomatic 滑窗语义，直观通用。
- `--method mott`：移植 `quality_trim_index`，可修中间低质量区，作为更精细的备选。

两者都以相似的方式结合 `--qual-threshold`/`--length-threshold` 与长度过滤。

### 4.3 值得借鉴的设计

1. **Cython/常数级优化**：cutadapt 用 Cython 手写质量修剪热路径（避免 Python 逐字符开销）。pgr 用 Rust 原生，天然高效，无需额外处理。
2. **5'/3' 独立 cutoff**：`-q 5,3` 允许两端不同阈值，比 sickle 的单阈值灵活。pgr 可支持 `--qual-front`/`--qual-back` 或 `--quality-cutoff FRONT,BACK`。
3. **`--quality-base` 可配置**：默认 33，但保留覆盖能力（Solexa +64 等）。
4. **NextSeq polyG 变体**：`nextseq_trim_index` 把 G 质量强制为 `cutoff-1`，简单优雅地处理 NovaSeq/NextSeq 的暗循环 G 尾巴。若 pgr 面向现代测序平台，值得作为 `--nextseq` 选项移植。
5. **期望错误数过滤**（`expected_errors`）：现代长读/高保真数据常用整条读的期望错误数而非逐碱基质量做过滤，可作为 `fq trim-q` 之外的 `fq` 过滤功能参考。

### 4.4 与 pgr 现有约束的应对

- `quality_trim_index` 对质量字符越界**不校验**（Cython 直接做 `qual[i] - base`，可能得负值），cutadapt 靠 `expected_errors` 里才校验。pgr 的硬约束是零 panic，移植时应在读取质量时校验字符落在 `[base, base+41]` 或报告错误，而非静默产生负质量。
- cutadapt 的模块是修改共享 `read` 对象（`read[start:stop]` 返回新切片），pgr 用 noodles 的 `SequenceRecord`，注意切片语义（左闭右开）一致即可。

---

*参考来源: [cutadapt GitHub](https://github.com/marcelm/cutadapt) | 本项目源码 `cutadapt-main/`（v5.2） | [sickle.md](./sickle.md)*