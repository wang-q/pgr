# cutadapt: Mott/BWA 式质量修剪与接头去除

> 整理于 2026-08，源自对 `cutadapt-main/`（v5.2, 2025-10-23）源码的分析。目的：承接 [sickle.md](./sickle.md) 中"现代质量修剪算法对比"的结论——滑窗仍是默认，cutadapt 的 `-q` 提供了**唯一真正不同的算法方向**（Mott/BWA 累积质量法）。本文聚焦**按质量分数修剪**（非接头去除），为 pgr 的 `fq trim-qual` 提供算法与 CLI 参考。接头/适配器部分仅作背景，pgr 当前只关心质量修剪。

## 1. 简介

`cutadapt`（Marcel Martin, 2011）是 FASTQ 预处理的事实标准工具，以**接头/引物去除**见长，同时提供质量修剪、长度过滤、poly-A 修剪、expected-errors 过滤等。本文只分析其质量修剪。

- **质量修剪入口**：CLI 的 `-q, --quality-cutoff [5'CUTOFF,]3'CUTOFF`（见 [cli.py](file:///home/wangq/Scripts/pgr/cutadapt-main/src/cutadapt/cli.py#L268)）。
- **核心实现**：`quality_trim_index`（[qualtrim.pyx](file:///home/wangq/Scripts/pgr/cutadapt-main/src/cutadapt/qualtrim.pyx#L22)），Cython 加速，单条读段 O(n)。
- **算法来源**：注释明确说明与 **BWA 的 `bwa_trim_read`** 相同（`qualtrim.pyx:29-33`）。这与经典 Mott 算法同源（BWA `-q` 即 Mott 算法），累计和取最小。
- **变体**：`nextseq_trim_index`（NextSeq polyG 暗循环）、`poly_a_trim_index`（poly-A/poly-T）、`expected_errors`（Edgar 2015 期望错误数）。

> **范围说明 / 落地状态**：本文初稿时 pgr 的 `fq trim-qual` 还是设计提案；现已实现于 `src/cmd_pgr/fq/trim_qual.rs` + `src/libs/fq/trim.rs`，`quality_trim_index`（Mott/BWA 累积质量法）作为 `--method mott`，与继承自 sickle 的滑窗（`--method sliding`）并存。§4 已从"设计提案"改写为"已实现对照 + 剩余借鉴点"。

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
2. **`poly_a_trim_index`**（`--poly-a`）：poly-A/poly-T 尾巴检测，'A'(或'T') 得 +1，其他碱基 −2，累计 score 最大处为切点。错误率上限 0.2 的校验是**位置相关的**：5' 端（polyT head）为 `errors * 5 <= i+1`、3' 端（polyA tail）为 `errors * 5 <= n-i`（`qualtrim.pyx:147,161`），即错误按"当前已扫到的尾巴长度"而非整条读长计。长度 < 3 的尾巴忽略（`best_index < 3` → 0；`best_index > n-3` → n）。
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

- `-q, --quality-cutoff [5'CUTOFF,]3'CUTOFF`：可指定单个（默认只修 3' 端）或 `5,3` 两个值。注：cutoff 为 `0` 时该端禁用修剪（`cli.py:1059` 的 `cutoff != "0"` 判断）。
- `--quality-base N`：Phred 偏移，默认 33（Sanger）。
- `-Q`：R2 的独立 cutoff（配对，默认继承 R1）。
- 双端时 `-q 5` 未给 R2 则 R2 复制 R1 的修剪器（`cli.py:1065`）。

**`parse_cutoffs`**（`cli.py:419`）解析 `"5"` → `(0,5)`（只修 3'），`"6,7"` → `(6,7)`。

### 3.3 执行顺序（pipeline）

cutadapt 的 read 修饰器按**固定顺序**依次对 read 生效，`make_pipeline_from_args` 依序 append（`cli.py:937-975`）。实际顺序是：

1. `--cut` 无条件切头尾（`UnconditionalCutter`）
2. `--nextseq-trim`（NextSeq polyG 修剪）
3. `-q/--quality-cutoff` 质量修剪（`QualityTrimmer`）
4. 去接头（`-a/-g/-b`，`AdapterCutter`）
5. `--poly-a` poly-A/poly-T 修剪
6. `-l/--length` 缩短（`Shortener`）
7. 重命名/后缀

**注意：质量修剪在去接头之前**（`cli.py:947-967`，`make_quality_trimmers` 先于 `make_adapter_cutter`），与直觉相反——先按质量把低质量区剪掉，再做接头比对，这样接头搜索不受低质量 3' 尾干扰。pgr 移植纯质量修剪时无此依赖，但若日后加 `--method` 组合（如质量修剪 + 去接头）应保持"先质量后接头"的顺序。

## 4. 对 pgr 的启示：`fq trim-qual` 与 `fq clean`

> **落地状态**：pgr 的 `fq trim-qual` 已实现（`src/cmd_pgr/fq/trim_qual.rs` + `src/libs/fq/trim.rs`），本文 §2 的 `quality_trim_index` 以 `Method::Mott` 落地为 `--method mott`。本节从初稿的"设计提案"改写为"已实现对照 + 剩余借鉴点"。

### 4.1 已落地：`mott_cut`（`libs/fq/trim.rs`）

`Method::Mott` 在 `trim_interval`（`libs/fq/trim.rs:263`）中调用 `mott_cut`（`libs/fq/trim.rs:192`），与 Cython 版 `quality_trim_index` 逐行对应，约 30 行，仍是"单遍累积 + 局部最大"：

```rust
fn mott_cut(qual: &[u8], base: u8, cutoff_front: f64, cutoff_back: f64) -> (usize, usize) {
    // 5' end: s += cutoff_front - (q - base); s<0 break; s>max_s => start=i+1
    // 3' end(reverse): 对称，s<0 break; s>max_s => stop=i
    if start >= stop { (0, 0) } else { (start, stop) }
}
```

与 Cython 原版的差异：

- **int → f64**：cutadapt 的 `-q` 是整数 cutoff；pgr 的 `--qual-threshold` 是 `f64`，score 累积用 f64，阈值语义更宽（可给小数）。数值上等价于整数版。
- **单阈值两端共用**：`trim_interval` 把 `--qual-threshold` 同时作为 front 与 back cutoff；`--no-fiveprime` 时 front 置 `0.0` 禁用 5' 修剪。**未实现** cutadapt 的 `5,3` 独立 cutoff（见 §4.3.1）。
- **零 panic 校验**：`process_record` 先调 `validate_quality`（`trim.rs:245`），把质量字符限制在 `[base, base+93]`，越界即 `bail!` 报错——正是初稿 §4.4 的建议，已落实。
- **区间判定一致**：`start >= stop → (0,0)`；`trim_interval` 再按 `--length-threshold` 决定丢弃整条。
- **`--polyg-right`（polyg_end, `trim.rs:227`）**：3' 端数连续 G 跑，达标才剪（见 §4.3.2 与 cutadapt 的差异）。

### 4.2 与滑窗共存（已实现）

CLI 已提供 `--method sliding`（默认）与 `--method mott`（`trim_qual.rs:71-77`），对应 §2.2 的两类算法，验证了 [sickle.md](./sickle.md) §4.5 的结论。两者共用 `--qual-threshold`/`--length-threshold` 与长度过滤。

### 4.3 其余值得借鉴但尚未落地的点

1. **5'/3' 独立 cutoff**：cutadapt `-q 5,3` 允许两端不同阈值；pgr `mott_cut` 目前单阈值（加 `--no-fiveprime` 整体禁 5'）。若需"前端高严格、后端低严格"，可加 `--qual-front/--qual-back` 或 `--quality-cutoff FRONT,BACK`。
2. **cutadapt 的 NextSeq polyG（`nextseq_trim_index`）未移植**：pgr 的 `--polyg-right`（`polyg_end`）是 **BBDuk 式简单实现**——从 3' 端数连续 G 跑，长度达标才剪；`fq clean` 的 `--trim-poly-g-right` 同理。cutadapt 则是把 G 的质量强制设为 `cutoff-1` 后**嵌入累积算法**，能处理"夹着少量错配/低质量非 G 段"的 G 尾。对"纯连续 G 尾"两者等价，对"含噪声的 G 尾"cutadapt 式更鲁棒。若 pgr 面向 NovaSeq/NextSeq 平台，可评估移植 nextseq 式。
3. **poly-A 位置相关错误率（`poly_a_trim_index`）**：cutadapt 按"已扫到的尾巴长度"限制错误率（`errors*5 <= i+1` / `n-i`），比"整条按错误数"更精细；pgr `clean` 的 `--trim-poly-a` 走 BBDuk `trimpolyA`，语义不同。
4. **`--max-ee` 期望错误数过滤未落地为独立选项**：pgr 的 `expected_errors`（`trim_adapter.rs:907`、`clump.rs:607`）是 BBTools 式（`sum(10^(-Q/10))`，`prob[128]` 查表），用于 dedup/clump 排序，**不是** cutadapt 的 Edgar 2015 用户级过滤。长读场景可把 `--max-ee` 作为 `fq filter` 的过滤选项参考。
5. **`--quality-base` 可配置**：pgr 的 `--quality-base`（`33/64/auto`）比 cutadapt 更强——除显式覆盖外，还通过 BBDuk flip-flop 启发式（`detect_quality_base`, `trim.rs:88`）自动探测编码。

### 4.4 与 `fq clean` 的关系（重要澄清）

pgr 目前有**两条质量修剪路径，算法来源不同**：

- `pgr fq trim-qual`：滑窗（sickle）+ Mott（cutadapt 本文算法），纯质量修剪，单/双端，可配 `--outfile-2/--outfile-single`。
- `pgr fq clean`：整体是 **BBDuk（bbduk.sh）** 移植（`libs/fq/trim_adapter.rs`），其 `--qtrim` 取 `r/l/rl/w/f`——这是 **BBDuk 的质量修剪模式**（`r`=右侧逐碱基、`w`=滑窗等，对应 `trimq`/`qtrim-window`），**不是** cutadapt 的 Mott。即 `clean` 的质量修剪不来自 cutadapt。

因此 cutadapt 对 pgr 的算法价值**已集中在 `trim-qual --method mott`**；接头/适配器去除部分 pgr 走 BBDuk k-mer 路线（`clean --ref` + `--k`），未采用 cutadapt 的 `-a/-g/-b` 适配器匹配（BWA/seed 半全局比对）。

### 4.5 边界与实现注意

- `quality_trim_index` 对质量越界**不校验**（Cython 直接 `qual[i] - base`，可能得负值），cutadapt 只在 `expected_errors` 里校验。pgr 的 `validate_quality` 已补齐（零 panic 硬约束）。
- 左闭右开切片 `read[start:stop]` 与 pgr `SeqRecord::sequence()[start..end]` 语义一致。
- cutadapt `-q 0` 禁用修剪（`cli.py:1059` `cutoff != "0"`）；pgr 以 `--no-fiveprime` / 阈值 0 近似，语义可对齐。
- pgr 用 `f64` 阈值而非整数，注意与 cutadapt 精确逐碱基复现时，`f64` 累积与整数累积在极端大 read 上可能有浮点尾差，但常规阈值下等价。

---

*参考来源: [cutadapt GitHub](https://github.com/marcelm/cutadapt) | 本项目源码 `cutadapt-main/`（v5.2） | [sickle.md](./sickle.md)*