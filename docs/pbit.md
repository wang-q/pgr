# pgr pbit

> **状态**：pbit 格式 v1010；多参考、PAF 驱动 CIGAR/Identity 编码、
> `to-paf` 无损还原均已落地，`create`/`append` 强制要求 PAF
> （设计决策见 `notes/design/pbit.md`）。

`pgr pbit` 用于管理 **pbit**（population 2bit + delta）归档文件。pbit 是 pgr 原生的群体基因组压缩格式，
它将参考基因组以标准 2bit 记录存储，并将每个样本以 **PAF 驱动的 CIGAR delta 编码**压缩存储
（纯匹配段为零开销 Identity 引用）；PAF 未覆盖或无匹配的段用 LZ-diff/Raw 兜底，
适合在保留随机访问能力的同时压缩大量同源样本。

## 核心定位

- **定位**：群体基因组序列压缩与随机访问工具。
- **输入**：参考 FASTA 文件、样本 FASTA 文件（支持 plain 或 `.gz`），**PAF 比对文件（必需）**。
- **输出**：`.pbit` 二进制归档文件，或从归档中提取的 FASTA 文件。
- **互补**：
  - 上游：`pgr fa`、`pgr fa gz`（准备 FASTA），`minimap2`、`wfmash`（生成 PAF）。
  - 下游：`pgr paf`、`pgr fa range`（按坐标提取序列）。

## 子命令分组

*   **build**: 创建或追加归档。
    *   `create`: 从参考 FASTA 和样本 FASTA 创建新的 pbit 归档。
    *   `append`: 向已有 pbit 归档追加新样本。
    *   `append-ref`: 向已有 pbit 归档追加新参考基因组。
*   **info**: 查看归档元信息。
    *   `stat`: 显示归档统计信息、样本列表、参考 contig 列表或样本 contig 列表。
*   **subset**: 按名称或坐标提取序列。
    *   `range`: 按坐标区间从所有样本中提取序列片段。
    *   `some`: 按 contig 名称列表从所有样本中提取完整序列。
*   **transform**: 格式转换。
    *   `to-fa`: 将归档中的样本序列导出为每个样本一个 FASTA 文件。

---

## build Commands

### create

从参考 FASTA 和一个或多个样本 FASTA 创建新的 pbit 归档。样本与参考的 contig 名需匹配；
不匹配的样本 contig 会被跳过。

```bash
pgr pbit create [OPTIONS] -r <ref.fa> -i <sample.fa>... -o <out.pbit>
```

#### Options

*   `-r, --ref <file>`: 参考 FASTA 文件（必需，可重复指定多个参考，支持 plain 或 `.gz`）。
*   `-i, --infile <file>`: 样本 FASTA 文件，可多次指定（与 `--name` 互斥）。
*   `--name <file>`: TSV 文件，每行格式为 `sample_name<TAB>fasta_path<TAB>paf_path`
    `[<TAB>ref_name]`（第 3 列 PAF **必需**；第 4 列选择参考，名或序号，默认 0；
    与 `-i/--infile` 和 `--paf` 互斥）。
*   `-p, --paf <file>`: 与 `-i` 顺序对应的 PAF 文件（**必需**，与 `--name` 互斥）；
    传**空 PAF 文件**可禁用 CIGAR 编码（所有段走 LZ-diff/Raw 兜底）。
*   `-s, --segment-size <int>`: 参考序列分段大小，默认 4096 bp。
*   `-k, --kmer-len <int>`: LZ-diff 哈希 k-mer 长度，默认 15。
*   `-l, --min-match-len <int>`: LZ-diff 最小匹配长度，默认 18。
*   `-o, --outfile <file>`: 输出文件名（必需）。

#### Notes

*   样本名默认取自输入 FASTA 的文件名（使用 `--name` 可覆盖）。
*   仅支持 `ACGTN` 字符；IUPAC 简并碱基（R、Y、S 等）会被有损映射为 `N`。
*   参考与样本均**存储遮蔽**（soft mask，小写区间，语义同 2bit `mask_blocks`）：
   输入 FASTA 的小写区域会保留，`to-fa`/`some`/`range` 提取时原样还原小写，
   存进存出一致（v1005+）。
*   `--paf` 与 `--name` 互斥；如需为每个样本指定不同 PAF，请使用 `--name` 的第三列。
*   PAF 文件需包含 `cg:Z:` CIGAR，建议使用 `--eqx` 输出（如 `minimap2 -cx asm20 --eqx`）。
*   无 `cg:Z` 的 PAF 记录跳过 CIGAR 编码、原样存行还原（`to-paf` 仍可输出）。
*   从 MAF 建立归档（chainnet 输出为 MAF）：
    `pgr maf to-paf in.maf | pgr pbit create -r ref.fa -i sample.fa -p stdin -o out.pbit`
    （`-p stdin` 从管道读 PAF）。
*   实测边际 delta（相对 gzip-9 压缩的样本）：近缘样本 ≈ 50–57%，分歧
    （~90% ANI）≈ 81%。归档后建议跑 `pbit to-fa` 做覆盖率质量门。

#### Examples

1.  **创建单样本归档**:
    ```bash
    minimap2 -cx asm20 --eqx ref.fa sample.fa > sample.paf
    pgr pbit create -r ref.fa -i sample.fa -p sample.paf -o out.pbit
    ```

2.  **创建多样本归档**:
    ```bash
    pgr pbit create -r ref.fa -i s1.fa -p s1.paf -i s2.fa -p s2.paf -o cohort.pbit
    ```

3.  **通过 TSV 批量指定样本和 PAF**:
    ```bash
    cat samples.tsv
    # sample1	/path/to/s1.fa	/path/to/s1.paf
    # sample2	/path/to/s2.fa	/path/to/empty.paf
    pgr pbit create -r ref.fa --name samples.tsv -o cohort.pbit
    ```

4.  **多参考归档**（样本经 TSV 第 4 列路由到参考）:
    ```bash
    cat samples.tsv
    # s1	/path/to/s1.fa	/path/to/s1.paf	mg1655
    # s2	/path/to/s2.fa	/path/to/s2.paf	1
    pgr pbit create -r ref1.fa -r ref2.fa --name samples.tsv -o cohort.pbit
    ```

### append-ref

向已有归档追加新的参考基因组（2bit 段），旧样本全部保留。

```bash
pgr pbit append-ref <archive.pbit> -r <ref.fa> [-o <out.pbit>]
```

---

### append

向已有 pbit 归档追加新的样本 FASTA。参考序列已嵌入归档，因此不需要 `-r`。

```bash
pgr pbit append [OPTIONS] <archive.pbit> -i <sample.fa>...
```

#### Options

*   `<archive.pbit>`: 现有 pbit 归档文件（必需）。
*   `-i, --infile <file>`: 要追加的样本 FASTA 文件，可多次指定。
*   `--name <file>`: TSV 文件，格式同 `create`（与 `-i` 互斥）。
*   `-p, --paf <file>`: 与 `-i` 顺序对应的 PAF 文件（**必需**；空 PAF 可禁用 CIGAR）。
*   `-o, --outfile <file>`: 输出文件名。省略则原地修改输入归档。

#### Notes

*   省略 `-o` 时，追加通过临时文件 + 原子重命名完成，失败不会损坏原归档。
*   `-o` 指定的路径不能与输入归档相同；如需原地修改，请省略 `-o`。

#### Examples

1.  **原地追加样本**:
    ```bash
    pgr pbit append archive.pbit -i new_sample.fa -p new_sample.paf
    ```

2.  **追加到新的归档**:
    ```bash
    pgr pbit append archive.pbit -i s1.fa -p s1.paf -i s2.fa -p s2.paf -o new_archive.pbit
    ```

---

## info Commands

### stat

显示 pbit 归档的统计信息或列表。

```bash
pgr pbit stat [OPTIONS] <infile>
```

#### Options

*   `--samples`: 列出所有样本名。
*   `--refs`: 列出参考 contig 及其 segment 数量。
*   `--contigs`: 列出每个样本包含的 contig；结合 `-s` 可只列单个样本。
*   `-s, --sample <name>`: 与 `--contigs` 联用，只列出指定样本的 contig。
*   `-o, --outfile <file>`: 输出文件名（默认 stdout）。

#### Examples

1.  **显示归档概览**:
    ```bash
    pgr pbit stat archive.pbit
    ```

2.  **列出所有样本**:
    ```bash
    pgr pbit stat archive.pbit --samples
    ```

3.  **列出参考 contig 的 segment 数**:
    ```bash
    pgr pbit stat archive.pbit --refs
    ```

4.  **列出某样本的所有 contig**:
    ```bash
    pgr pbit stat archive.pbit --contigs -s sample1
    ```

---

## subset Commands

### range

按坐标区间从归档的所有样本中提取序列片段。每个匹配区间会输出一条 FASTA 记录。

```bash
pgr pbit range [OPTIONS] <infile> [ranges]...
```

#### Arguments

*   `[ranges]...`: 区间列表，格式为 `seq_name(strand):start-end`。
    *   `seq_name`: contig 名称（必需）。
    *   `strand`: 可选，`+`（默认）或 `-`。
    *   `start-end`: 1-based 闭区间坐标。
    *   只写 `seq_name` 表示提取整个 contig。

#### Options

*   `-r, --rgfile <file>`: 从文件读取区间，每行一个。
*   `-o, --outfile <file>`: 输出文件名（默认 stdout）。

#### Notes

*   坐标基于正链；指定 `-` 链时输出序列会被反向互补。
*   pbit 文件需要随机访问，不支持 stdin 或 `.gz` 输入。

#### Examples

1.  **提取单个区间**:
    ```bash
    pgr pbit range archive.pbit "chr1:1-1000" -o out.fa
    ```

2.  **提取多个区间**:
    ```bash
    pgr pbit range archive.pbit chr1:1-100 chr2:1-100 chr3 -o out.fa
    ```

3.  **从文件读取区间**:
    ```bash
    pgr pbit range archive.pbit -r ranges.txt -o out.fa
    ```

4.  **提取负链区间**:
    ```bash
    pgr pbit range archive.pbit "chr1(-):1-1000" -o out.fa
    ```

---

### some

按 contig 名称列表从归档的所有样本中提取完整序列。

```bash
pgr pbit some [OPTIONS] <infile> <list.txt>
```

#### Options

*   `-i, --invert`: 反向选择，输出不在列表中的 contig。
*   `-o, --outfile <file>`: 输出文件名（默认 stdout）。

#### Notes

*   列表文件每行一个 contig 名，空行和 `#` 开头的行被忽略。
*   名称匹配区分大小写。

#### Examples

1.  **提取列表中的 contig**:
    ```bash
    pgr pbit some archive.pbit list.txt -o out.fa
    ```

2.  **提取不在列表中的 contig**:
    ```bash
    pgr pbit some archive.pbit list.txt -i -o out.fa
    ```

---

## transform Commands

### to-fa

将归档中的所有样本序列导出为 FASTA 文件，每个样本一个文件。

```bash
pgr pbit to-fa [OPTIONS] <infile> -o <outdir>
```

#### Options

*   `-o, --outdir <dir>`: 输出目录（必需）。
*   `-s, --sample <name>`: 只导出指定样本。

#### Notes

*   输出文件为 `{outdir}/{sample_name}.fa`。
*   序列行宽固定为 60 bp。
*   样本名不能为空、不能包含 `/`、`\`，且不能为 `.`/`..`。

#### Examples

1.  **导出所有样本**:
    ```bash
    pgr pbit to-fa archive.pbit -o outdir/
    ```

2.  **只导出指定样本**:
    ```bash
    pgr pbit to-fa archive.pbit -o outdir/ -s sample1
    ```

### to-paf

将归档中内嵌的比对（PAF）导出为标准 PAF（12 列 + `cg:Z` 等标签）。

```bash
pgr pbit to-paf [OPTIONS] <infile> -o <out.paf>
```

#### Options

*   `-o, --outfile <file>`: 输出文件名（必需）。
*   `-s, --sample <name>`: 只导出指定样本的比对。

#### Notes

*   大链（主链）从归档 CIGAR 重建（`cg/cs/gi/bi/ms` 重算）；碎链和无
    `cg:Z` 记录原样还原输入行。
*   PAF 逐条可还原（"存进去什么，出来什么"）；坐标/链/CIGAR 均可访问。

#### Examples

1.  **导出全部样本的比对**:
    ```bash
    pgr pbit to-paf archive.pbit -o out.paf
    ```

---

## 典型工作流

### 场景 A：高质量参考 + 同源样本

```bash
# 1. 创建归档
minimap2 -cx asm20 --eqx ref.fa sample.fa > sample.paf
pgr pbit create -r ref.fa -i sample.fa -p sample.paf -o cohort.pbit

# 2. 查看统计
pgr pbit stat cohort.pbit

# 3. 导出为 FASTA
pgr pbit to-fa cohort.pbit -o outdir/
```

### 场景 B：从 MAF 建立归档（chainnet 链路）

```bash
# 1. 比对 + 链化（chainnet 输出 MAF），转 PAF 并管道进 pbit
pgr align pgi ref.fa sample.fa -o ref_vs_sample.psl
pgr pl chainnet ref.fa sample.fa ref_vs_sample.psl --t-name '' --q-name '' -o chain_out
pgr maf to-paf chain_out/*.maf | pgr pbit create -r ref.fa -i sample.fa -p stdin -o out.pbit

# 2. 追加更多样本
pgr pbit append cohort.pbit -i sample2.fa -p sample2.paf
```

### 场景 C：按坐标提取变异区域

```bash
# 提取目标区域在所有样本中的同源序列
pgr pbit range cohort.pbit "chr1:10000-11000" -o region.fa

# 进一步用于多序列比对或变异检测
pgr fas consensus region.fa
```

---

## 输入输出格式

### `--name` TSV 格式

`--name` 文件每行三到四列（第 3 列 PAF 必需，第 4 列可选），制表符分隔：

```text
sample_name<TAB>fasta_path<TAB>paf_path[<TAB>ref_name]
```

*   空行和 `#` 开头的行被忽略。
*   第三列为 PAF 文件路径（**必需**）；传空 PAF 文件可禁用 CIGAR 编码
   （该样本全部段走 LZ-diff/Raw 兜底）。
*   第四列（多参考时）选择该样本所属参考，可为参考名或序号；缺省时路由到参考 0。

### pbit 文件限制

*   二进制格式，需要随机访问（seek），不支持 stdin 或 `.gz` 输入。
*   仅支持 `ACGTN` 字符。
*   PAF 驱动的 CIGAR 编码不要求样本 contig 名与参考同名（按 PAF 记录映射）；
   无 PAF 覆盖的段走 LZ-diff 内容匹配 / Raw 兜底。

## 参考

*   设计与实现细节：`notes/design/pbit.md`
*   PAF 驱动模式背景：`notes/paf-pangenome.md`
