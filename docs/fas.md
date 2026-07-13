# pgr fas

`pgr fas` 提供了一套用于操作 **block FA** 文件的工具。block FA 是一种表示多序列比对（MSA）的格式，其中每个“block”包含来自不同物种或基因组区域的比对序列。

## 子命令

子命令按以下类别组织：

*   **信息（Info）**：从 block FA 文件中提取信息或统计量。
    *   `check`：根据参考基因组 FA 文件检查基因组位置。
    *   `cover`：计算染色体上的覆盖区域。
    *   `link`：提取双边或多边区间链接。
    *   `name`：列出文件中出现的物种名。
    *   `stat`：计算比对统计量（长度、差异等）。
*   **子集（Subset）**：筛选并提取数据的特定部分。
    *   `filter`：按物种存在与否或序列长度过滤 block。
    *   `slice`：使用 runlist 提取特定的比对切片。
    *   `subset`：从 block 中提取物种子集。
*   **转换（Transform）**：修改或合并 block FA 文件。
    *   `concat`：连接同一物种的序列片段。
    *   `consensus`：使用 POA（偏序比对）生成一致性序列。
    *   `join`：基于共同的目标序列合并多个文件。
    *   `multiz`：使用类 multiz 的带状动态规划算法合并 block FA 文件。
    *   `refine`：使用内置或外部工具对 block 内的序列进行重新比对。
    *   `replace`：使用映射文件替换序列头。
*   **文件（File）**：创建或拆分 block FA 文件。
    *   `create`：根据区间链接创建 block FA 文件。
    *   `separate`：按物种将 block 拆分为独立文件。
    *   `split`：按比对块或染色体拆分 block FA 文件。
*   **变异（Variation）**：从比对中 calling 变异。
    *   `to-vcf`：将替换（SNP）导出为 VCF 格式。
    *   `to-xlsx`：将替换和 indel 导出为 Excel 文件。
    *   `variation`：以 TSV 格式列出变异（替换）。

---

## 信息命令

### check

检查 block 头中指定的基因组位置是否对参考基因组 FA 文件有效。

```bash
pgr fas check [OPTIONS] --genome <genome> <infiles>...
```

*   `-g, --genome <path>`：参考基因组 FA 文件路径（必填）。
*   `-n, --name <name>`：仅检查特定物种的序列。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

输出格式（制表符分隔）：每行包含条目区间及其状态（`OK` 或 `FAILED`）。

### cover

输出比对在染色体上的覆盖区域，格式为 JSON。

```bash
pgr fas cover [OPTIONS] <infiles>...
```

*   `-n, --name <name>`：仅输出该物种的覆盖区域。
*   `--trim <int>`：将比对边界向内修剪 N 个碱基以避免重叠（对 lastz 结果有用）。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

### link

输出比对 block 中区间（基因组坐标）之间的链接。

```bash
pgr fas link [OPTIONS] <infiles>...
```

*   `--pair`：输出双边（成对）链接。
*   `--best`：基于序列距离输出最近邻双边链接（去重）。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

注意：

*   `--pair` 与 `--best` 互斥。

输出格式：每行以制表符分隔。默认情况下，一个 block 中的所有区间输出在同一行。使用 `--pair` 或 `--best` 时，每行包含两个区间。

### name

提取 block FA 文件中的所有物种名。

```bash
pgr fas name [OPTIONS] <infiles>...
```

*   `-C, --count`：同时输出每个物种名的出现次数。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

### stat

计算每个比对 block 的基本统计量（长度、可比较碱基数、差异、gap 等）。

```bash
pgr fas stat [OPTIONS] <infiles>...
```

*   `--outgroup`：将每个 block 的最后一条序列视为外群。`length` 列始终反映完整比对长度（含外群）；其他统计量（`comparable`、`difference`、`gap`、`ambiguous`、`D`、`indel`）均基于排除外群后的序列计算。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

输出列（制表符分隔）：

*   `target`：block 的目标区间。
*   `length`：包含 gap 的比对长度。
*   `comparable`：所有序列均为非模糊碱基的位置数。
*   `difference`：可比较碱基中的多态位置数。
*   `gap`：所有序列均含 gap 的位置数。
*   `ambiguous`：至少含一个模糊碱基且不含 gap 的位置数。
*   `D`：所有序列对之间的平均成对分歧度。
*   `indel`：所有 indel 区域的总跨度。

---

## 子集命令

### filter

根据物种存在与否和序列长度过滤 block，并可选择性地格式化序列。

```bash
pgr fas filter [OPTIONS] <infiles>...
```

*   `-n, --name <name>`：用于长度过滤的物种。不包含该物种的 block 会被跳过。默认使用每个 block 的第一个物种。
*   `--min-len <int>`：保留所选物种比对长度（含 gap）大于等于该值的 block。
*   `--max-len <int>`：保留所选物种比对长度（含 gap）小于等于该值的 block。
*   `--upper`：将序列转换为大写。
*   `--dash`：从序列中移除 dash（gap）。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

### slice

根据提供的 runlist（JSON）提取特定的比对切片。

```bash
pgr fas slice [OPTIONS] --runlist <runlist.json> <infiles>...
```

*   `--runlist <file>`：描述要提取区间的 JSON 文件（必填）。
*   `-n, --name <name>`：参考物种名（默认：第一个 block 的第一个物种）。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

### subset

从比对 block 中提取物种子集。

```bash
pgr fas subset [OPTIONS] --required <name.lst> <infiles>...
```

*   `-R, --required <file>`：包含要保留物种名的文件，每行一个（必填）。
*   `--strict`：跳过不包含所有必需物种的 block。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

---

## 转换命令

### concat

连接多个 block 中同一物种的序列片段。

```bash
pgr fas concat [OPTIONS] --required <name.lst> <infiles>...
```

*   `-R, --required <file>`：包含要保留/排序物种名的文件，每行一个（必填）。
*   `--phylip`：以 relaxed PHYLIP 格式输出，而非 FASTA。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

注意：

*   缺失的序列用 gap（`-`）填充。

### consensus

使用偏序比对（POA）为每个 block 生成一致性序列。

```bash
pgr fas consensus [OPTIONS] <infiles>...
```

*   `--engine <builtin|spoa>`：使用的 POA 引擎（默认：builtin）。
*   `--match <int>`：匹配碱基的得分（默认：5）。
*   `--mismatch <int>`：不匹配碱基的得分（默认：-4）。
*   `--gap-open <int>`：gap 开放罚分（默认：-8）。
*   `--gap-extend <int>`：gap 延伸罚分（默认：-6）。
*   `--align-mode <local|global|semi_global>`：比对模式（默认：global）。
*   `--consensus-name <name>`：一致性序列的名称（默认：consensus）。
*   `--outgroup`：表示每个 block 的最后一条序列为外群。外群不参与一致性计算，但会保留在输出 block 中。
*   `-p, --parallel <int>`：线程数（默认：1）。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

### join

基于共同的目标序列合并多个 block FA 文件。

```bash
pgr fas join [OPTIONS] <infiles>...
```

*   `-n, --name <name>`：目标物种名。默认使用第一个 block 的第一个物种，并作为所有 block 的共同目标。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

### multiz

在共享的参考坐标系下，使用类 multiz 的带状动态规划算法合并多个 block FA 文件。

```bash
pgr fas multiz [OPTIONS] --ref-name <name> <infiles>...
```

*   `-r, --ref-name <name>`：所有输入中都存在的参考序列名（必填）。
*   `--radius <int>`：参考对角线周围的带状 DP 半径（默认：30）。
*   `--min-width <int>`：参与合并的最小窗口宽度（默认：1）。
*   `--mode <core|union>`：合并模式（默认：core）。
*   `--score-scheme <file>`：评分方案文件（LASTZ 格式）或预设名（如 `hoxd55`）。
*   `--gap-model <constant|medium|loose>`：gap 模型（默认：medium）。
*   `--align-gap-open <int>`：比对 gap 开放罚分，覆盖 `--gap-model` 的默认值。
*   `--align-gap-extend <int>`：比对 gap 延伸罚分，覆盖 `--gap-model` 的默认值。
*   `--match-score <int>`：匹配得分（默认：2）。
*   `--mismatch-score <int>`：不匹配罚分（默认：-1）。
*   `--gap-score <int>`：gap 罚分（默认：-2）。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

### refine

使用内置或外部 MSA 工具对 block 内的序列进行重新比对。

```bash
pgr fas refine [OPTIONS] <infiles>...
```

*   `--engine <program>`：比对程序：`builtin`（默认）、`clustalw`、`mafft`、`muscle`、`spoa`、`none`。
*   `--outgroup`：表示存在外群。
*   `--chop <usize>`：修剪头部和尾部的 indel（默认：0，即禁用）。
*   `--quick`：快速模式，仅比对 indel 邻近区域。
*   `--indel-pad <int>`：快速模式下，扩大 indel 区域（默认：50）。
*   `--fill <int>`：快速模式下，填充 indel 之间的空洞（默认：50）。
*   `-p, --parallel <int>`：线程数（默认：1）。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

### replace

使用映射文件替换 block FA 文件中的序列头。

```bash
pgr fas replace [OPTIONS] --replace-tsv <replace.tsv> <infiles>...
```

*   `--replace-tsv <file>`：包含替换规则的 TSV 文件（必填）。每行是一个制表符分隔的列表：
    *   一个字段：如果该名唯一匹配 block 中的一个头，则丢弃整个 block。
    *   两个字段：`original_name<TAB>new_name`，替换匹配的头。
    *   三个或更多字段：对第一个字段之后的每个替换名复制一次 block。
    *   如果一个 block 包含多个匹配头，则保持 block 不变并发出警告。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

同一个 block 内出现多次的头也视为多个匹配头，该 block 将保持不变。

---

## 文件命令

### create

根据区间链接（例如来自 `pgr fas link`）创建 block FA 文件。

```bash
pgr fas create [OPTIONS] --genome <genome> <infiles>...
```

*   `-g, --genome <file>`：参考基因组 FA 文件路径（必填）。
*   `-n, --name <name>`：为区间设置物种名（默认：从区间字符串中的物种名推断）。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

注意：

*   参考基因组 FA 文件可以是纯文本或 bgzip 压缩。

### separate

将 block FA 文件按物种拆分为独立文件。

```bash
pgr fas separate [OPTIONS] <infiles>...
```

*   `-s, --suffix <string>`：输出文件扩展名（默认：.fasta）。
*   `--rc`：如果链为 '-'，则对序列进行反向互补。
*   `-o, --outdir <dir>`：输出目录（默认：stdout）。

注意：

*   序列中的 dash 会被移除。
*   已存在的输出文件会被覆盖。

### split

将 block FA 文件按比对块或染色体拆分为文件。

```bash
pgr fas split [OPTIONS] <infiles>...
```

*   `--chr`：按染色体拆分文件。
*   `--simple`：简化头信息，仅保留物种名。同时作用于 stdout 和按文件输出。
*   `-s, --suffix <string>`：输出文件扩展名（默认：.fas）。
*   `-o, --outdir <dir>`：输出目录（默认：stdout）。

---

## 变异命令

### to-vcf

将替换（SNP）导出为 VCF 格式。

```bash
pgr fas to-vcf [OPTIONS] <infiles>...
```

*   `--sizes <file>`：染色体长度文件，用于输出 `##contig` 头。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

注意：

*   所有 block 必须包含相同物种且顺序一致，因为 VCF 使用固定的样本头。

### to-xlsx

将变异（替换和 indel）导出为带格式的 Excel 文件。

```bash
pgr fas to-xlsx [OPTIONS] <infiles>...
```

*   `--indel`：包含 indel。
*   `--outgroup`：表示存在外群。
*   `--no-single`：省略单例变异。
*   `--no-complex`：省略复杂变异。
*   `--min-freq <float>`：最小频率。
*   `--max-freq <float>`：最大频率。
*   `--wrap <int>`：可视化换行长度（默认：50）。
*   `-o, --outfile <file>`：输出文件名（默认：variations.xlsx）。

### variation

以 TSV 格式列出变异（替换）。

```bash
pgr fas variation [OPTIONS] <infiles>...
```

*   `--outgroup`：表示存在外群。
*   `-o, --outfile <file>`：输出文件名（默认：stdout）。

注意：

*   `--outgroup` 要求每个 block 至少包含 2 条序列。

输出列（制表符分隔）：

*   `#target`：block 的目标区间。
*   `chr`：染色体名。
*   `chr_pos`：染色体上的位置。
*   `range`：格式为 `chr:pos` 的染色体位置。
*   `pos`：比对内的位置（1-based）。
*   `tbase`、`qbase`、`bases`、`mutant_to`、`freq`、`pattern`、`obase`：来自替换记录的字段（见 `pgr::libs::alignment::Substitution`）。
