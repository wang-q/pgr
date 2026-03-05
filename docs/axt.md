# axt

`pgr axt` 模块用于处理 **AXT** 格式的成对基因组比对（Pairwise Genomic Alignments）。该格式通常由 `lastz` 或 `blastz` 等比对工具生成。

## 核心定位

- **定位**：AXT 格式的操作与转换工具。
- **输入**：AXT 格式文件（支持 gzip）。
- **输出**：排序后的 AXT，或转换为其它标准格式（FASTA, MAF, PSL）。
- **互补**：
  - 上游：`lastz`/`blastz` (生成 AXT)。
  - 下游：`pgr maf`, `pgr psl` (进一步处理转换后的格式)。

## 子命令详解

### 1. `pgr axt sort`: 排序与重编号
对 AXT 文件进行排序，支持多种排序键值。

- **排序模式**:
  - **默认**: 按 Target (Reference) 的名称和起始位置排序。
  - `--query`: 按 Query 的名称和起始位置排序。
  - `--by-score`: 按比对得分（Score）降序排列。
- **重编号 (`--renumber`/`-r`)**:
  - 排序后重新分配 ID（从 0 开始），确保 ID 的唯一性和顺序性。
  - 类似 UCSC `axtSort` 的行为。

### 2. `pgr axt to-fas`: 转换为 Block FASTA
将 AXT 记录转换为 Block FASTA (每条比对记录对应一对序列块)。

- **用途**: 用于需要逐块序列分析的场景。
- **坐标处理**:
  - 需要提供 Query 基因组的染色体大小文件 (`chr.sizes`)，以正确处理负链坐标。
  - 输出格式为 `>Target:Start-End` 和 `>Query:Start-End`。
- **参数**:
  - `chr.sizes`: 必需。Query 基因组的大小文件。
  - `--tname`/`--qname`: 自定义输出 FASTA 头中的物种名称前缀。

### 3. `pgr axt to-maf`: 转换为 MAF 格式
将 AXT 转换为 Multiple Alignment Format (MAF)。

- **用途**: MAF 是基因组多重比对的标准格式，支持更丰富的元数据。
- **必需参数**:
  - `-t`/`--t-sizes`: Target 基因组大小文件。
  - `-q`/`--q-sizes`: Query 基因组大小文件。
- **高级选项**:
  - `--t-split`: 按 Target 序列名称将输出拆分到不同文件（输出为目录）。
  - `--t-prefix`/`--q-prefix`: 为序列名称添加前缀（如 `hg38.chr1`）。

### 4. `pgr axt to-psl`: 转换为 PSL 格式
将 AXT 转换为 PSL (UCSC BLAT) 格式。

- **用途**: PSL 格式常用于基因组浏览器的可视化。
- **坐标转换**: 自动处理 AXT (负链相对坐标) 到 PSL (正链坐标 + Strand 标记) 的转换。
- **必需参数**:
  - `-t`/`--t-sizes`: Target 基因组大小文件。
  - `-q`/`--q-sizes`: Query 基因组大小文件。

## 典型用法

### 场景 A：排序并转换为 MAF

```bash
# 1. 排序 AXT (按 Target)
pgr axt sort raw.axt -o sorted.axt

# 2. 转换为 MAF (添加物种前缀)
pgr axt to-maf sorted.axt \
    -t target.sizes -q query.sizes \
    --t-prefix "hg38." --q-prefix "mm10." \
    -o out.maf
```

### 场景 B：转换为 Block FASTA 用于分析

```bash
# 转换，指定物种名称
pgr axt to-fas query.sizes raw.axt \
    --tname Human --qname Mouse \
    -o blocks.fas
```

### 场景 C：按 Target 拆分 MAF

```bash
# 输出到目录 split_mafs/，每个 Target 序列一个文件
pgr axt to-maf sorted.axt \
    -t target.sizes -q query.sizes \
    --t-split -o split_mafs/
```
