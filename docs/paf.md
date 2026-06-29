# pgr paf

`pgr paf` 用于处理**PAF** (Pairwise mApping Format) 文件。它把成对比对当作一张隐式泛基因组图，
将目标区间通过比对网络投影出去，回答"哪些序列的哪些区域与该位点同源"这一问题。查询类子命令 （`query`
/ `to-bed` / `to-maf` / `to-vcf` / `to-gfa`）按需遍历这张隐式图，不物化整张图；而 `graph` 和 `stat`
子命令会物化一张**粗粒度**的全基因组图，只在 ≥ `--min-var-len`（默认 100 bp）的结构变异处切分节点，
小 indel 留在节点内部。

## 为什么在泛基因组场景使用 PAF

泛基因组分析的起点是**全两两比对 (all-versus-all)**：每一对基因组（或单倍型）都要互相比对。对于*N*
个基因组，比对记录数为*O(N²)*，单个位点可能出现在其中*N−1*条记录里。

绝大多数传统比对格式都把**序列字符串**嵌入每条记录：

- **SAM/BAM**每条记录都有 `SEQ` 字段。
- **MAF**在每个 block 中把每条序列的碱基（含 gap）全部打印出来。
- **AXT**把 query 和 target 的序列原样保存。

在 all-vs-all 场景下这是灾难性的。同一段 10 kb 区域可能被存储几十次（每条触及它的成对比对都存一份），
因此比对文件常常**比输入 FASTA 还大一两个数量级**。

**PAF 采取了完全不同的策略**：它只存**坐标和匹配统计**（12 列加可选 tag 如 `cg:Z:` CIGAR），
**完全不存序列**。这就把比对索引与序列仓库解耦了：

- 序列只存一份在 BGZF 压缩的 FASTA 中，需要时按坐标随机 seek。
- 比对只存坐标，不论该位点被复用多少次，PAF 中都只出现一次。

结果是 all-vs-all 的 PAF 通常**比输入 FASTA 还小**，即使*N*增长到几百、几千仍然可控。
序列物化被推迟到查询时：`pgr paf` 只取回请求区域真正需要的那几个切片，从不加载整张比对。

PAF 适合这一角色的其他原因：

- **多个主流比对器都能直接输出 PAF**，无需格式转换即可作为 `pgr paf` 的输入：
    - **minimap2**（`minimap2 -c --eqx`）：泛基因组与读长比对事实上的标准工具，输出带 `cg:Z:` CIGAR
      的 PAF。
    - **wfmash**：seqwish / PGGB / impg 流程的默认 aligner，由 MashMap 找同源块、WFA 做碱基级比对，
      常用于 all-vs-all 场景。
    - **FastGA**：impg 支持的备选 aligner，磁盘占用更低，适合内存受限场景。
    - **MashMap / MashMap2 / MashMap3**：基于 k-mer 距离定位加间隙式比对，输出 PAF-like 格式。
    - **paftools.js**（minimap2 自带）：用于后处理和过滤 minimap2 的 PAF。注：`pgr paf`
      本身不绑定任何 aligner，只要 PAF 行带 `cg:Z:` CIGAR 且 `=`/`X` 操作符齐全即可
      （`minimap2 --eqx` 与 `wfmash` 均满足此要求）。
- **Tab 分隔、流式友好**：无 header，不需要 seek，易于管道化、切分和并行。
- **`cg:Z:` CIGAR 足以按需重建任意成对比对**（`pgr paf to-maf` 默认模式即采用此方式），
  因此不存序列也不会丢失信息。
- **坐标-only 记录天然可传递组合**：`pgr paf query -t` 能链接 A→B 与 B→C，而无需物化 A↔C 的成对比对，
  这是泛基因组遍历可扩展的关键。

## 功能概览

给定一组成对比对（带 `cg:Z:` CIGAR 的 PAF）和一组 BGZF 压缩的 FASTA，`pgr paf` 能够：

- **Index**：把比对构建成可复用的区间树索引（`.paf.idx`）。
- **Query**：查询一个目标区域，把坐标投影到所有比对上的 query，可选传递式 BFS 走多跳同源链。
- **Export**：把收集到的同源片段导出为 BED、pairwise / multi-way MAF、局部 GFA 或 multi-way VCF。
- **Build**：构建粗粒度的全局泛基因组图（seqwish 风格 DSU），或输出其拓扑统计报告。

## 工作原理

`pgr paf` 使用 per-target 的区间树做快速范围查找，CIGAR 字符串以紧凑 delta 形式存储。
对于传递式查询，它构建一个双向 mirror 索引：一条 `A→B` 的 PAF 记录同时也启用 `B→A` 的遍历，
然后走 BFS 直到 `--max-depth` 跳。对于 multi-way 输出（`to-maf --msa` / `to-gfa` / `to-vcf`），
它把一个区域的所有同源片段喂给 Partial Order Alignment (POA) 引擎产出多序列比对，再由 MSA
推导出对应的图或变异。

## Quick start

下面的示例是自包含的，逐块复制粘贴即可构建一个 3 基因组的小 demo（1 个 SNP + 1 个 2bp 插入），
并端到端跑通每一个 `pgr paf` 子命令。

### 1. 准备 demo 数据

```bash
mkdir -p pgr-paf-demo && cd pgr-paf-demo

# 三个基因组：A (参考), B (pos 12 处 1 个 SNP), C (pos 10 处 2bp 插入)
printf '>A\nACGTACGTACGTACGTACGTACGTACGTAC\n>B\nACGTACGTACGTGCGTACGTACGTACGTAC\n>C\nACGTACGTACTTGTACGTACGTACGTACGTAC\n' > refs.fa

# BGZF 压缩 FASTA（-f TSV 要求，用于随机访问）
pgr fa gz refs.fa          # 产出 refs.fa.gz + refs.fa.gz.gzi

# FASTA TSV: genome_name<TAB>bgzf_path (这里 3 个基因组共用一个文件)
printf 'A\trefs.fa.gz\nB\trefs.fa.gz\nC\trefs.fa.gz\n' > genomes.tsv

# PAF 比对 (query<TAB>target=A): B 有 1 SNP, C 有 2bp 插入
printf 'B\t30\t0\t30\t+\tA\t30\t0\t30\t29\t30\t255\tcg:Z:12=1X17=\nC\t32\t0\t32\t+\tA\t30\t0\t30\t30\t30\t255\tcg:Z:10=2I20=\n' > aln.paf
```

### 2. 建索引

```bash
pgr paf index aln.paf -o aln.paf.idx
```

### 3. 查询区域

```bash
# 单跳：找出所有与 A:0-30 比对上的序列
pgr paf query aln.paf.idx A:0-30

# 传递式 BFS：也会走 mirror 索引的 A→B、A→C
pgr paf query aln.paf.idx A:0-30 -t
```

### 4. 导出同源片段

```bash
# Pairwise MAF（每条比对一个 block，由 CIGAR 重建）
pgr paf to-maf aln.paf.idx A:0-30 -f genomes.tsv

# Multi-way MSA MAF（POA 把所有序列合并成一个 block）
pgr paf to-maf aln.paf.idx A:0-30 -t --msa -f genomes.tsv

# 局部 GFA 图（SNP bubble + indel bubble）
pgr paf to-gfa aln.paf.idx A:0-30 -t -f genomes.tsv

# Multi-way VCF（pos 10 处 1 个 INS，pos 13 处 1 个 SNP）
pgr paf to-vcf aln.paf.idx A:0-30 -t -f genomes.tsv

# 仅 BED3 坐标
pgr paf to-bed aln.paf.idx A:0-30 -t
```

### 5. 构建粗粒度全局图与报告

```bash
# 粗粒度泛基因组 GFA（默认 < 100bp 的 indel 留在节点内部）
pgr paf graph aln.paf -f refs.fa -o graph.gfa

# 拓扑报告（约 25 个维度，TSV）
pgr paf stat aln.paf -f refs.fa
```

## 子命令

### `index` — 构建可复用索引

```bash
pgr paf index <infiles>... [-o <idx>]
```

从一个或多个 PAF 文件构建 per-target 的区间树索引。多个文件会被合并成统一索引；
同名序列在不同文件中共享同一个内部 ID。使用 `-o` 可将索引持久化到 `.paf.idx`，后续查询启动即可用。

- 支持纯文本和 gzip（`.gz` / BGZF）PAF。
- BGZF 输入启用 CIGAR 懒加载（CIGAR 按需通过虚拟文件偏移读取，降低内存）。

```bash
pgr paf index a.paf b.paf -o merged.paf.idx
```

### `query` — 通过比对投影区域

```bash
pgr paf query <infile> <region> [options]
pgr paf query <infile> -b <bed> [options]
```

两种模式：

- **默认**（单跳）：找出所有 target 区间与查询区域重叠的 PAF 记录，把坐标 lift 到 query 序列上。
- **`-t/--transitive`**（多跳 BFS）：通过中间序列迭代投影，直到 `--max-depth` 跳。双向 mirror
  索引允许双向遍历。

输出为 PAF（12 列加 `gi`/`bi`/`cg` tags）。要输出 BED / MAF / GFA / VCF，请分别使用 `to-bed` /
`to-maf` / `to-gfa` / `to-vcf`。

```bash
# 从 PAF 文件单跳查询（索引即时构建）
pgr paf query aln.paf A:0-30

# 带一致性过滤的传递式 BFS
pgr paf query aln.paf A:0-30 -t --min-identity 0.8

# 从 BED 文件批量查询
pgr paf query aln.paf.idx -b regions.bed
```

### `to-bed` — 输出 BED3 坐标

```bash
pgr paf to-bed <infile> <region> [options]
```

`pgr paf query` 的管道友好、仅坐标视图。每条查询结果输出 `name<TAB>start<TAB>end`。所有 query 选项
（区域、`--transitive`、过滤）都支持。

```bash
pgr paf to-bed aln.paf A:0-30 -t
```

### `to-maf` — 输出 pairwise 或 multi-way MAF

```bash
pgr paf to-maf <infile> <region> -f <tsv> [options]
```

- **默认**（pairwise）：每条查询结果变成一个 2 序列 MAF block，直接由 CIGAR 重建。假定比对已被
  chain/net 精炼过，不再做 POA 再精炼。
- **`--msa`**（multi-way）：把每个区域的所有查询结果合并成单个多序列 MAF block，由 POA 完成。CIGAR
  被忽略；序列（target 在前，然后每个 query，`-` 链反向互补）喂给 POA 引擎。建议配合 `--transitive`
  使用。

`-f/--fasta-tsv`（必填）为两列 TSV，即 `genome_name<TAB>bgzf_fasta_path`。每个 genome name 必须与
PAF 索引中的 query/target 名字匹配。要求 BGZF 格式（由 `pgr fa gz` 产出），以便随机访问。

```bash
# Pairwise MAF
pgr paf to-maf aln.paf A:0-30 -f genomes.tsv

# 带传递式 BFS 的 multi-way MSA
pgr paf to-maf aln.paf A:0-30 -t --msa -f genomes.tsv
```

### `to-vcf` — 输出 multi-way VCF

```bash
pgr paf to-vcf <infile> <region> -f <tsv> [options]
```

对每个区域的所有同源片段运行 POA 多序列比对（无需 `--msa` 标志，`to-vcf` 始终走 POA 路径），再从
MSA 中调用替换和 indel，产出 VCF。三类变异：

- **SNP**：单个 target 非 gap 列上有 ≥1 个 query 不同。
- **INS**：连续的 target gap 列（锚点为 gap 之前的 target 碱基）。
- **DEL**：连续的 target 非 gap 列上有 ≥1 个 query 是 gap。

Indel 会相对参考做左对齐：当锚点之前的参考碱基与每个非空 indel 序列的最后一个碱基相等时，锚点左移。
GT 字段编码每个样本的等位基因（0=REF, 1..=N=ALT 索引, `.`=gap 或非 ACGT）。REF 为 target 序列。

```bash
pgr paf to-vcf aln.paf A:0-30 -t -f genomes.tsv -o out.vcf
```

### `to-gfa` — 通过 POA 图输出局部 GFA

```bash
pgr paf to-gfa <infile> <region> -f <tsv> [options]
```

从每个区域的 POA 多序列比对中导出局部 GFA (v1.0) 图。POA 图的节点是碱基、边是邻接关系、
路径追踪每条输入序列；经过压缩（单碱基节点的线性段合并为多碱基 segment）后导出为 GFA S/L/P 行。
SNP/indel bubble 保留为图分支。

- **`--crush`**：可选的 impg `crush` 风格后处理，压缩 SNP bubble（共享相同入/出邻居的节点）
  为单个节点，保留权重最高的等位基因。会丢失碱基级 ALT 信息，仅用于 SV 概览图。
- 每个区域产出独立的 GFA block（节点 ID 从 1 重新计数）。多个区域之间用 `# region: <name>`
  注释分隔。

```bash
# 带传递式 BFS 的局部图
pgr paf to-gfa aln.paf A:0-30 -t -f genomes.tsv

# bubble 压缩后的 SV 概览图
pgr paf to-gfa aln.paf A:0-30 -t -f genomes.tsv --crush
```

### `graph` — 构建粗粒度全局 GFA

```bash
pgr paf graph <infile> -f <fasta>... [--min-var-len <n>] [-o <gfa>]
```

采用 seqwish 风格的 segment 级 DSU 算法，从所有 PAF 比对和一组 FASTA 序列构建粗粒度泛基因组图（GFA
v1.0）：

1. 在 ≥ `--min-var-len` 的 indel 处把每条比对切分成 match segment。
2. 通过并查集把对齐的 segment 合并（传递闭包）。
3. 推导出图节点（DSU 类）、边（路径邻接）和新颖 segment（未对齐 gap），输出 S/L/P 行。

只有大的结构变异（≥ `--min-var-len`，默认 100）会切分节点；小 indel 留在节点内部。S 行带 rGFA tag
（`SN:Z` 源序列, `SO:i` 0-based 起点, `SR:i:0`）。

```bash
# 默认 SV 阈值（100bp）
pgr paf graph aln.paf -f refs.fa -o graph.gfa

# 更严格阈值（只有 ≥500bp 的 SV 才切分节点）
pgr paf graph aln.paf -f refs.fa --min-var-len 500 -o graph.gfa
```

### `stat` — 报告图拓扑指标

```bash
pgr paf stat <infile> -f <fasta>... [--min-var-len <n>] [-o <tsv>]
```

在 PAF 比对诱导出的粗粒度泛基因组图上计算拓扑报告（TSV: `key<TAB>value`），构建路径与
`pgr paf graph` 相同。约 25 个维度：

- **基础**：segments, links, paths, path_steps, total_segment_bp
- **节点长度分布**：min/mean/median/max
- **节点覆盖**：mean/median, singleton_nodes, reused_nodes, reused_nodes_cross_path
- **结构**：components, largest_component_nodes, tips, isolated_nodes, self_loop_edges
- **路径长度分布**：min/median/max（按步数和按 bp）

可在不物化 GFA 的情况下评估图质量，或对比不同 `--min-var-len` 阈值下构建的图。

```bash
pgr paf stat aln.paf -f refs.fa -o report.tsv

# 对比不同阈值
pgr paf stat aln.paf -f refs.fa --min-var-len 50  -o r50.tsv
pgr paf stat aln.paf -f refs.fa --min-var-len 500 -o r500.tsv
```

## 通用 query 选项

`query`、`to-bed`、`to-maf`、`to-vcf`、`to-gfa` 子命令共享以下选项（由 `add_query_args` 定义）：

| 选项                            | 默认值 | 说明                                           |
|---------------------------------|--------|------------------------------------------------|
| `<region>`                      | —      | 目标区域（如 `chr1:1000-5000`）                |
| `-b, --bed-regions <file>`      | —      | 含多个区域的 BED 文件                          |
| `-t, --transitive`              | off    | 启用传递式 BFS 遍历                            |
| `-m, --max-depth <n>`           | 2      | BFS 最大深度（0 = 无限）                       |
| `--min-len <n>`                 | 10     | 传递的最小区间长度                             |
| `--min-dist <n>`                | 10     | 合并相邻区间的最小距离                         |
| `--min-identity <f>`            | 0.0    | 最小 gap-compressed 一致性 (0.0-1.0)           |
| `--min-output-len <n>`          | 0      | 输出区间最小长度（0 = 不过滤）                 |
| `--merge-distance <n>`          | 0      | 在此距离内合并相邻输出区间（0 = 关闭）         |
| `--min-degree <n>`              | 0      | 每个区域最少 distinct query 序列数（0 = 关闭） |
| `--min-chain-length <n>`        | 0      | 每个 query 最小总比对长度（0 = 关闭）          |
| `--subset-sequence-list <file>` | —      | 要包含的序列名列表文件                         |
| `--syntenic-filter <chain>`     | —      | UCSC chain 文件，用于丢弃非共线性的 query 结果 |

## FASTA TSV 格式

`-f/--fasta-tsv` 参数（`to-maf`、`to-vcf`、`to-gfa` 必填）是两列 tab 分隔文件：

```
genome_name<TAB>bgzf_fasta_path
```

- 每个 `genome_name` 必须与 PAF 索引中的 query/target 名字匹配。
- 一个 FASTA 文件可以被多个 genome name 引用（多染色体）。
- PAF 索引中的所有 genome name 都必须出现在 TSV 中（严格校验，缺失项会硬报错）。
- FASTA 文件必须是 BGZF 压缩（由 `pgr fa gz` 产出），以通过 `.gzi` 索引实现随机访问。

`genomes.tsv` 示例：

```
A    /data/genomes/A.fa.gz
B    /data/genomes/B.fa.gz
C    /data/genomes/C.fa.gz
```

## 注意事项

- 输入 PAF 文件应包含 `cg:Z:` tag，以便准确做坐标投影和图切分。缺少 CIGAR 时部分子命令会降级或报错。
- `index` / `query` / `to-bed` 支持 PAF 与 `.paf.idx` 输入；`to-maf` / `to-vcf` / `to-gfa`
  同样支持两种输入，但都额外要求 `-f/--fasta-tsv` 提供 BGZF FASTA。
- `graph` / `stat` 的 `-f` 参数直接接一个或多个 FASTA 文件（不是 TSV），用于提供节点序列和长度。
- 支持纯文本和 gzip（`.gz`）文件（含 BGZF）。BGZF 输入启用 CIGAR 懒加载，降低内存占用。
- 输入文件为 `stdin` 时从 stdin 读取 PAF。

## 另见

- [`pgr maf to-paf`](maf.md)：把 MAF 转成 PAF 作为 `pgr paf` 的输入。
- [`pgr fa gz`](fa.md)：BGZF 压缩 FASTA（`-f` TSV 要求）。
- [`pgr fa range`](fa.md)：按坐标从 BGZF FASTA 提取子序列。

