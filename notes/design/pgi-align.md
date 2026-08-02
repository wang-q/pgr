# pgi 两索引归并比对（设计定稿：v1 最小闭环）

> 定位：`.pgi` 的第一个比对消费者。输入两个已构建的 .pgi，
> 输出 PSL 块，喂给 `pgr pl chainnet`（UCSC 链化由 pgr 承担，见
> [[fastga.md]] §12.3 决策 3）。
> 状态：2026-08-02 定稿，随实现迭代。

## 1. 范围（v1 做什么、不做什么）

**做**：

1. 两个排序 .pgi 流的线性归并 → 种子命中（频率过滤）；
2. anti-diagonal 空间贪心链化 → 链（tube 语义简化）；
3. 每条链输出一个 PSL 块（链覆盖区间，无 CIGAR 细节）；
4. `pgr pgi align` CLI + 集成测试 + E. coli 验证。

**不做（v2 及以后）**：

- 局部扩展（banded/POA DP 细化 CIGAR，需要序列输入）——v1 链块直接由
  `pgr psl to_chain` / chainnet 接手（块结构即输入，UCSC 链化自算 score）；
- lcp 连续传播（adaptamer 变长种子）——固定 k 种子先行；
- pbit 内嵌索引段消费（待 v1002 落地后接线）；
- mmap/流式读取（E. coli 规模直接整体载入内存）。

## 2. 数据流

```
ref.pgi + query.pgi（v1 无需序列文件：链块只依赖索引的 k-mer/位置/contig 表）
  │ 1. 排序流归并（频率过滤）→ 种子命中
  ▼
hit = (key, a_contig, a_pos, a_strand, b_contig, b_pos, b_strand)
  │ 2. 方向解析 + diag/anti 坐标变换
  ▼
fwd/rev 两个空间：(pos_a, pos_b'), diag = pos_a − pos_b', anti = pos_a + pos_b'
  │ 3. 按 (contig_a, contig_b, 方向) 贪心链化
  ▼
chain = 区间 + 对角线带（间距/带宽容忍，跨度过滤）
  │ 4. PSL 块输出（q=query, t=ref）
  ▼
out.psl → pgr psl to_chain → pgr pl chainnet（现有字节级验证的链化主场）
```

## 3. 关键语义

### 3.1 种子方向（pgi 双链条目）

`.pgi` 每条位置带 strand 标记（0=正向 k-mer，1=RC k-mer）。两索引归并时
key 相等按 (a_strand, b_strand) 解析方向：

| (a, b) strand | 含义 | 变换 |
|---|---|---|
| (0,0) / (1,1) | 正链命中（两侧实际窗口相等） | `pos_b' = pos_b` |
| (0,1) / (1,0) | 负链命中（b 窗口 = RC(a 窗口)） | `pos_b' = b_len − k − pos_b` |

负链 `pos_b'` 是 RC(b) 空间坐标；PSL 输出时用
`reverse_range(q_start, q_end, q_size)` 还原原始坐标、strand 记 `-`。

### 3.2 链化（简化 tube）

按 (contig_a, contig_b, 方向) 分组，命中按 (diag, pos_a) 排序后贪心：

- 当前链从首个命中开始；
- 延伸条件：`|diag − 链均线| ≤ band` 且 `Δpos_a ≤ max_gap` 且 `Δpos_b ≤ max_gap`；
- 否则收链：两侧轴向上种子跨度 `(last + k − first) ≥ min_span` 才保留。

参数默认对齐 FastGA（-f 10 / -c 85 / -s 1000，本实现为未加倍空间）：

| 参数 | 默认 | FastGA 对应 |
|---|---|---|
| `--freq` | 10 | `-f`（任一侧频率超限即跳过） |
| `--min-span` | 85 | `-c`（CHAIN_MIN） |
| `--max-gap` | 1000 | `-s`（CHAIN_BREAK/2） |
| `--band` | 128 | 对角线带宽（tube dgmin..dgmax） |

### 3.3 PSL 输出

- q = query 基因组，t = ref 基因组（`pgr pgi align <ref> <query>`）；
- 每条链 = 一个块：q 区间 `[min_pos_b'..max_pos_b'+k)`，t 区间 `[min_a..max_a+k)`；
- 负链块 q_start/q_end 走 RC 空间 → `reverse_range` 还原；
- match/mismatch 计数置 0（链化只消费块结构；真实计数待 v2 扩展）。

## 4. CLI

```
pgr pgi align <ref.pgi> <query.pgi> -o out.psl [--freq 10] [--min-span 85] [--max-gap 1000] [--band 128]
```

- 参数默认对齐 FastGA（见 §3.2）；两侧索引参数必须一致（复用 `dist pgi` 校验）。
- v2 局部扩展阶段再引入序列文件（2bit 随机访问 / FASTA 整载）。

## 5. 验证

1. 单元测试：归并方向解析（fwd/rev）、频率过滤、链化（band/gap/span 边界）；
2. 集成测试：小序列人工 PSL 形状断言；
3. E. coli：MG1655 自比对（`pgr pgi align`）全基因组覆盖 ~100%；
   MG1655 vs 近缘株与 `FastGA -psl` 结果做覆盖/共线性结构对照（非字节一致）。

## 5.1 实测结果（2026-08-02，release，默认参数）

**MG1655 自比对**（`mg1655.pgi` vs 自身）：

- 745 块；其中**正链主链 = 1 块覆盖 4,641,650 / 4,641,652 bp（99.9999%）**；
- 其余小块为真实重复结构（rRNA 操纵子、IS 元件等），负链 186 块对应反向重复
  （rRNA 操纵子的反向拷贝）。

**MG1655 vs Sakai（O157）**：

| 指标 | pgr pgi align | FastGA -psl |
|---|---|---|
| 块数 | 1019 | 701 |
| 正/负链 | 924 / 95 | 553 / 148 |
| 查询覆盖总和 | 4.44 Mb（95.7%） | 4.63 Mb（99.7%） |
| 最大块 | 58 kb | 108 kb |

差异符合 v1 预期：固定 k=40 精确种子在 ~5% 分歧区段断链，FastGA 的
adaptamer（lcp 扩展）+ wave 对齐补平间隙。两者 PSL 均可被
`pgr psl to-chain` 正常消费（pgr 2038 行 / FastGA 12007 行 chain）。

## 5.2 v2 局部扩展实测（2026-08-02，banded SW + `--ref-seq/--query-seq`）

- 1019 块中 **1015 块被 banded 局部比对细化**（>30 kb 链回退为块，4 块）；
- 扩展身份率 **98.41%**（4.20M match / 68k mismatch），FastGA 为 **97.83%**
  （4.52M / 100k）——同区间参数下高度一致（差异来自块边界与覆盖范围）；
- 并行化：扩展链 rayon `par_iter`，1019 块 20.5s → 2.0s（32 核），
  输出与单线程逐字节一致。

**v2 边界与已知取舍**：

- >30 kb 链不扩展（内存/时间上限），保持块输出——大共线性区块的真实 CIGAR
  留待 wavefront 或分块扩展；
- banded 局部比对用线性 gap（`AlignmentParams` 默认 +5/-4/-8/-6 的 open 部分），
  CIGAR 的 gap 精度不如 affine/wave；
- 块内多段同源时只取最佳局部段（FastGA 同语义）。

## 5.3 株系验证（2026-08-02，MG1655 vs 三株）

| query | pgr 扩展身份率 | FastGA 身份率 | 结构 |
|---|---|---|---|
| nissle1917（Nissle，重排株） | 97.62% | 97.09% | 双方均为 ~30–70 kb 共线性碎片块（大规模倒位） |
| sakai（O157:H7） | 98.41% | 97.83% | 最大块 pgr 58 kb / FastGA 108 kb |
| ec958（UPEC） | 97.60% | — | 最大块 ~30 kb |

pgr 身份率稳定高于 FastGA ~0.5%（banded 局部比对取精确匹配核心；
FastGA 的 wave 延伸进分歧区、覆盖更广）。注意身份率只统计 ≤30 kb
的扩展块（>30 kb 的链回退为无计数块）。

## 6. 相关文档

- 索引格式与消费者规划：[[pbit-index-extension.md]] §3.1
- FastGA 管线与简化移植评估：[[fastga.md]] §11/§12
- 泛基因组场景：[[ecoli-cohort.md]]、[[paf-pangenome.md]]
