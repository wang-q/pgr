# 规模延迟实测（#8）与 PBit 端到端试点（#14）

> 日期：2026-08-08。对应 `design/genome-nn-query.md` §7.4 #8/#14。

## #8 E. coli NR 规模：真实 HV 检索延迟

方法：E. coli NR 抽 494 个基因组（每 4 个取 1，跳过缺失组装），
`pgr pgi build + to-hv`（syncmer 40/8/5，D=4096），L2 归一化；
`bench hv_ann_real` 测精确扫描与 hnsw_rs HNSW 的延迟与 recall_HV@10
（无 ANI 真值，只看 HV 图检索层）。构建耗时 ≈ 0.7 s/基因组（含 pgi）。

| 规模 | 精确扫描 µs/查询 | HNSW ef=10 | HNSW ef=50 |
|---|---|---|---|
| 135 | 331 | 177 µs（recall_HV 0.993） | 312 µs（0.996） |
| 494 | 1,165 | 306 µs（0.985） | 616 µs（0.997） |

线性外推：2,115 NR ≈ 精确 5 ms、HNSW 1.3–2.6 ms/查询；15,574 全 NR
≈ 精确 37 ms、HNSW 10–20 ms。内存：4096 维 i32 HV = 16 KB/基因组，
15k ≈ 240 MB 向量 + 图边。结论：**≤10k 规模精确扫描仍在 ms 级，
维持 §6.4 建议**；HNSW 在 500 规模给出 ~3.8×（ef=10）且不损失 HV
排序（recall_HV≥0.985）。

## #14 PBit 试点：1 参考 + 5 近缘样本 E. coli

方法：`pgr pbit create -r ref -i 5 samples`（默认 segment 4096、
k=15、min-match 18）。样本为 NR 近缘株（ANI>99%）。

> **重要更正（2026-08-08）**：首次试点报出的"边际 delta ≈30 B/株、
> ~55,000× 压缩"是**假象**——`pbit create` 按 **contig 名**匹配参考，
> 而样本是 draft 组装（NZ_* contig 名），与参考染色体名不同，导致
> **全部 contig 被跳过**（create 日志逐条 WARN "not found in reference;
> skipping"），归档只存了参考，`pbit to-fa` 重建为 0 bp。因此 naive
> create 对跨组装样本**不存任何样本内容**，压缩比无意义。

正确用法：样本须经比对（PAF，CIGAR 编码）或 contig 名与参考一致。
核心场景（align pgi → pbit）中 PAF 是**必需**输入，不是可选项。

| 指标 | 值 |
|---|---|
| naive create（无 PAF） | 归档≈参考 2bit（1.36 MB）；样本全跳过，to-fa 重建 0 bp |
| create --paf（需 cg:Z PAF） | 见下节：bug 已修 + 约束明确；真实压缩率待完整管线 |

结论（修正后）：**pbit 的压缩前提是样本已被比对到参考**（PAF/CIGAR
或同名 contig）；对 draft 组装直接 create 会静默丢数据。核心场景
（聚类选参考 → align pgi → pbit）中比对步骤不可省略，且归档前应
校验样本重建覆盖率（`pbit to-fa` + 长度对比）作为质量门。

### #14b PAF 路径深挖：一个真 bug + 三条设计约束（2026-08-08）

**Bug（已修 + 回归测试）**：`append_sample_with_paf` 在尝试 CIGAR 编码
前先用**样本 contig 名**查**参考**键表（`contig_ref_groups`），跨组装
命名（draft NZ_* vs 参考）的样本永远被跳过，PAF 白给。修复：CIGAR
优先（走 PAF 索引，不要求名字一致），LZ-diff 兜底才按名字匹配；跳过
时不再注册空 contig。新增回归测试
`test_append_sample_with_paf_cross_assembly_names`，全部 617 个测试通过。

**约束 1（集成缺口）**：pbit 需要带 `cg:Z`（且为 `=/X/I/D` 语义，M 被
拒）的 PAF；`pgr psl to-paf` 只输出 12 列、无 cg:Z——pgr 自己的
align→pbit 闭环缺一环（需实现 PSL→PAF+cg 转换，块内按序列比对生成
`=/X`）。

**约束 2（段覆盖）**：CIGAR 编码要求每个 4096 bp 查询段被**单条 PAF
记录全覆盖**；draft 组装短 contig/碎片化比对（单条记录只盖几百 bp）
全部落回 LZ（名字不匹配则跳过）→ 重建 0。

**约束 3（段内目标）**：段的 target 区间必须落在**单个**参考 4096 段内
（`t_seg_idx_start == t_seg_idx_end`）；完整 E. coli 基因组间的重排/
大 gap 使 pgi 块散布，绝大多数段的目标跨边界 → 被拒。实测两个完整
E. coli（ANI≈99.9%）pair：1159 条 PSL 链去重后仍大量重叠，贪心单调链
生成的长 CIGAR 中 D gap 达 3 Mb，重建仍为 0。

**对设计的含义**：pbit 的 CIGAR 路径面向**近同、共线性、无大重排**的
样本（完整染色体 vs 完整染色体）；真实菌株间的反转/质粒差异会显著
降低覆盖率。选参考时应优先与样本共线性的参考；归档质量门
（`pbit to-fa` 覆盖率）应进工作流。真实压缩率数字需要先把约束 1 的
转换器实现（或接入 minimap2），列为后续任务。

## 复现

```bash
# #8: PGR_HV_REAL_DIR=/tmp/hv_calib/hv500 ... cargo bench --bench hv_ann_real
# #14: pgr pbit create -r ref.fa -i s1.fa ... -o a.pbit
#      用 1/2/5 样本归档大小差算边际成本
```
