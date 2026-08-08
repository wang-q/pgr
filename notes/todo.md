# pgr 近期待办

> 依据 `project-understanding.md`、各设计笔记与审计记录整理。
> 功能层基本齐备，近期大头是验证与数据驱动的扩展。
> 已完成条目只留一行结论，细节见链接文档。

## 1. 手头数据就能做

- [x] **dist 命令族拆分 + Mash 兼容收尾**（2026-08-08）：`dist seq` 删除，
      拆为 mini/mash/frac；denom/containment 语义、流式内存、空集、文件级
      并行 + bottom-k 预筛全部对齐 Mash（详见
      `benchmarks/bench-dist-mash-compat.md`、docs/dist.md）。
- [x] **FracMinHash 采样器落地**（2026-08-08）：独立命令 `dist frac` + `--ci`，
      无偏验证 Spearman 1.0（详见 `benchmarks/dist-cohort-validation.md`）。
- [x] **pgi 引擎 SD 灵敏度优化**（2026-08-06）：默认 freq=50/k=31，漏检率
      13.1%→0.26%（详见 `design/sd.md` §4.9）。
- [x] **repeat masking 闭环**（2026-08-06）：先遮蔽再 `sd search --engine pgi`
      （详见 `design/sd.md` §4.10）。
- [x] **tube 工作流失效重测**（2026-08-06）：未复现，无需改动
      （来源：`audit/audit-rept-sd.md`）。

## 2. 等数据/场景到位再启动

- [x] **syncmer 采样位置偏差**（2026-08-08 量化收尾）：pgi（k=40）与
      `dist hv --sampler syncmer`（s-mer=8）同为 syncmer 位置漂移导致的交集
      低估，数值不可信、排序弱于 frac；**不可修复**（pgi syncmer 是 align
      锚点基石）。分工：pgi+syncmer=比对/排序、数值 ANI=`dist frac`；
      体验入口 = `dist hv --sampler syncmer` / `dist pgi`
      （详见 `benchmarks/dist-cohort-validation.md`）。
- [ ] 稀疏 s=1 完整 45 对 cohort 复测：等 10 株数据（5 株 10 对已验 ρ=0.988）。
- [x] **FracMinHash containment 偏差评估**（2026-08-08）：原 10% 偏差是
      k 不匹配假象，同 k 偏差 1–2%（采样方差）；Hera 校正对大基因组无效，
      **不实现**（详见 `benchmarks/bench-dist-mash-compat.md` 同 k 对照节）。
- [x] **`.hv` 稀疏选择**（2026-08-08 已落地 s=1）：理论补齐 + cohort 复测
      通过；剩余实施方案定夺（用户）（详见 `design/hv.md` §2.7）。
- [x] **FASTA `dist hv` 量纲问题**（2026-08-08 已修复）：改用 `hash_hv_bit`
      （详见 `design/hv.md` §3.4）。
- [ ] 4 万 E. coli cohort 端到端：核心步骤就绪，等真实 cohort 数据
      （来源：`ecoli-cohort.md`）。
- [ ] 人类规模（GRCh38/CHM13）验证：`.pgi` 字段上限、内存/耗时与 FastGA
      对照（来源：`design/pgi-align.md` §7.2）。
- [ ] pbit 自动路由：等多样性 cohort 数据证明收益（来源：`design/pbit.md`）。
- [ ] pbit HV sketch 内嵌（决策 B）：触发 = 无源 FASTA 归档需距离粗筛。
- [ ] `--sym` 场景开关：先量化方向偏差再实现（来源：`design/pgi-align.md` §7.4.1）。
- [ ] 完整 adaptamer（变长种子 >k）：前置 lcp 已落地，只差立项
      （来源：`design/pgi-query-layer.md`）。
- [ ] `dist mash` 序列级并行：等单文件多 contig 大规模场景（文件级并行
      已覆盖多文件，见 `benchmarks/bench-dist-mash-compat.md` 性能节）。

## 3. 低风险审计记录项（可顺手修）

- [x] minimizer 默认参数对 DNA 不充分（2026-08-08 已修复）：DNA 默认
      k=21/w=5（`resolve_kmer_window`）。
- [ ] PAF `query_length`/`target_length` 恒 0：**待决策**——改 `.paf.idx`
      持久化 src_size，影响索引格式兼容性（来源：`audit/audit-paf.md`）。
- [ ] `syncmer.rs` 重复发射同一位置：**暂缓**——消费方已去重，收益小风险高
      （来源：`audit/audit-rept-sd.md`）。

## 4. 技术债（有空再议）

- [x] HV 矢量化提速不明显（2026-08-07 已解决）：AVX2 为主，bit ±1 ~4.8×
      （详见 `benchmarks/bench-simd-hv-jaccard.md`）。
- [ ] `fas` 模块职责过重（20 子命令）考虑拆分（来源：§8.6）。

## 5. 明确不做（避免重复立项）

- Gap_Improver、完整 LCP、`.1aln`、trace points、ALNchain、GDB/GIX 分片
  （`design/pgi-align.md` §6）；多 mask union（§7.5）；`-S` 对称 adaptamer
  （§7.4.1）；hybrid 逻辑留 `cmd_pgr/`（commit `d5281bc`，有意为之）。
- sd 边界扩展（BISER MAX_EXTEND）：已评估不做（`design/sd.md` §4.8/§4.9）。
- `dist mash` murmur3 SIMD：已评估不做（1.3× 差距可接受，风险大于收益）。
- 命令树 dispatch 宏简化：已评估不做（用户裁定——宏抹消对注册/分发代码
  的显式理解；`tests/cli_consistency.rs` 已核对注册一致性）。
