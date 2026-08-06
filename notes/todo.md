# pgr 近期待办

> 依据 `project-understanding.md`、各设计笔记与审计记录整理。
> 功能层基本齐备，近期大头是验证与数据驱动的扩展。
> 已完成的条目（fill/rest 拆分、wide 迁移、审计项批量修复等）详见
> `design/pgi-lastz-hybrid.md`、`benchmarks/bench-simd-hv-jaccard.md` 与
> `project-understanding.md`，此处只留未完成项。

## 1. 手头数据就能做

- [ ] repeat masking 闭环：e2e 回归已固化（`tests/cli_rept.rs` 的
      `command_rept_e_kmer_end_to_end`，用 `tests/pgr/tncentral.fa.gz`，已实测通过）；
      三库备于 `~/data/repeats/`，本机手跑示例见 `docs/rept.md`（大库不进 CI）；
      剩余：遮蔽版（`pgi build --mask`）按 `design/repeat-masking.md` §2.5 验证。
- [ ] sd 边界扩展（可选）：BISER MAX_EXTEND 未移植，`sd search` 块边界比真实短 1–11 bp，
      检出后向两侧扩展再进 chain/net（来源：`design/sd.md` §4.8）。
- [ ] tube 工作流"库 vs 基因组"失效重测：原结论基于修复前代码，syncmer/排序键修复后
      用真实数据重测（来源：`audit/audit-rept-sd.md`）。

## 2. 等数据/场景到位再启动

- [ ] 4 万 E. coli cohort 端到端：pgr 核心步骤（PAF 索引/查询/图）已就绪，等真实 cohort 数据；
      到位后重跑 4 万规模基准，按 `paf-pangenome.md` §5.2 判断标准选优化项，再做应用层。
      上游去冗余/sparsify 在 pgr 外（Mash + FastGA）；远期可封装 `pgr pl dedup` / `sparsify`
      （来源：`ecoli-cohort.md`）。
- [ ] 人类规模（GRCh38/CHM13）验证：核对 `.pgi` 字段上限、构建/比对内存与耗时，与 FastGA
      对照；同时决定 `--sym` 在真核场景是否值得做（来源：`design/pgi-align.md` §7.2）。
- [ ] pbit 自动路由：留待多样性 cohort 数据证明收益（来源：`design/pbit.md` 顶部决策）。
- [ ] pbit HV sketch 内嵌（决策 B）：触发条件 = "无源 FASTA、仅归档、需距离粗筛"的真实工作流。
- [ ] `--sym` 场景开关：先量化方向偏差（45 对 cohort 抽 5–10 对跑 (A,B) vs (B,A)），
      并确认 `sd cross` 是否需要双向全量；量化完再实现（来源：`design/pgi-align.md` §7.4.1）。
- [ ] 完整 adaptamer（变长种子 >k）：前置 lcp 已落地，只差立项；当前非优先级
      （来源：`design/pgi-query-layer.md` 目标 4）。

## 3. 低风险审计记录项（可顺手修）

- [ ] PAF 输出 `query_length` / `target_length` 恒 0：需改索引格式持久化 `src_size`，
      属跨格式变更，**待决策**：改 `.paf.idx` 持久化 src_size 才能填充，
      影响索引格式兼容性（来源：`audit/audit-paf.md`）。
- [ ] `syncmer.rs` 参考实现与 `collect_one_contig` 重复发射同一位置，消费方已去重，
      可合并——**暂缓**：涉及 pgi build 种子发射核心，收益小风险高，消费方已去重
      （来源：`audit/audit-rept-sd.md`）。

## 4. 技术债（有空再议）

- [ ] HV 矢量化提速不明显（作者长期疑虑，2026-08-06 记录）：`hash_hv_bit` /
      `hash_hv_i8` 用 SIMD 后相对标量只有 ~1.5–2×，8-lane 向量理论上应有更高
      收益——疑点可能在 RNG 生成、字节→向量转换或内存带宽，待深挖
      （基准见 `benchmarks/bench-simd-hv-jaccard.md` §2）。
- [ ] 命令树三跳 dispatch 宏简化，防新增命令漏注册（来源：§8.3）。
- [ ] `fas` 模块职责过重（20 子命令），`fas multiz` 等复杂逻辑考虑拆分（来源：§8.6）。

## 5. 明确不做（避免重复立项）

- Gap_Improver、完整 LCP、`.1aln`、trace points、ALNchain、GDB/GIX 分片：已裁定不做
  （来源：`design/pgi-align.md` §6）。
- 多 mask union：暂缓，低优先便利（§7.5）。
- `-S` 对称 adaptamer：默认不做，仅场景开关候选（§7.4.1）。
- hybrid 逻辑留在 `cmd_pgr/`（commit `d5281bc`，作者 2026-08-06 确认有意为之，
  消除跨模块导入开销、代码仅该命令使用；不回迁 libs，不再作为技术债讨论）。
