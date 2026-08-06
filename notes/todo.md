# pgr 近期待办

> 依据 `project-understanding.md`、各设计笔记与审计记录整理。
> 功能层基本齐备，近期大头是验证与数据驱动的扩展。
> 已完成的条目（fill/rest 拆分、wide 迁移、审计项批量修复等）详见
> `design/pgi-lastz-hybrid.md`、`benchmarks/bench-simd-hv-jaccard.md` 与
> `project-understanding.md`，此处只留未完成项。

## 1. 手头数据就能做

- [x] **pgi 引擎 SD 灵敏度优化**（2026-08-06 已实现）：`sd search` 透传
      pgi 参数，默认 `freq=50, kmer=31`——10 个 E. coli 整体漏检率
      **13.1% → 0.26%**（e2348_69 高拷贝重复 562→0、sakai/e24377a 的
      90–93% 分歧漏检归零）；剩余 11 个边缘个案（低复杂度结构）记为
      已知限制（详见 `design/sd.md` §4.9）。
- [x] repeat masking 闭环（2026-08-06）：e2e 回归已固化 + **遮蔽版验证完成**——
      10 个 E. coli 基因组跑 `rept e-kmer` 三库 → `fa mask --hard` → `sd search`：
      遮蔽 ~1.2%（IS 元件主导），遮蔽后 pgi/lastz 引擎 SD 检出收敛
      （互相漏检 3.2%/6.0%，未遮蔽时 pgi 漏检 13.1%）；标准 SD 流程 =
      先遮蔽再 `sd search --engine pgi`（默认 freq=50/k=31）
      （详见 `design/sd.md` §4.10）。
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
- sd 边界扩展（BISER MAX_EXTEND 移植）：**已评估后不做**（2026-08-06，与 lastz
  比较：pgi hit 左边界仅短 2–6 bp、右边界一致，收益边际；且灵敏度优化
  freq=50/k=31 已解决更实质的漏检问题——详见 `design/sd.md` §4.8/§4.9）。
