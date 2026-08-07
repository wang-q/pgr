# pgr 近期待办

> 依据 `project-understanding.md`、各设计笔记与审计记录整理。
> 功能层基本齐备，近期大头是验证与数据驱动的扩展。
> 已完成的条目（fill/rest 拆分、wide 迁移、审计项批量修复等）详见
> `design/pgi-lastz-hybrid.md`、`benchmarks/bench-simd-hv-jaccard.md` 与
> `project-understanding.md`，此处只留未完成项。

## 1. 手头数据就能做

- [x] **FracMinHash 采样器落地**（2026-08-08）：`dist seq --sampler frachash`
      （canonical k-mer 保留 hash < u64::MAX/scale，`--scale` 默认 1000）
      + `--ci` 输出 ANI 95% 置信区间（正态近似）。无偏验证：5 株 × 10 对
      与全 k=40 集合真值排序 **Spearman 1.0**、Jaccard 0.417 vs 真值
      0.451（`dist pgi` 仅 0.095——syncmer 偏差）；MG1655×Sakai ANI
      97.7% vs 真值 97.3%（详见 `benchmarks/dist-cohort-validation.md`）。
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
- [x] tube 工作流"库 vs 基因组"失效重测（2026-08-06）：greedy 已移除、tube
      唯一流程；MG1655 vs TnCentral `rept e-align` 正常检出（71.6 kb，
      79% 与 e-kmer 重叠），**失效未复现**（随 syncmer/排序键修复消失），
      无需改动（来源：`audit/audit-rept-sd.md`）。

## 2. 等数据/场景到位再启动

- [ ] `dist pgi` 的 syncmer 8/5 采样位置偏差（2026-08-08 发现，
      `benchmarks/dist-cohort-validation.md`）：近缘/含共享重复元件的
      基因组间 Jaccard 严重低估（e2348×cft073：0.095 vs 全 k=40 真值
      0.451）。**注意（2026-08-08 用户提醒）**：这是采样集合交集的固有
      行为，不是可"修复"的 bug——pgi 的 syncmer 是 align pgi（FastGA
      兼容比对）的锚点基石，不能换 FracMinHash（双选概率 1/scale、
      scale=1000 时锚点只剩 1/1000 且无窗口保证，链化失效）。正确分工：
      pgi+syncmer = 比对/排序；**无偏数值 ANI 用 `dist seq --sampler
      frachash`**（已实现 + CI）。补充（2026-08-08）：containment 略稳
      于 jaccard（相对偏差 34% vs 39%、排序 ρ 0.66 vs 0.52，但仍不可靠）；
      `.pgi` 约 FASTA 的 27×（37.6MB vs 1.4MB），**距离计算不应为它建
      索引**（初筛用 `.hv`、数值用 frachash，见 docs/dist.md 分层建议）。
- [ ] 稀疏 s=1 的完整 45 对 cohort 复测（2026-08-08 仅 5 株 10 对验证
      ρ=0.988；原 s=3 是 10 株 45 对）：等 10 株数据到位后补全，确认
      s=1 排序一致性（理论 + 10 对实测已支持，此为收尾验证）。
- [ ] FracMinHash containment/ANI 偏差校正（2026-08-08 立项，等 Hera
      论文：Hera et al. 2023, *Genome Res* 33(7):1061–1068,
      doi:10.1101/gr.277651.123）：当前 `dist seq --sampler frachash`
      的 containment 有 ~10% Jensen 偏差（实测与 scale 无关、增大采样
      无效）；实现 Hera 校正公式消除到一阶/二阶，使 containment/ANI
      数值精确（来源：`benchmarks/dist-cohort-validation.md`
      FracMinHash containment 小节）。
- [x] 重新审视 `.hv` 路径的稀疏选择（2026-08-08 已分析 + s=1 已落地）：
      稀疏是历史性能优化产物（非有意设计）但理论已补齐（§2.7/§6.5/§6.6）；
      权衡 = 编码瓶颈场景用稀疏 s=1（~50× 快 + 大 D 免费精度）、小规模
      用稠密 bit；s=1 已落地且 cohort 复测（5 株 10 对 Spearman 0.988）
      通过。剩余：实施方案层面的最终定夺（用户）。
- [x] FASTA `dist hv` 量纲问题（2026-08-08 已修复，`design/hv.md` §3.4）：
      `load_hv_from_fasta` / `load_hv_from_fasta_syncmer` 改用 `hash_hv_bit`；
      模拟 Jaccard 0.102 vs 真值 0.091、两株 E. coli 实测输出合理
      （测试 `test_hash_hv_bit_jaccard_accurate`）。
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

- [x] `dist seq` / `dist hv` 的 minimizer 默认参数对 DNA 不充分
      （2026-08-08 已修复）：默认 k=7 的 2-bit 编码仅 2^14 空间 →
      unique minimizer ≈ 16383（4.6Mb 基因组饱和），采样严重不足；
      `resolve_kmer_window` 现对 DNA minimizer 默认 k=21/w=5（与文档
      docs/dist.md 一致），实测 MG1655×Sakai mash 0.023 → ANI 97.7%
      （真值 97.3%）。
- [ ] PAF 输出 `query_length` / `target_length` 恒 0：需改索引格式持久化 `src_size`，
      属跨格式变更，**待决策**：改 `.paf.idx` 持久化 src_size 才能填充，
      影响索引格式兼容性（来源：`audit/audit-paf.md`）。
- [ ] `syncmer.rs` 参考实现与 `collect_one_contig` 重复发射同一位置，消费方已去重，
      可合并——**暂缓**：涉及 pgi build 种子发射核心，收益小风险高，消费方已去重
      （来源：`audit/audit-rept-sd.md`）。

## 4. 技术债（有空再议）

- [x] HV 矢量化提速不明显（2026-08-07 已解决）：`hash_hv_bit` / `hash_hv_i8`
      以 AVX2（256-bit）为主实现（跳步 RapidRng + 块主序 + 位展开），
      bit ±1 编码实测 ~4.8×（相对旧 bit）/ ~3.1×（相对旧 i8）、
      i8 保语义 ~2.1×；AVX-512 仅作基准参考不参与分派，无 AVX2 自动降级
      （基准见 `benchmarks/bench-simd-hv-jaccard.md` §2/§5）。
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
