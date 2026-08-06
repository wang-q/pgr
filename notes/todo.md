# pgr 近期待办

> 依据 `project-understanding.md`、各设计笔记与审计记录整理。
> 功能层基本齐备，近期大头是验证与数据驱动的扩展。

## 1. 手头数据就能做

- [x] `align hybrid` 三路验证（2026-08-06，`scripts/verify-hybrid-real.sh`）：
      MG1655 vs Sakai 实测覆盖 pgi 90.74% / hybrid **93.08%** / lastz 93.11%，
      耗时 1.2 s / 116.6 s / 135 s——补集方案覆盖几乎追平 lastz（差 0.03 pp）、
      无碎片化（MAF 565 块 < pgi 582）；模拟灵敏度 hybrid 256/600 = lastz
      （补集语义下 gapfill 只是子集，作者 2026-08-06 定稿，实现为
      `compute_holes` + hole × 全 query LASTZ job，`--min-gap` 移除、
      `--max-gap` 默认不限制）（来源：`design/pgi-lastz-hybrid.md` §5.2）。
- [x] `align hybrid` 拆分 → **已实现**（2026-08-06）：`pgr align fill`（2D
      gap fill）+ `pgr align rest`（两侧 trim→excise→holes 一维补集，ref
      holes × 整套 query holes 多序列）+ 复用 `pgr psl lift`（含
      `chr(+):` range 名修复）；**rest 已加采样预筛配对**（默认 syncmer
      s17/w5/ms1：rest 6.5 s/91.81%（PSL 并集口径；**chainnet --syn
      MAF 口径下预筛 vs 全量差仅 0.012 pp**，syntenic 场景无实际损失），
      比全量快 3.7×；--sampler none 全量 92.20%，--smer 15 折中
      92.00%）；集成
      测试 fill 7 例 + rest 4 例，全量 1341 测试通过；模拟灵敏度
      rest 255/600 ≈ lastz 256/600；
      脚本 `scripts/verify-align-fill-rest.sh`（五路，含 fill+rest 组合）+ 改造后的
      `verify-hybrid-sensitivity.sh`（来源：`design/pgi-lastz-hybrid.md` 前半）。
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

- [ ] `align lastz --lastz-args` 帮助文本提示 `--lastz-args=<val>` 写法（值以 `-` 开头时）。
- [ ] `align pgi --self` 校验 `--ref-seq` / `--query-seq` 一致（现在只校验索引 ref==query）。
- [ ] PAF 输出 `query_length` / `target_length` 恒 0：需改索引格式持久化 `src_size`，
      属跨格式变更，明确决策做不做（来源：`audit/audit-paf.md`）。
- [ ] `runlist split -o stdout` 丢弃键名，无法区分同值不同键。
- [ ] `run_lastz` self 模式 n×n job 列表，大目录可提前过滤。
- [ ] `syncmer.rs` 参考实现与 `collect_one_contig` 重复发射同一位置，消费方已去重，可合并。
- [ ] `fa size --no-ns` 把非 IUPAC 字符计为"有效碱基"。

## 4. 技术债（有空再议）

- [ ] HV 矢量化提速不明显（作者长期疑虑，2026-08-06 记录）：`hash_hv_bit` /
      `hash_hv_i8` 用 SIMD 后相对标量只有 ~1.5–2×，8-lane 向量理论上应有更高
      收益——疑点可能在 RNG 生成、字节→向量转换或内存带宽，待深挖
      （基准见 `benchmarks/bench-simd-hv-jaccard.md` §2）。
- [x] nightly `portable_simd` 依赖 → 已迁移 `wide` 1.6.0（2026-08-06）：std::simd 仍未稳定
      （tracking issue #86656 未完成），wide 已 1.x（~2590 万下载、284 个运行时依赖，
      维护活跃）；`hv.rs` 位操作（u32x8 shift/and、u8→i8→i32 转换）与 `linalg.rs`
      （f32x8 reduce_sum→reduce_add、simd_min/max→min/max）已逐项核对迁移，
      `rust-toolchain.toml` 切 stable 1.97；顺手清理 stable 新增 clippy lint 12 处
      （byte_char_slices / collapsible_match）；hnsm 三个基准已迁移并复跑
      （norm ~7.8×、HV i8 快于 bit ~1.5×、rapidhash Jaccard 最快，见
      `benchmarks/bench-simd-hv-jaccard.md`；迁移中修复 hash_hv_i8 标量转换
      导致的 4.3× 退化）（来源：`project-understanding.md` §8.1）。
- [ ] 命令树三跳 dispatch 宏简化，防新增命令漏注册（来源：§8.3）。
- [ ] `fas` 模块职责过重（20 子命令），`fas multiz` 等复杂逻辑考虑拆分（来源：§8.6）。
- [x] 分层一致性：commit `d5281bc` 将 hybrid 逻辑迁入 `cmd_pgr/align/hybrid.rs`，
      与 AGENTS.md"复杂逻辑放 libs"相悖——2026-08-06 作者确认**有意为之**
      （消除跨模块导入开销，代码仅此命令使用），维持现状不回迁；此条不再作为技术债。

## 5. 明确不做（避免重复立项）

- Gap_Improver、完整 LCP、`.1aln`、trace points、ALNchain、GDB/GIX 分片：已裁定不做
  （来源：`design/pgi-align.md` §6）。
- 多 mask union：暂缓，低优先便利（§7.5）。
- `-S` 对称 adaptamer：默认不做，仅场景开关候选（§7.4.1）。
