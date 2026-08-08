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
- [ ] 物种内聚类选参考 + PBit 归档（核心用例，UI 待讨论）：方法调研 +
      场景工作流已记录（`design/genome-nn-query.md` §5）：dist 输出
      pair TSV → **Necom**（`~/Scripts/necom`，聚类/构树/剪枝）→ 选参考
      → align pgi → pbit create；缺口 = pgr 侧 `dist` 输出与 Necom 格式
      对齐、参考挑选、pbit 自动路由（`design/pbit.md` 决策点 1 的触发
      场景）；HV 最近邻存储方案已调研（§6：≤10 万 SQLite + sqlite-vec，
      或零依赖自研 SIMD 扫描；usearch 为 C++ FFI 负担大暂不引入）；ANN
      ANN 召回已实测（rust-cv `hnsw`：`benchmarks/bench-hv-ann-recall.md`；
      GSearch 同源 `hnsw_rs` + HubNSW 单层：`benchmarks/bench-hv-ann-hubnsw.md`，
      2026-08-08）：rust-cv 在 30k 召回上限 0.92 但快 20–46×；hnsw_rs
      召回 0.99 但只快 2.3–6.5×，HubNSW 单层仅微改善；结论 = ≤10k 精确
      扫描即可，10k–30k 依召回/速度偏好选实现，>30k 先降维再评估 ANN；
      sqlite-vec 4096 维延迟与降维路线待实测）。
- [x] **HV/Mash vs ANI 标定（P1 主体，2026-08-08）**：135 个真实基因组
      （E. coli NR + 其他 Escherichia + Yersinia，全部挂靠 pass.lst/NR.lst/
      genome.taxon.tsv）上 HV（D=4096/16384）与 Mash vs skani ANI 的
      Spearman/RMSE/recall@10（详见 `benchmarks/bench-hv-ani-calibration.md`）。
      结论：HV 仅 ANI 90–98% 中等可靠，≥99%/＜85% 失效，D 不救近缘；
      Mash 同种内 ρ=−0.97。
- [x] **标定剩余（P1，2026-08-08）**：frac/mini 同 cohort 标定
      （frac≈Mash，mini≈HV 近缘失效）；完整度鲁棒性（近缘对 50% 完整度
      HV +43%/Mash +84%，ANI 稳定）；sampler/k/D 扫描（frac 默认参数
      已近上限，mini k/w/hasher 影响小）；组装质量（碎片化驱动误差）
      （`benchmarks/bench-parameter-scan.md`、`bench-completeness-robustness.md`）。
- [x] **ANI 真值下的图检索召回（P1，2026-08-08）**：真实 HV 上全局
      HNSW recall_HV≥0.993（图检索误差可忽略），recall_ANI 0.664=精确；
      物种硬路由在小 clade 上反而有害（R=1 跌到 0.70），clade 需 ≥K
      成员（`benchmarks/bench-hv-ann-real.md`）。
- [ ] **标定/检索剩余**：sqlite-vec 真实 HV 延迟（等安装）；E. coli NR
      全量（2,115）与全 NR（15,574）实跑（已 494 规模外推）；pbit
      多参考/高分歧样本验证；bac120 标记基因路由准确率（§7.4 #10/#9/#14/#19）。
- [x] **bac120 标记基因路由（#19，2026-08-08）**：8 个保守标记蛋白 aa
      最近邻路由准确率 0.756（ANI 金标准上限 0.800，HV 路由 0.822）
      （`benchmarks/bench-marker-routing.md`）。
- [ ] **pbit PAF 闭环（#14，2026-08-08 部分）**：已修 `--paf` 跨组装
      命名 bug（回归测试过）；还需实现 `psl to-paf` 的 cg:Z（=/X）输出
      或接入 minimap2，才能拿到真实压缩率（约束见
      `benchmarks/bench-scale-and-pbit.md` #14b）。
- [ ] **聚类/选参考/PBit 端到端验证（P2）**：Necom 聚类 vs GTDB 标签
      一致率；参考策略（中心/最长/随机）→ pgi → pbit 压缩率对比
      （§7.2③–⑤，数据需求：一个物种的真实 cohort）。

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
