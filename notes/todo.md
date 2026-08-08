# pgr 近期待办

> 依据 `project-understanding.md`、各设计笔记与审计记录整理。
> 功能层基本齐备，近期大头是验证与数据驱动的扩展。
> 已完成条目只留一行结论，细节见链接文档。

## 0. 会话交接（2026-08-09，pbit v1007–v1010 + 主链/碎链判定）

**已完成（本会话）**：

- pbit 编码演进 v1006–v1010：LZ 内容匹配化 → 跨相位 CIGAR（段级混合）→
  `to-paf` 无损还原输入 PAF（809/809 逐字段一致）→ Identity 零载荷，
  详见 `design/pbit.md` §PAF 驱动编码的演进 + `bench-scale-and-pbit.md`
  #14f/#14i–l；
- **主链/碎链判定（定稿 + 落地）**：链级贪心（query 覆盖段数 → 相似度 →
  输入序，`BIG_CHAIN_MIN_LEN` 退役）；碎链也可编码其覆盖的段；无 `cg:Z`
  记录存行还原；真实数据还原 100%、编码分布与旧行为一致（详见
  `design/pbit.md` §大链与碎链）；
- **`--paf` 强制（A）**：CLI 校验已实现（`collect_samples_from_args` 拒绝
  无 PAF 样本，空 PAF 可禁用 CIGAR）；无 PAF 独立路径退役；
- **设计基础 7 条定稿**（样本复用参考 / 参考坐标系 / 无损 / PAF 存储 /
  可计算的不存储 / PAF 路由载体 / 匹配区间可访问性），见 `design/pbit.md`
  顶部。

**未提交**：pbit 相关代码、测试与文档（本会话改动，含 v1010 判定 + 强制
PAF）；并行审计另有 6 文件（fas_xlsx.rs、fas/subset.rs、cli_fas_vars.rs、
paf-pangenome.md、project-understanding.md、impg.md），**勿碰**。提交前跑
fmt/clippy/test 确认。

**挂账/待决**：

1. **pgi 长链链化**——长期挂账（pbit 路线 3，依赖对齐器）。
2. **碎链 cg 位打包**——暂缓（用户明确：别纠结碎链，整行压缩量可接受）。
3. ~~`--paf` 强制 / Identity 优化 / AGC 分析文档化 / 大链碎链判定~~ →
   已完成（见上；cg:Z 生产者明确不做，todo §4）。

## 1. 等数据/场景到位再启动

- [ ] 稀疏 s=1 完整 45 对 cohort 复测：等 10 株数据（5 株 10 对已验 ρ=0.988）。
- [ ] 4 万 E. coli cohort 端到端：核心步骤就绪，等真实 cohort 数据
      （来源：`ecoli-cohort.md`）。
- [ ] 人类规模（GRCh38/CHM13）验证：`.pgi` 字段上限、内存/耗时与 FastGA
      对照（来源：`design/pgi-align.md` §7.2）。
- [ ] pbit 自动路由：等多样性 cohort 数据证明收益（来源：`design/pbit.md`）。
- [ ] `--sym` 场景开关：先量化方向偏差再实现（来源：`design/pgi-align.md` §7.4.1）。
- [ ] 完整 adaptamer（变长种子 >k）：前置 lcp 已落地，只差立项
      （来源：`design/pgi-query-layer.md`）。
- [ ] `dist mash` 序列级并行：等单文件多 contig 大规模场景（文件级并行
      已覆盖多文件，见 `benchmarks/bench-dist-mash-compat.md` 性能节）。
- [ ] 物种内聚类选参考 + PBit 归档（核心用例，UI 待讨论）：场景工作流见
      `design/genome-nn-query.md` §5；缺口 = `dist` 输出与 Necom 格式对齐、
      参考挑选、pbit 自动路由（`design/pbit.md` 决策点 1 触发场景）；
      HV 最近邻/SQLite/ANN 调研结论见 §6 与 `bench-hv-ann-*.md`、
      `bench-scale-and-pbit.md` #10b。
- [ ] 标定/检索剩余：E. coli NR 全量（15,574）实跑（2,088 已实测，
      `bench-scale-and-pbit.md` #8b）；pbit 多参考/高分歧样本验证（#14 路线）；
      §7.4 #10/#9/#14/#19 状态见 `genome-nn-query.md`。

## 2. 低风险审计记录项（可顺手修）

- [ ] PAF `query_length`/`target_length` 恒 0：**待决策**——改 `.paf.idx`
      持久化 src_size，影响索引格式兼容性（来源：`audit/audit-paf.md`）。
- [ ] `syncmer.rs` 重复发射同一位置：**暂缓**——消费方已去重，收益小风险高
      （来源：`audit/audit-rept-sd.md`）。

## 3. 技术债（有空再议）

- [ ] `fas` 模块职责过重（20 子命令）考虑拆分（来源：`genome-nn-query.md` §8.6）。

## 4. 明确不做（避免重复立项）

- Gap_Improver、完整 LCP、`.1aln`、trace points、ALNchain、GDB/GIX 分片
  （`design/pgi-align.md` §6）；多 mask union（§7.5）；`-S` 对称 adaptamer
  （§7.4.1）；hybrid 逻辑留 `cmd_pgr/`（commit `d5281bc`，有意为之）。
- sd 边界扩展（BISER MAX_EXTEND）：已评估不做（`design/sd.md` §4.8/§4.9）。
- `dist mash` murmur3 SIMD：已评估不做（1.3× 差距可接受，风险大于收益）。
- 命令树 dispatch 宏简化：已评估不做（用户裁定——宏抹消对注册/分发代码
  的显式理解；`tests/cli_consistency.rs` 已核对注册一致性）。
- **pbit HV sketch 内嵌（决策 B）：明确不做**（2026-08-09 用户裁定——HV
  评测未达预期，后续换其他形式；原暂缓触发条件 = 无源 FASTA 归档需距离
  粗筛，见 `design/pbit.md` 决策 B）。
- **链级 `cg:Z` 生产者：明确不做**（2026-08-09 用户裁定——推荐链路
  `chainnet → maf to-paf` 自带 cg:Z；`psl to-paf` 无 cg:Z 记录走
  "跳过编码、存行还原"）。
