# pgr 近期待办

> 依据 `project-understanding.md`、各设计笔记与审计记录整理。
> 功能层基本齐备，近期大头是验证与数据驱动的扩展。
> 已完成条目只留一行结论，细节见链接文档。

## 0. 会话交接（2026-08-09，pbit v1007–v1010 + 主链/碎链判定）

**已完成（本会话）**：

- pbit 编码演进 v1006–v1010：LZ 内容匹配化 → 跨相位 CIGAR（段级混合）→
  `to-paf` 无损还原输入 PAF（809/809 逐字段一致）→ Identity 零载荷，
  详见 `design/pbit.md` §PAF 驱动编码的演进 + `benchmarks/bench-scale-and-pbit.md`
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
- **POA SIMD 加速（2026-08-09）**：`libs/poa/simd.rs` 垂直并行（lane = 序列
  位置），分派沿用 HV 式（AVX2 手写 + `wide` 回退，无 SSE4.1 中间档/SIMDe）；
  `Poa` 默认引擎切换，标量保留测试对照；基准 120 bp ~8.7× / 600 bp ~12.3×
  （`benches/poa_benchmark.rs`）；参考分析见 `references/spoa.md`（新增）。
- **fa 逐字节统计矢量化（2026-08-09）**：`libs/nt_simd.rs` 新增
  `count_valid`（`fa size --no-ns`）/ `count_n`（N 家族含 IUPAC）/
  `masked_bitmap`（`fa masked` 位图，区域合并扫描在 `fmt/fa.rs`），分派沿用
  HV 式（AVX2 手写 + `wide` 回退，N 家族非 AVX2 回退标量）；基准
  count_valid ~14×、count_n ~6.5×、masked_bitmap ~15×
  （`benches/byte_stat_benchmark.rs`，详见 `benchmarks/bench-nt-simd.md`）。

**未提交**：pbit 相关代码、测试与文档（本会话改动，含 v1010 判定 + 强制
PAF）；并行审计另有 6 文件（fas_xlsx.rs、fas/subset.rs、cli_fas_vars.rs、
paf-pangenome.md、project-understanding.md、impg.md），**勿碰**。提交前跑
fmt/clippy/test 确认。

**挂账/待决**：

1. **碎链 cg 位打包**——暂缓（用户明确：别纠结碎链，整行压缩量可接受）。
2. ~~`--paf` 强制 / Identity 优化 / AGC 分析文档化 / 大链碎链判定~~ →
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
      已覆盖多文件，见 `design/hv.md` 性能节）。
- [ ] 物种内聚类选参考 + PBit 归档（核心用例，UI 待讨论）：场景工作流见
      `design/genome-nn-query.md` §5；缺口 = `dist` 输出与 Necom 格式对齐、
      参考挑选、pbit 自动路由（`design/pbit.md` 决策点 1 触发场景）；
      HV 最近邻/SQLite/ANN 调研结论见 §6 与 `design/hv.md`（证据附录）、
      `benchmarks/bench-scale-and-pbit.md` #10b。
- [ ] 标定/检索剩余：E. coli NR 全量（15,574）实跑（2,088 已实测，
      `benchmarks/bench-scale-and-pbit.md` #8b）；pbit **多参考**样本验证待做
      （高分歧已在 #14g 实测：E. albertii ANI≈90%，delta/gzip 78–81%）；
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
- **pgi 长链链化（pbit 路线 3）：明确不做**（2026-08-09 用户裁定——项目
  优势 = 引入 UCSC chainnet 经典链化管线，自研 chain 效果始终不如它；
  链化依赖由 chainnet 承担，见 `design/pbit.md` §PAF 驱动编码的演进）。
- **gzip 并行解压 / zlib-ng / libdeflate：明确不做**（2026-08-09 用户裁定——
  程序常被 shell 包裹并行执行，pgr 侧 `fa` 保持单线程；inflate 内部已是
  zlib-rs AVX2，见 `benchmarks/bench-profile-hotspots.md` 场景 1）。

## 5. 待实现 / 待决策（2026-08-09 文档扫描补充）

- [x] ~~**`spanr cover` 名字截断问题**~~ → **已完成（代码核对，2026-08-09）**：
      pgr 内建 runlist 区间操作替代外部 spanr；`rept/trf.rs` 已实现"带点
      contig 名映射 → span 处理后恢复"的名字映射规避（原待决策选项③）。
- [ ] **repeat masking：pgi 参数标定 + 真核验证**：CLI 透传已实现
      （`align pgi` 的 `-f/--min-shared/-k/--smer/--window`），但默认值
      未按 §2.5 调整（`--freq` 10 → 100、`--min-shared` 12 → 16 待验证）；
      真核（拟南芥/玉米等转座子丰富）与 RepeatMasker masked 输出对比
      recall（E. coli 无转座子无参考价值）；polyA/卫星低复杂度缺口由
      `rept trf` 兜底（来源：`design/repeat-masking.md` §2.4/§2.5）。
- [x] ~~**wide 128-bit 化：linalg / poa / hv**~~ → **已完成（2026-08-09）**：
      全部改 128-bit（`f32x4`/`i32x4`/`u32x4`，SSE2/NEON 原生，编译 avx2
      无关）。**linalg 需双累加器**（8 元素块拆两个 128-bit，否则单累加器
      依赖链慢 2×：norm 0.72→1.45µs，双累加恢复 747ns）；**poa** 的 avx2
      模块遮蔽 `LANES=8`、`build_root/profile` 泛型化、avx2 自算 n_vec
      （基准 12.1×/8.1× 不变）；**hv** 处 2 重写为双 `i32x4`（RNG 调用数
      不变，AVX2 主路径不变）。三种编译配置（默认/+avx2/aarch64）全验证。
- [x] ~~**SIMD 第二梯队 `paf::cigar`（来源 `design/simd-optimization.md` §6）**~~ →
      **已实现（2026-08-09）**：`classify_alignment` SIMD 分类掩码一次扫描 +
      `scan_cigar_ops`/`scan_cs` 共享 + 位运算 run 跳扫；`maf to-paf` 40 M
      列 0.55 s → 0.347 s（~37%），输出逐字节一致。`rev_comp`/`complement`
      已评估暂缓（`fa rc` 实测仅 ~14%，memset/IO 主导）。
      `count_bases`（第一梯队）已实现（wide 7.5× / AVX2 47×）。
- [x] ~~**hv_benchmark 拆分**~~ → **已完成（2026-08-09）**：现役核心
      （`hash_hv_bit`/`i8`/`sparse`，16 组）留在 `benches/hv_benchmark.rs`
      （全跑 ~2.5 min）；历史对照（AVX-512 ref、RNG 候选、i16/pshufb、
      encode 变体、哈希吞吐）移入 `benches/hv_benchmark_ref.rs`（按需
      filter 跑）。原 54 组全跑 10–20 min 的根因是组数多，非单组数据量。
- [ ] **paf 查询层扩展（待实现）**：`--min-tree-coverage`（Caf Tree
      Coverage 过滤维度，查询时无法全图计算，作传递闭包后处理过滤）；
      `--end-trim` 推迟（需 per-interval 修剪 CIGAR，待序列输出引入时
      一并处理）（来源：`paf-pangenome.md` §Caf 过滤维度对照表）。
- [x] ~~**用户文档改动清单落地（#22）**~~ → **已完成（2026-08-09）**：
      dist.md（frac 无偏/ANI 推荐、hv 粗分层、mini 排序用）此前已落地；
      pbit.md（强制 PAF、空 PAF 禁用 CIGAR、to-paf 命令、MAF 管道、边际
      delta、无 cg:Z 存行）、align-pgi.md（链粒度 ~1 kb vs pbit 4 kb 段，
      建议 chainnet 链路）、pgi.md（近缘距离弱 ρ≈−0.71）本轮补齐
      （来源：`genome-nn-query.md` §8.6）。
- [ ] **chain 算法待验证（低优先）**：KD-tree 已实现并用于 `psl chain`
      （`libs/ds/kdtree.rs`）；`best_crossover` 已接入 `fas_multiz` merge
      （`libs/ds/crossover.rs`）——两者的**真实数据验证**待做；KD-tree
      用于 PAF 链式化 / POA 排序仍待评估（PAF 当前未明确需要链式化）
      （来源：`chain-algorithms.md` §12.3）。
- [x] **FastK/Profex 原生迁移（已实现 2026-08-09）**：`rept s-kmer`/`e-kmer`
      已原生化（`libs/kmer/` 计数 + profile + run 提取，`--keep-index` 缓存
      升级为单文件 `.pgrk`；不做 super-mer/磁盘分桶/外部格式兼容），
      设计 = `design/kmer.md`。
