# 替换 anchr merge 流水线中的 BBTools（分析 + 迁移计划）

> 2026-08-11 整理。anchr 的 `templates/merge.tera.sh` 是 read 纠错 + 合并
> 流水线（源自 BBTools `assemblyPipeline.sh`），与 trim 流程不同：涉及
> `bbmerge`/`tadpole` 等新工具，pgr 目前没有对应命令。本文拆解流程、
> 对照 pgr 现状、分析 BBTools 源码，给出迁移计划。
>
> 参考：`merge.tera.sh`（anchr 项目）+ BBTools-40.01 源码。

## 0. anchr BBTools 依赖全景（信息补齐，2026-08-11）

merge 只是 anchr 若干流程之一。全模板扫描（`templates/*.tera.sh`）后，
BBTools 工具依赖如下：

| 流程 | BBTools 工具 | pgr 状态 |
|---|---|---|
| `trim` | bbduk / bbnorm / clumpify / kmercountexact / reformat / repair | ✅ 全迁移 |
| `merge` | bbduk / clumpify / repair ✅ + **bbmerge / tadpole** ❌ | 部分（本文） |
| `2_insert_size` | **tadpole**（组装 contig）/ **bbmap**（比对）/ reformat（ihist） | ❌ 三个新 |
| `unitigs` | **tadpole**（可选组装器，与 bcalm/bifrost 并列） | ❌ 新 |
| `anchors` | **bbwrap/bbmap**（reads 比对到参考） | ❌ 新 |
| `8_spades` / `8_mr_spades` | repair | ✅ `fq split --repair` |

**结论**：anchr 的 BBTools 迁移面有三个新领域——

1. **tadpole**：用途最广（merge 的 ecc/extend、2_insert_size 的组装、
   unitigs 的组装器），本质是 **kmer 图组装器**（contig/unitigs 级），
   merge 里的 ecc/extend 只是它能力的一部分；
2. **bbmap/bbwrap**：短读比对工具（anchors、2_insert_size 用）——
   **pgr 完全没有比对工具**，这是全新领域（与 paf 的隐式图不同，
   bbmap 是参考序列比对）；
3. **bbmerge**：双端 overlap 合并（merge 流程核心）。

因此 merge 迁移（本文）是更大拼图的一部分；tadpole 的迁移价值
**远超 merge 流程本身**（三个流程共用），bbmap 则是另一个独立大项。

### 0.1 BBTools 代码量拆解（50 万行从哪里来，2026-08-11）

`tokei --no-ignore`（BBTools-40.01）：

* 全仓库 Code ≈ **501,813 行**：Java 469,880 + Shell 29,355 + C 2,300；
* Java 总行数 728,209，其中注释 165,927 行（占 23%）；
* 仓库还提交了 2,027 个 `.class` 编译产物（17.6 MB），不算 LOC 但目录显重。

主要构成（Java Code 行）：jgi 72,297（JGI 遗留：Dedupe/BBMerge/Seal/
RQCFilter）、align2 42,756（bbmap 比对核心）、cardinality 41,807
（177 个文件，LogLog 位宽/压缩变体）、stream 30,320（Read/SamLine 等 IO）、
bbduk 24,954（19 个 BBDuk* 文件）、sketch 13,318、var2 13,013、
idaligner 12,156、bin 12,152、aligner 11,571、bloom 10,535、tax 9,678、
ml 9,406、map 9,048、clade 8,996、assemble 8,971（tadpole）。

膨胀原因：

1. **300+ 命令的家族**：369 个 `.sh` 包装 ≈ 305 个入口类，每个都是完整
   独立 CLI（参数解析 + 帮助文本 + 校验），`driver/` 包 7,259 行；
2. **历史累积不删旧代码**：BBDuk 19 个版本并存（`jgi/BBDuk` →
   `bbduk/BBDuk2..6`/`BBDukS`/`BBDukProcessor*`），Dedupe/Dedupe2 同理；
3. **模板/变体膨胀**：cardinality 大量同算法不同位宽/压缩参数的变体类、
   `SIMDByte256`/`SIMDByte256_only` 双份、`template/` 代码生成；
4. **与 pgr 无关的领域**：var2（变异）、tax、ml、sketch、idaligner、
   clade、synth、barcode、hiseq、covid、server 等合计十几万行；
5. **Java 风格**：一类一文件、超大工具类（`Tools.java` 190 KB、
   `Read.java` 128 KB）、注释占 23%。

结论（与本文关系）：merge 流程真正相关的源码量很小——BBMerge 3,340 +
BBMergeOverlapper 1,503 + Tadpole 7,563 + Clumpify ~700 + bbduk qtrim
相关 ≈ 2 万行内，且多为参数枚举与边界处理。**决定不碰 bbmap 后**，
align2 42,756 + aligner 11,571 + map 9,048 ≈ 6.3 万行整体划出范围，
剩下与 anchr 相关的核心就是 tadpole + bbmerge。pgr 用 9 万行 Rust 覆盖
同样功能面是合理的：50 万里大部分是"壳 + 历史版本 + 无关领域 + 变体"。

## 1. 流程拆解（merge.tera.sh）

anchr `merge` 命令生成 `merge.sh`，输入 R1/R2（可选 SE），流程如下：

| # | 步骤 | 工具 | 作用 | pgr 现状 |
|---|---|---|---|---|
| 1 | clumpify | `clumpify.sh dedupe dupesubs=0` | 排序 + 精确整对去重（SE 追加） | `fq clump --dedupe` ✅ |
| 2a | EC phase 1 | `bbmerge.sh ecco mix vstrict` + `ihist` | overlap 区域纠错（不合并） | ❌ 新 |
| 2b | EC phase 2 | `clumpify.sh passes=4 ecc unpair repair` | clump 内共识纠错（单端） | ❌ 新（**可跳过**，见 §3.3） |
| 2c | EC phase 3 | `tadpole.sh ecc tossjunk tossdepth=2 tossuncorrectable` | kmer 图纠错 + 丢弃坏 read | ❌ 新 |
| 3 | Read extension | `tadpole.sh mode=extend el=20 er=20 k=62` | 3'/5' 端扩展 reads | ❌ 新 |
| 4 | Read merging | `bbmerge-auto.sh strict k=81 extend2=80` + `ihist` | overlap 合并 + insert size 直方图 | ❌ 新 |
| 5 | Dedupe merged | `clumpify.sh dedupe dupesubs=0` | 合并后 reads 去重 | `fq clump --dedupe` ✅ |
| 6 | bbduk qtrim | `bbduk.sh qtrim=r trimq=... minlen=...` | 未合并 reads 质量修剪 | `fq clean`（无 ref）✅ |
| 7 | repair | `repair.sh repair` | 拆分 R1/R2/singles | `fq split --repair` ✅ |

> 步骤 2a/2b/2c 由 `opt.ecphase` 控制（anchr 默认 `"1 2 3"` 全开）；
> `opt.prefilter` 默认 0（不开，内存够时避免 countmin 两次耗时）。
> **anchr `template.rs` 的 `--ecphase` 帮助明确 "Phase 2 can be skipped"**，
> 且用户反馈 phase 2 以前经常卡住（见 §3.3）。

## 2. 工具对照

### 已有（无需迁移）

* `fq clump --dedupe`（步骤 1/5）：`dupesubs=0` 精确整对去重，已 golden 对照；
* `fq clean`（步骤 6）：无 ref 纯 qtrim（merge 用 `qtrim=r trimq minlen`，
  与 trim.era.sh 的纯 qtrim 一致，已解决）；
* `fq split --repair`（步骤 7）：repair.sh rp 模式按名字配对，已对照。

### 需迁移（新命令/功能）

| 工具 | 功能点 | 规模评估 |
|---|---|---|
| **bbmerge** | overlap merge + ecco + ihist | 核心算法 `mateByOverlapRatio` 已有简化版（trim_adapter 的 tbo）；merge/ecco/ihist 需完整移植，~3000 行源码 |
| **tadpole** | kmer 图纠错（ecc）+ read 扩展（extend） | 最大工程：基于 kmer 图的组装式算法，~8000 行源码；pgr 有精确 KmerTable（`libs/kmer`）可作基础 |
| **clumpify ecc** | `passes=4 ecc unpair repair` 模式 | **建议不做**（与 tadpole ecc 冗余，anchr 可跳过，见 §3.3） |

## 3. 新工具源码分析（BBTools-40.01）

### 3.1 BBMerge（`jgi/BBMerge.java` 3340 行 + `BBMergeOverlapper.java`）

* **overlap 检测**：`mateByOverlapRatio`（ratio 打分：找双端重叠最优
  insert size；`mateByOverlapJava_WithQualities`/`mateByOverlapJava` 纯
  Java 路径，JNI 在两版本中禁用）。**pgr trim_adapter 的
  `mate_by_overlap_ratio` 已移植无质量路径（tbo 用）**——merge 需要
  质量路径（`_WithQualities`）与完整参数（minoverlap/mininsert/maxratio/
  margin/offset）。
* **ecco 模式**（`ecco=t`）：overlap 检测后只纠错重叠区域、不合并输出
  （`mix=t` 时合并+未合并同文件输出）。纠错 = 重叠区碱基按质量/一致取
  共识。
* **merge 模式**（bbmerge-auto 的 `strict k=81 extend2=80`）：overlap 后
  合并成单条 read（含 kmer 校验 `kfilter`、`extend2` 失败后扩展）。
* **ihist**：insert size 直方图（overlap 检测的 insert 分布）——顺带解决
  todo 里 2_insert_size 的 ihist 缺口（reformat ihist 是另一个实现）。
* 严格度：`strict`/`verystrict`/`vstrict`（降低 FP）+ `loose` 族。
* `prefilter`：countmin sketch 预过滤低深度 kmer（防内存爆炸）。

### 3.2 Tadpole（`assemble/Tadpole.java` + `Tadpole1/2` + `TadpoleWrapper`，~8000 行）

* **kmer 图**：建 kmer 表（默认 31，extend 用 k=62）→ 每个 read 沿图
  扩展/纠错（类似 SPAdes/velvet 的 de Bruijn 图路径）。
* **ecc 模式**（`mode=correct`/`ecc`）：逐 read 走图，低覆盖路径判错修正；
  纠错策略多（`ecctail`/`eccpincer`/`eccreassemble`/`eccrollback`/
  `requirebidirectional`），anchr 用默认。
* **toss 参数**：`tossjunk`（无效字符）、`tossdepth=2`（低深度 kmer 比例
  超限丢弃）、`tossuncorrectable`（无法纠错丢弃）——"丢弃坏 read"语义
  与 `s-filter` 一致（可复用 QualityTable/判定逻辑）。
* **extend 模式**（`mode=extend el=20 er=20`）：从 read 端沿 kmer 图延伸
  el/er 个碱基（3'/5'）。
* **insert 模式**（`mode=insert` + `ihist`）：tadpole 也能检测 insert
  size（`findOverlap(r1,r2)` 找双端重叠算分布）——但 merge.tera.sh 的
  ihist 用的是 **bbmerge 的 ihist**（overlap 合并时统计），tadpole
  insertMode 不在流程必需范围（可选，做完整 tadpole 时再考虑）。
* `prefilter`：同 bbmerge（countmin sketch）。

### 3.3 clumpify ecc（`clumpify.sh passes=4 ecc unpair repair`）

* ecc：clump 内共识纠错（排序后共享 kmer 的 reads 组内纠错，`Clumpify.java`
  把 `ecco` 转发给内部处理）；`unpair`（单端化）；`repair`（修复配对）；
  `passes=4`（多轮）。pgr `fq clump` 已实现排序/去重，缺这三个模式。
* **可跳过（建议不做）**：三个纠错阶段的分工是——phase 1（bbmerge ecco）
  负责配对 overlap 特有纠错（tadpole 单端做不了）；phase 3（tadpole ecc）
  负责单端 kmer 图纠错且**更强**（专门的组装式纠错器，多策略）；phase 2
  （clumpify ecc）与 phase 3 功能重叠但更弱（clump 共识 vs kmer 图），
  anchr 文档明示可跳过，且 `unpair + passes=4` 多轮建图是内存/时间
  卡点（用户反馈）。**跳过 phase 2 不影响 phase 1/3 的纠错效果**。

## 4. 迁移策略（建议分阶段）

### 阶段 A：bbmerge 核心（merge + ecco + ihist）

* 新命令（命名待定，候选 `fq merge`）：overlap merge（质量路径
  `mateByOverlapRatio_WithQualities`）+ `--ecco` + `--ihist`；
* 复用 trim_adapter 的 `mate_by_overlap_ratio`（补质量路径）；
* 验证：Lambda 双端 reads 与 `bbmerge.sh` 黑盒对照（merge 率/insert 分布）。

### 阶段 B：clumpify ecc + tadpole ecc

* **clumpify ecc（phase 2）建议不做**（与 phase 3 冗余、anchr 可跳过、
  用户反馈卡住）——除非要逐字节对齐默认 `"1 2 3"` 全流程；
* 新命令 `fq ecc`（或并入）：tadpole 式 kmer 图纠错 + toss 参数
  （tossjunk/tossdepth/tossuncorrectable）——可复用 `s-filter` 的
  QualityTable 与判定语义；
* 验证：与 `clumpify.sh ecc`/`tadpole.sh ecc` 黑盒对照。

### 阶段 C：tadpole extend + prefilter

* `fq extend`：read 端扩展（kmer 图沿伸）；
* `prefilter`（countmin sketch）：内存保护，大数据才需要；
* 验证：合成 reads 扩展长度对照。

## 5. 风险与决策点

1. **迁移范围**：tadpole 是 ~8000 行组装式算法，完整移植工作量大。建议
   先做 bbmerge（阶段 A，复用已有 overlap 逻辑）；tadpole 的 ecc/extend
   按"丢弃坏 read"优先（toss 语义与 s-filter 一致）。
   **clumpify ecc（phase 2）不做**（与 phase 3 冗余 + anchr 可跳过 +
   用户反馈卡住），单端纠错由 tadpole ecc 承担。
2. **命名**：merge 流程的核心命令候选 `fq merge`（bbmerge 等价）、
   `fq ecc` / `fq extend`（tadpole 等价）；沿用 pgr 长名风格 + bbduk 标注。
3. **ihist 归属**：insert 检测在 BBTools 有三处——bbmerge ihist（merge
   流程用，随阶段 A 迁移）、tadpole `mode=insert`（不在流程、可选）、
   reformat ihist（2_insert_size 用，todo 挂账）。merge 流程的 ihist 由
   bbmerge 满足；reformat ihist 仍是独立缺口。
4. **验证**：所有迁移命令与 BBTools 39.38/40.01 黑盒对照（Lambda 双端
   reads），沿用 trim 流程的 golden 思路（不新增测试料）。
5. **prefilter 优先级**：默认关（anchr prefilter=0），阶段 C 再考虑。

---

*参考来源: [merge.tera.sh](../../../anchr/templates/merge.tera.sh)（anchr 项目，只读） |
BBTools-40.01 源码（`jgi/BBMerge*`、`assemble/Tadpole*`、`jgi/Clumpify*`）*
