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
| `merge` | bbduk / clumpify / repair ✅ + **bbmerge / tadpole** ✅ | 全迁移（本文） |
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

### 0.2 bcalm / Bifrost（unitigs 流程的 unitigger，2026-08-11）

anchr `unitigs.tera.sh` 支持 4 种 unitigger：tadpole / bcalm / bifrost /
superreads（MaSuRCA 的 `create_k_unitigs_large_k`），`--unitigger` 选择；
**template 默认 bcalm**，独立 `anchr unitigs` 命令默认 superreads（两处
默认不一致，属 anchr 侧问题）。输入均为 merge 流产物 pe.cor.fa。

实际调用（unitigs.tera.sh）：

* bcalm：
  `bcalm -in pe.cor.fa -kmer-size K -abundance-min 3 -verbose 0 \
   -nb-cores N -out K{K}`，产物 `K{K}.unitigs.fa` → `unitigs_K{K}.fasta`；
* Bifrost：
  `Bifrost build --input-seq-file pe.cor.fa --kmer-length K \
   --clip-tips --del-isolated --threads N --fasta --no-compress-out \
   --output-file unitigs_K{K}`；
* 产出 unitigs 后接 `anchr contained/orient/merge`（anchr 自有命令，
  无需迁移）。

**两者都不是 BBTools**，anchr 直接调用外部工具（check_dep 必装、
install_dep 安装；本机已装 bcalm + Bifrost 1.3.5）→ **unitigs 流程本身
没有 BBTools 依赖**，唯一 BBTools 分支是 `--unitigger tadpole`（可选）。
pgr 迁移面不含 unitig 组装；merge 流程的 tadpole 只取其 ecc/extend。

**算法对照**：

* **BCALM2**（Chikhi et al., ISMB 2016）：精确 cDBG。先建 solid kmer 集合
  （`-abundance-min` 阈值），再压实成 unitig（= 最大无分支路径，确定性）；
  用 minimizer 把 kmer 分区到不同 worker，低内存并行，输出 FASTA 头带
  `LN`（unitig 长度）/`KC`（kmer 深度）/`km` 和左右边信息。本机版本
  BCALM2 1.3.1 支持 minimizer/bloom（`neighbor` + cascading debloom）/
  branching 节点存储等高级选项；kmer 转 canonical。
* **Bifrost**（Holley & Melsted, 2020）：直接构建压缩 de Bruijn 图，
  最小完美哈希 + Bloom filter（默认 24 bits/kmer），支持彩色图（可增量
  update/query）。`--clip-tips` 剪长度 < k 的 tip，`--del-isolated` 删
  长度 < k 的孤立 unitig；kmer 唯一出现 1 次默认直接丢弃（与 bcalm 的
  abundance 阈值思路类似）。
* **Tadpole** 是贪心近似：不建显式图，kmer 表上沿 path 延伸，branch/
  dead-end 停，深度比可穿强分支；作为 unitigger 结果不确定（ecc/extend
  场景没问题，unitig 组装不如前两者）。

**结论**：若将来 pgr 要接管 unitigs（`pgr unitigs`），参考实现应是
BCALM 式 minimizer 分区 + cDBG 压实（确定性 unitig，可并行、内存可控），
而非 Tadpole 克隆；彩色图（Bifrost 的差异化能力）暂不需要。本次 merge
迁移不动 unitigs，bcalm/Bifrost 维持外部依赖即可。

## 1. 流程拆解（merge.tera.sh）

anchr `merge` 命令生成 `merge.sh`，输入 R1/R2（可选 SE），流程如下：

| # | 步骤 | 工具 | 作用 | pgr 现状 |
|---|---|---|---|---|
| 1 | clumpify | `clumpify.sh dedupe dupesubs=0` | 排序 + 精确整对去重（SE 追加） | `fq clump --dedupe` ✅ |
| 2a | EC phase 1 | `bbmerge.sh ecco mix vstrict` + `ihist` | overlap 区域纠错（不合并） | `fq ec-overlap` ✅ |
| 2b | EC phase 2 | `clumpify.sh passes=4 ecc unpair repair` | clump 内共识纠错（单端） | 跳过（见 §3.3） |
| 2c | EC phase 3 | `tadpole.sh ecc tossjunk tossdepth=2 tossuncorrectable` | kmer 图纠错 + 丢弃坏 read | `fq ec-kmer` ✅ |
| 3 | Read extension | `tadpole.sh mode=extend el=20 er=20 k=62` | 3'/5' 端扩展 reads | `fq extend` ✅ |
| 4 | Read merging | `bbmerge-auto.sh strict k=81 extend2=80` + `ihist` | overlap 合并 + insert size 直方图 | `fq merge --extend2 --rem` ✅ |
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

### 已迁移（2026-08-11 全部完成，逐字节 golden 对照）

| 工具 | pgr 命令 | 对照 |
|---|---|---|
| **bbmerge** | `fq merge`（--strict / --no-make-vector / --extend2 --rem） | `merge.*` + `merge4.*` golden，net/classic/ecco/ihist 全一致 |
| **tadpole ecc** | `fq ec-kmer`（tossjunk/tossdepth/tossuncorrectable） | `ecct_sub.fq.gz` golden，丢弃判定一致 |
| **tadpole extend** | `fq extend`（k=62 el/er） | `ext_sub.fq.gz` golden，扩展碱基一致 |
| **clumpify ecc** | 跳过 | 与 tadpole ecc 冗余，anchr 可跳过（§3.3） |

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

**算法复杂度（原作者自述）**：Brian Bushnell 2015 年发布帖（SEQanswers
"Introducing Tadpole"）明说 Tadpole "is only a contig-builder"——把 kmer
组装成 contig，遇到 branch 或 dead-end 就截断；**不建显式 de Bruijn 图、
不消除杂合 bubble、不做完美遍历、不做 scaffolding**（"It does not generate
the explicit DeBruijn graph and try to remove heterozygous bubbles, or find
a perfect traversal"）。设计初衷是给 BBMerge 做 read 延伸/纠错（"my primary
design goals were for read extension and error-correction"）。官方 guide 同：
"does not do any complicated graph analysis"。纠错方式 = "组装穿过错误"——
用图路径穿过错误处、以组装出的碱基替换错误碱基。

实际算法就四块：kmer 计数表 → 贪心延伸（查 4 个后继深度，branch 判定用
**深度比自适应阈值** `branchmult1=20`/`branchmult2=3`/`branchlower=3`，非绝对
深度）→ 纠错（`reassemble` 默认 + `pincer`/`tail` 可选）→ toss（低深度 kmer
占比超限丢弃）。shave/rinse/pop 是可选增强、默认关，为组装连续性服务，
merge 流程不涉及。

**实现策略（建议）**：算法照原版（保守 + 深度比分支 + 组装式纠错），代码不
照抄——砍掉 Tadpole1/Tadpole2 双实现（kmer 编码统一即可）、几百个选项
（流程只用 `k`/`el`/`er`/`ecc`/三个 toss 参数 + 分支判定默认值）、线程/IO
管道（rayon + pgr fq reader）、shave/rinse/pop/prefilter。核心逻辑估
~2-3 千行，大头是与原版黑盒对照（Lambda golden）。需保留行为参数：bm1/bm2/
blc 默认值、ecc 默认策略（reassemble + rollback + rbi）、toss 语义
（`tossdepth=2` 时 pair 任一 read 失败即丢）。

**同类思路：MaSuRCA super-reads（Zimin et al. 2013, Bioinformatics）**：
super-read 把短 reads 用 kmer 查找表在 5'/3' 两端、**延伸唯一时**逐碱基延成
更长的伪长读：先纠错（QuORUM）简化图，再沿 k-unitig（无分支最大路径，即
唯一延伸路径）延伸，把 50-100× 覆盖压成 2-3× 喂给 OLC（Celera Assembler）。
与 Tadpole 同源思路（kmer 唯一路径延伸 + 保守），差异：MaSuRCA 先独立纠错
再延伸、严格"遇 branch 即停"；Tadpole 把纠错融入延伸（组装式）、分支处按
深度比自适应决策、不接 OLC。侧面印证"唯一路径延伸"是成熟简单的设计模式，
pgr 简化实现可行。

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
* 新命令 `fq ec-kmer`（或并入）：tadpole 式 kmer 图纠错 + toss 参数
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
   `fq ec-kmer` / `fq extend`（tadpole 等价）；沿用 pgr 长名风格 + bbduk 标注。
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

## 6. 实现状态（bbmerge 迁移完成，2026-08-11）

**`pgr fq merge` 已实现并通过 BBTools 40.01 黑盒逐字节对照**（Lambda 40k
pairs，`tests/cli_fq_merge.rs`）：

* `--ecco --mix --vstrict --net bbmerge.bbnet` ≡ `bbmerge.sh ecco mix vstrict`
  （anchr merge phase 1，ihist 也一致）；
* `--strict --net` ≡ `bbmerge.sh strict`（merged + unmerged + ihist 一致）；
* `--no-make-vector` ≡ `bbmerge.sh ... makevector=f`（经典 efilter/pfilter
  路径，vstrict/strict、ecco/join 全一致）。

关键发现（BBTools 40.01）：

1. **`BBMerge.main()` 无条件置 `MAKE_VECTOR=true`**（除非用 tadpole），把
   ratio 预筛 `maxratio` 强制成 0.7，并跳过 ambig/pfilter 拒绝；最终合并
   与否由 **bbmerge.bbnet 神经网络**（23 维特征，6 层稠密网络，`##ctf`
   阈值）决定。pgr 移植了该 net 的推理（`libs/fq/bbnet.rs`，含 SIG/TANH/
   MSIG/RSLOG 激活与 SIMD.fma 点积语义）；
2. **质量值在解析时转 phred**（`applyQualOffset`，-33；N 碱基置 0，ACGT
   至少 2），输出写回 +33 —— 未合并 reads 的碱基/质量也会被规范化；
3. no-quality 路径的 `bestGood/secondBestGood` **永远为 0**（40.01 未赋值），
   是 net 特征的重要输入；
4. `fq merge` 默认走 net（与 bbmerge.sh 一致），`--net` 必填；classic 路径
   用 `--no-make-vector`。

## 7. 实现状态（tadpole ecc/extend 迁移完成，2026-08-11）

**`pgr fq ec-kmer` / `pgr fq extend` 已实现并通过 BBTools 40.01 黑盒逐字节对照**
（Lambda 40k pairs 全量 + 2k pairs 子集 golden，`tests/cli_fq_ecc.rs` /
`tests/cli_fq_extend.rs`，fmt/clippy clean，全量测试绿）：

* `fq ec-kmer ... --toss-junk --toss-depth 2 --toss-uncorrectable` ≡
  `tadpole.sh ecc tossjunk tossdepth=2 tossuncorrectable`（phase 3，
  丢弃判定 1702 对完全一致）；
* `fq extend -k 62 --el 20 --er 20` ≡ `tadpole.sh mode=extend el=20 er=20
  k=62`（read extension 步骤，k>31 走 Tadpole2 路径，扩展碱基 1,444,344
  完全一致）。

### 7.1 关键发现（Tadpole 源码语义，逐条复刻）

1. **N 会重置 kmer 窗口与 minprob 乘积**：`AminoAcid.baseToNumber` 静态
   初始化先 `Arrays.fill(..., -1)` 再写 ACGT(U)，所以 N 是 -1 而非 0；
   计数表永远不含跨 N 的窗口（表构建、fillKmers、hasKmersAtOrBelow、
   isJunk、reassemble_inner、regenerateCounts 全部一致重置）。
2. **计数数组带符号**：缺失 kmer 的 `getCount` 返回 -1（NOT_PRESENT），
   N 窗口在 fillCounts 里是 0；isError/isSimilar 全部按 Java 有符号语义
   （low=-1 时 `low*em1<high` 恒真等）。
3. **修正判定 `num==rightMax` 是"碱基编码 == 计数"**（不是位置索引）：
   当碱基编码恰好等于深度时 BBTools 视作已修正跳过 —— 移植必须复刻。
4. **extend 模式不做纠错**：`mode=extend` 时 Java `ecc_=false`（correctMode
   才置 true），`processRead` 直接扩展；pgr 用 `TadpoleOptions.ecc` 开关
   区分（`fq ec-kmer` 置 true、`fq extend` 默认 false）。
5. **扩展不做左分支检查**：`ExtendThread.leftCounts` 从未初始化（null），
   `extendToRight2_inner` 只判右 junction；pgr `extend_to_right2` 的
   `use_left=false`。
6. **junction base 追加条件随 Tadpole 版本翻转**：Tadpole1（k≤31）
   `kmer>rkmer`；Tadpole2（k>31，canonical 为 MIN）`kmer.key()==array2()`
   ⇔ `kmer<rkmer`。
7. **扩展 seed 含 N 直接失败**：`rightmostKmer` 对末尾 k 碱基含 N 返回
   null → 该方向不扩展（pgr seed 构建带 N 重置）。
8. **fromRight 的 clearWindow2 在反向方向做**（Java `reverseInPlace` +
   reverse quals），正向方向 + 反向 quals 会把端头修正误清。
9. **similar_range 负 loc2 是空区间 → true**（Java clamp 到 -1 循环不执行），
   Rust usize 转换会包成巨大值，必须显式判负。

### 7.2 merge phase 4（bbmerge-auto extend2/rem，2026-08-11 完成）

**`pgr fq merge --strict --no-make-vector --extend2 80 --rem` 已实现并通过
BBTools 40.01 黑盒逐字节对照**（`ext_sub.fq.gz` 2000 对，merged/unmerged/
ihist 全一致，`merge4.*` golden + `cli_fq_merge.rs` 测试）。

实现要点（逐条对照 Java `BBMerge.processReadPair`）：

1. **触发条件**：`rem`（requireExtensionMatch）时**每个 pair** 都走扩展
   块（条件 `requireExtensionMatch || AMBIG || NO_SOLUTION`），不是只对
   AMBIG/NO_SOLUTION——初始 Merged 的 pair 也会扩展后重检；
2. **扩展调用**：`tadpole.extendToRight2(..., includeJunctionBase=false)`，
   与 `fq extend` 的 `true` 不同；`leftCounts` 仍为 null（
   `extendThroughLeftJunctions` 默认 true），左分支检查关闭；k=81 走
   Tadpole2 路径（多字 kmer，见下）；
3. **迭代次数**：`extendIterations` 默认 1——只扩展一轮（每 read 至多
   extend2 碱基），不是无限循环；
4. **rem 接受规则**：`lengthSum` 用**扩展前** reads 长度（
   `approxMaxOverlappingInsert = lengthSum - 26`）；`minExt =
   min(12, extend2*2)`；只有"扩展前无 overlap 且扩展后 insert 超过
   approxMax 且 extension >= minExt"才接受扩展合并；
5. **unmerged 输出用原始 reads**：BBMerge 在 `findOverlapInThread` 快照
   原始 reads（`originals`），未合并 pair 写 outu 前恢复——扩展后的 reads
   不进 outu；
6. **histogram 只记最终结果**：扩展块结束后重跑 overlap 检测（即使
   e1=e2=0），按最终 outcome 记 stats（BBMerge 每 pair 只计一次）；
7. **多字 kmer 支持（k>64）**：tadpole 的 kmer 从 u128 重构为
   `Vec<u64>` 多字表示（`Kmer`），`push_right` 的最高字掩码按 `2k%64`
   取位（k=62 时 60 位，旧 `!(3<<62)` 会留垃圾位）、`push_left` 的进位
   方向（低字顶部落入高字底部）、junction 方向判定（`is_lt`）三个隐蔽
   bug 修复后，k=62/81 均与 BBTools 逐字节一致；
8. **bbmerge-auto.sh 本身**只是内存自动检测的包装，参数原样传给
   `jgi.BBMerge`，无需单独迁移。

至此 anchr merge 流程 7 步全部有 pgr 等价命令（clumpify 去重 / bbmerge
ecco / clumpify ecc（跳过）/ tadpole ecc / tadpole extend / bbmerge-auto
extend2+rem / clumpify 去重 / bbduk qtrim / repair）。

## 8. 代码审核修复（2026-08-11）

对 bbmerge/tadpole 迁移代码做了一轮完整审核（对照 BBTools-40.01 源码 +
本机 Java 实测），修复如下：

1. **panic 修复**：`--kmer 0` 在 clap 层校验（`RangedU64ValueParser`，
   `libs::fq::tadpole::run` 另有 `k>=1` 防御）；FASTA（无质量）输入在
   `count_errors`/`count_errors_from` 传 `None`（BBTools null quality 用
   固定 q=20），不再越界。
2. **`--ecco` 默认 mix**：对齐 `bbmerge.sh`（`ecco && !setMix` 自动
   `MIX_BAD_AND_GOOD=true`）；新增 `--no-mix` 表达 `mix=f`。此前
   `--ecco` 不带 `--mix` 会静默丢弃未合并 reads（golden 一直显式带
   `--mix`，未暴露）。
3. **efilter 从死选项变为生效**：classic 模式按 Java 语义先跑
   expected-error filter（触发则跳过 pfilter，保护"观测 bad 与质量预期
   一致"的 pair 不被 pfilter 丢弃）；`extraMult` 在 make-vector 模式
   对齐 Java 的 4.0（原恒 1.2，低危保真度偏差）。合成用例与 Java
   `strict makevector=f` 逐字节一致。
4. **清理**：`count_read_kmers` 里残留的 `chain12_hits` cfg(test) 探针、
   `extend_read`/`process_read` 的 `let _ =` 死语句、libs 里不可达的
   `--extend2 requires --no-make-vector` bail（cmd 层与 Java 自身都会
   强制 make_vector=false）均删除。
5. **`--parallel`**：ecc/extend 接受但忽略，改为校验值合法性
   （`parse_parallel_auto`），帮助文本注明"为 tadpole.sh CLI 兼容保留、
   忽略（确定性单线程）"。
6. **文档**：`docs/fq.md` 补上 merge/ecc/extend 三个子命令章节
   （此前全 fq 子命令唯独这三个没有文档）。
7. **测试**：新增 7 个回归测试——ecco 缺省 mix 与 golden 逐字节一致、
   `--no-mix` 只输出可纠正对、efilter 保护合成用例、`--kmer 0` 报错不
   panic、FASTA 输入不 panic、`--parallel` 校验。

**复核更正**：审核初稿把"双 1bp reads"列为 RET_BAD/RET_AMBIG 偏差，
重读 `findOverlapInThread` 写出逻辑后确认 Java 的 RET_BAD 同样写入
outu（`else if(listb!=null)`），pgr 语义一致，无需修复；
`from_phred` 的 saturating 与 Java byte 溢出的差异仅 phred>122 可达，
内部质量上限（merge≤50、tadpole≤32）使该分支实际不可达，保持现状。
