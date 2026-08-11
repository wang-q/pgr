# `pgr fq assemble`：tadpole contigMode 迁移（设计）

> 2026-08-11。目标：替代 anchr 中 tadpole 的**组装用途**（contigMode）：
> 2_insert_size 流程（硬依赖）与 unitigs 流程（`--unitigger tadpole` 可选
> 分支）。ecc/extend 已由 `fq ecc`/`fq extend` 覆盖（见
> [anchr-merge-replace.md](anchr-merge-replace.md) §6-7）。
> 参考：BBTools-40.01 `assemble/Tadpole*.java`。

## 1. anchr 中的调用点（tadpole 组装）

| 流程 | 调用 | 参数 |
|---|---|---|
| `2_insert_size.tera.sh` | `tadpole.sh in=R1 in2=R2 out=PREFIX.tadpole.contig.fasta threads=N overwrite [prefilter]` | 默认 k=31 |
| `unitigs.tera.sh` | `tadpole.sh in=pe.cor.fa out=unitigs_K{k}.fasta threads=N k={k} overwrite` | k ∈ opt.kmer（如 "31 81"） |

两者都不传 `mode=` → **默认 contigMode**（`Tadpole.java`：无 ecc/extend/
toss 标志时 `processingMode=contigMode`）。下游：2_insert_size 接 bbmap +
reformat-ihist（另计）；unitigs 接 `anchr contained/orient/merge`（anchr
自有，不迁移）。**prefilter 默认关**（anchr `opt.prefilter=0`），本命令
不实现（同 ecc/extend）。

## 2. contigMode 语义（源码确认，逐条移植）

### 2.1 建表与输入

- 与 ecc/extend 相同的 canonical kmer 计数（`TadpoleTable`），minprob
  质量门控（FASTA 无质量 → 不过滤）。
- `kmerRangeMin=0`/`kmerRangeMax=MAX` 默认不做 kmer 范围过滤；
  `removeBubbles=false`/`removeDeadEnds=false`（shave/rinse 默认关，
  跳过）；`processContigs=false`（不建图/不 pop bubble）。

### 2.2 多轮种子（contigPasses=16，contigPassMult=1.7）

`BuildThread.run`（Tadpole2，k>31）：

```
for i = 15 .. 1:
    minCountSeedCurrent = max(3+i, floor(3 * 1.7^i * 0.92 - 0.25))
    扫描全表：count >= 阈值 且 未被 claim 的 kmer → processKmer
最终轮：minCountSeedCurrent = 3，再扫一遍
```

Tadpole1（k≤31）同构（`Tadpole1.BuildThread`）。**单线程下就是 16 轮
全表扫描**，每轮对未认领的种子 kmer 建 contig。

### 2.3 认领（ownership，contigMode 下 useOwnership=true）

- 每个 kmer 有 owner（-1 未认领 / 0..N-1 线程 id）。单线程 id=0：
  认领集合 = HashSet<Kmer>。
- `processCell`：count < 阈值跳过；已认领跳过；认领后 `processKmer`。
- 行走中每个新 kmer 也要认领；`owner==id`（本线程已认领）→ 环形检测，
  返回 `fbranch ? F_BRANCH : LOOP`；被其它线程认领 → BAD_OWNER（单线程
  不会发生）。
- `leftCounts` 在 BuildThread 里**非空**（区别于 extend 模式），因此
  行走时启用左 junction 与隐藏分支检查（`leftMaxPos != evicted` → 停）。

### 2.4 行走（extendToRight，contigMode 版）

入口：count(minCountSeed)=3、owner 检查 → 左/右计数（4 桶）→
`rightMax < minCountExtend(2)` → DEAD_END；`isJunction(rightMax,
rightSecond)` → `isJunction(leftMax,leftSecond) ? D_BRANCH : F_BRANCH`；
`isJunction(leftMax,leftSecond)` → B_BRANCH。循环：取 rightMaxPos 碱基
追加 → 新 kmer 认领/环形检测 → 若 `bbranch` 返回 `fbranch?D_BRANCH:
B_BRANCH`；`hbranch`（leftMaxPos!=evicted 且 branchMult1>0）同；追加后
`fbranch` → F_BRANCH；`rightMax<2` → DEAD_END。

isJunction（深度比阈值，与 ecc/extend 共用）：`second<1 || second*20<
max || (second<=3 && max>=max(2, second*3))` → false（非 junction）。

makeContig（Tadpole2）：种子 kmer → extendToRight → `reverseComplement`
再 extendToRight（即向左延伸）→ `doubleClaim`（单线程恒真）→
`trimEnds(0)` → 长度 `>= k+minExtension(2)` 且 `>= minContigLen` →
生成 Contig（`leftCode/rightCode/leftRatio/rightRatio`，canonical 方向）。
Tadpole1（k≤31）路径同构，junction 方向判定按
`kmer>rkmer`（Tadpole1）/`kmer<rkmer`（Tadpole2）翻转——与现有
`extend_to_right2` 的 `canonical_is_rc = k > 31` 一致。

### 2.5 覆盖度与输出

- `calcCoverage`：contig 每个 kmer（canonical count）均值
  `sum/(float)kmers`，并记 `minCov`/`maxCov`；`coverage<minCoverage(1)`
  或 `>maxCoverage` → 丢弃该 contig（默认不丢）。
- `minContigLen` 默认 `max(124, 2k)`（k=31→124，k=81→162）；`minExtension=2`。
- 输出 FASTA：`>contig_{id},len={len},cov={cov:.1},gc={gc:.3},min={min},
  max={max},hh={hh:.3},caga={caga:.3},left={code},right={code}` + 序列。
  gc/hh/caga 由 `calcScalarsFast` 计算（实现时逐行对照）。
- 排序：length 降序 → coverage 降序 → 序列字典序 → id。
- **顺序无关性**：行走受 leftCounts 隐藏分支约束，contig 集合由
  （图 + 阈值 + 认领）唯一确定，与种子扫描顺序无关 → pgr 的
  HashMap 迭代顺序不影响输出（黑盒验证确认后定案；排序含序列字典序
  兜底，id 只在完全重复时生效）。

## 3. 命令形状

```
pgr fq assemble [OPTIONS] <infiles>...
  -k, --kmer <int>        默认 31（2_insert_size 用；unitigs 由模板按 k 循环调用）
  -o, --outfile <file>    输出 FASTA（默认 stdout）
  -p, --parallel <int|auto>  兼容参数，校验但不启用（确定性单线程，同 ecc/extend）
  --min-contig-len <int>  默认 auto = max(124, 2k)
```

- 单 k 每次调用；unitigs 模板按 k 循环（与现模板逐 k 调用 tadpole 一致）。
- 单线程确定性 = 与 `tadpole.sh threads=1` 黑盒对照的前提。

## 4. 验证

- Lambda `tests/bbtools/Lambda/pe.cor.fa.gz`（或 ecco 输出）k=31/81，
  `tadpole.sh threads=1` 生成 golden（contig FASTA），逐字节对照。
- 关键正确性锚点：contig 集合 + 序列 + 头字段（len/cov/gc/min/max/hh/
  caga/left/right）+ 排序。
- 已知风险：HashMap 迭代顺序与 Java 哈希表不同——若黑盒出现顺序相关
  差异，需复刻 Java 表布局或改确定性扫描顺序（先以实验定论）。

## 5. 不做

- prefilter（anchr 默认 0）；shave/rinse/pop/bubble；`mode=insert`；
  多线程；`trimCircular`/`trimEnds`（默认 0）；kmer 范围过滤。
- 2_insert_size 的 bbmap + reformat-ihist 仍属独立缺口（todo 挂账）。

## 6. 实现状态与已知偏差（2026-08-11 定案）

`pgr fq assemble` 已实现：contig 构建（多轮种子/行走/认领）+ contig 图 +
BubblePopper + 排序重编号 + 输出，全部确定性与单线程等价。

### 已验证（逐字节）

- **pre-pop contig 集合 89/89 与 `tadpole.sh threads=1` 逐字节一致**
  （含短 contig，`popbubbles=f mincontiglen=1` 对照）。
- k=31 左端边行走在 rc 空间（Tadpole1 `processContigLeft` 交换
  kmer/rkmer 语义）——修复了"行走绕回自身生成自环边"的 bug。

### 已知偏差（bubble 解析顺序）

tadpole 的气泡消除**顺序相关**：expand 顺序决定重叠气泡中"谁吸收谁"，
进而决定链式合并。实验证据：

- 把 Java `popBubbles` 的 expand 迭代倒序，输出从 67 变 66 contig；
- Java 的哈希表 cell 顺序随 `-Xmx` 变化（Xmx1g/3g/8g → prime
  228983/213973/194057），但构建顺序 68/89、71/89 个 id 不同——即
  "逐字节一致"本身只对特定内存参数成立；
- pgr 用确定性扫描顺序（canonical kmer 排序）代替 Java 的哈希 cell
  顺序，输出确定且跨运行稳定，但 bubble 解析结果与 tadpole 有少量
  差异（Lambda 2000 对：67 vs 66 contig，总碱基差 ≤100，序列集合
  重合 ≥90%）。

**决策（用户确认，2026-08-11）**：不做哈希表布局复刻，接受确定性输出
+ 文档化偏差。理由：逐字节一致需复刻 BBTools 内存模型（`-Xmx` 相关
prime）+ 开放寻址插入顺序 + 溢出树，约几百行无生物学价值的"镜像
Java 内存布局"代码，且结果脆弱（换 `-Xmx` 即失效）。pgr 的 contig
集合/总碱基与 tadpole 一致，差异只是少数气泡的"走哪条路径"（两条
都是合法组装选择，序列质量等价），对 anchr 用途（insert-size 参考、
unitigs 组装）影响可忽略。

回归测试：`tests/cli_fq_assemble.rs` + golden
`tests/bbtools/Lambda/golden/tadpole_contigs31.fasta.gz`（tadpole
默认输出 67 contig），断言确定性、总碱基差 ≤100、序列集合重合 ≥90%。
