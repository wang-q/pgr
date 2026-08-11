# `pgr asm map`：完美匹配 reads 映射（anchors 流程 bbmap 替代，设计）

> 2026-08-11。目标：替代 anchr `anchors.tera.sh` 中的 `bbwrap.sh`
> `perfectmode maxindel=0 strictmaxindel` 调用。需求 = **exact string
> match**（read 整条完美匹配参考，无错配、无 gap），不是完整比对器。

## 1. 需求来源（anchors.tera.sh，已核对本地 anchr 模板）

```bash
bbwrap.sh maxindel=0 strictmaxindel perfectmode \
    threads=... ambiguous=all nodisk append \
    ref=UT.fasta in=R1,R2 \
    outm=mapped.sam outu=unmapped.sam basecov=basecov.txt
```

下游实际只用：

1. `mapped.sam` / `unmapped.sam`——`wc -l` 算 mapped 比例，随即删除；
2. `basecov.txt`（`RefName Pos Coverage` 三列，Pos 0-based）——按覆盖度
   中位数/MAD 找正常覆盖区间 → spanr fill/excise → anchor 区域 → 从
   UT.fasta 切出锚定序列。

BBTools 源码确认 `perfectmode` 语义（`AbstractMapThread.java:1371`
`maxMismatches=0`；`BBMap.java setPerfectMode` `MINIMUM_ALIGNMENT_SCORE_RATIO=1.0`；
`maxindel=0 strictmaxindel` 禁 indel）→ read 必须与参考某位置整条一致。

**2026-08-11 更新**：basecov 不再由 `asm map` 在内存里累计，改为从
mapped SAM 派生（`pgr sam to-rg` → `pgr rg coverage`）。理由：职责分离
（map 只出 SAM），覆盖度从产物推导与内存累计数学等价（完美匹配 CIGAR
恒为 `L M`，正/反链都只是一个区间，重复区域的 ambiguous=all 每个 hit 一
条记录，两侧语义自动一致），且省掉 map 内 4 B/bp 的 AtomicU32 数组。
模板侧的两个消费者都能用 rg/runlist 组合表达：median/MAD 对 detailed
JSON 按层加权；covered 区域 = `rg coverage -d` → `runlist some`（lower..upper
层名）→ `runlist combine --op union`。

## 2. 算法

### 2.1 参考索引（一次构建）

- 输入：UT.fasta（组装结果，通常 Mb 级）。
- 对每个 k-mer 窗口（默认 k=31，`--kmer`，上限 64 = u128 key）用
  `libs::kmer::canonical_keys` 收集 canonical key + 位置。
- 记录 `(key, contig_id, pos)`，按 key 排序（`libs::ds::radix_sort_u128`
  Myers American-flag），key 相同的候选位置连续 → 二分查找区间。
- 不需要存 strand：验证阶段同时比较 read 与其 rc。

### 2.2 read 匹配（rayon 并行，逐 read 独立）

1. read 长度 < k → 直接 unmapped（完美匹配语义下无解）。
2. 种子 = read 首 k-mer 的 canonical key（完美匹配时必在索引中；
   其 canonical 与 read 实际匹配位置的 canonical 相同，正/反链皆然）。
3. 对索引中该 key 的全部候选位置 `(cid, pos)`（重复区域多候选，
   `ambiguous=all` 语义）：
   - `pos + L <= contig.len()` 且 `contig[pos..pos+L] == read` → forward；
   - 或 `contig[pos..pos+L] == rc(read)`（rc 只算一次）→ reverse。
4. 收集全部通过验证的位置 → mapped（可能多个）；全失败 → unmapped。

### 2.3 流式处理与 SAM 头

- reads 按 100k 条分块读入（有界内存），每块 rayon 并行匹配后按输入序
  写出。
- `--outm`/`--outu` 都写相同的 SAM 头（bbmap 的 outm/outu 共享
  `useSharedHeader`，头对称保证模板 `wc -l` 比例不偏）。

### 2.4 输出

- **SAM**：标准头（`@HD` + 每个 contig 的 `@SQ`，与 bbmap 一致，模板
  `wc -l` 行为接近原版）；mapped 行 FLAG 0/16、POS 1-based、CIGAR
  `L M`、MAPQ 255；unmapped 行 FLAG 4、`*`/0。按输入 read 顺序输出。

### 2.5 覆盖度派生（组合管道，不在 map 内）

```
pgr asm map UT.fasta R1.fq.gz R2.fq.gz --outm mapped.sam --outu unmapped.sam
pgr sam to-rg mapped.sam > mapped.rg        # 每条 mapped 记录一行 chr:start-end
pgr rg coverage mapped.rg -m 2 -o cov.json  # 逐位深度（runlist JSON）
```

- `sam to-rg`：跳过 `@` 头与 unmapped（FLAG 0x4 / RNAME `*` / POS 0）；
  CIGAR 的 M/D/N/=/X 计入跨度，I/S/H/P 不计（本命令只产出 `L M`，
  解析器按通用 SAM 写）。
- `rg coverage -d`（detailed）输出"每深度一层"的 runlist JSON；模板的
  median/MAD 按层加权（run 长度 × 深度），covered 区域 = 层名在
  [lower, upper] 的层 `runlist some` + `runlist combine --op union`。
- 内存画像：map 内零覆盖度内存；rg 侧事件数组约 32 B/read，与参考大小
  无关（与原来的 4 B/bp 恰好相反，锚定场景两者都只是几十~几百 MB）。
- 代价：mapped SAM 必须保留并重读一遍（模板本来就要 cat 一次数行数，
  可合并）。百万级 reads 时是秒级 I/O，可接受。

## 3. 命令形状

```
pgr asm map [OPTIONS] <ref.fa> <reads.fq...>
  -k, --kmer <int>      种子 k-mer 长度，默认 31，上限 64
  --outm <file>         完美匹配 reads 的 SAM（mapped.sam 兼容）
  --outu <file>         未匹配 reads 的 SAM（unmapped.sam 兼容）
  -p, --parallel <int>  真实并行（rayon，与 assemble 的单线程确定性不同）

pgr sam to-rg [OPTIONS] <infile>            # SAM → .rg（stdin 可用）
```

- reads 接受 1+ 个文件（R1/R2 或多个单端文件），FASTQ 或 FASTA。
- 输出按 read 顺序确定（并行分块收集后按序写出），跨运行逐字节一致。
- 覆盖度走 §2.5 的组合管道，map 不再有 `--basecov`。

## 4. 验证

- 合成参考 + reads：
  - 精确正向/反向匹配 → mapped，位置/方向正确；
  - 1 个错配 / 1 个 gap → unmapped（完美语义）；
  - 重复区域 read → 多个 mapped 记录（ambiguous=all）；
  - read 长度 < k → unmapped。
- Lambda 数据 sanity：mapped 比例合理、确定性（两次运行逐字节一致）、
  `sam to-rg` → `rg coverage` 的深度与 mapped SAM 一致（合成测试：两条
  50 bp read 重叠 40 bp → `rg coverage -m 2` 输出恰为 `21-60`）。
- 与 bbmap 黑盒对照暂不做（本机 Java 配对读 gz 失败；语义由 perfectmode
  源码 + 合成测试锚定）。

## 5. 不做

- 错配/gap 比对、MAPQ 模型、paired 利用（2_insert_size 的 bbmap 是
  完整比对需求，属另一项，暂不做）。
- 种子多策略（首 k-mer 已充分：完美匹配下任意 k-mer 必命中）。
- SAM/BAM 排序、bam 输出（模板只要计数 + 派生覆盖度）。

## 6. 基准与优化评估（2026-08-11）

`benches/asm_map_benchmark.rs`（Lambda 2k reads 对 tadpole contig golden，
release）基线：

| 项 | 耗时 |
|---|---|
| build_index（67 contigs 参考） | ~1 ms |
| map 2k reads（无输出） | ~4 ms |

2026-08-11：basecov 已整体移出 map（见 §1 更新），不再有内存累计的
优化议题；派生管道的额外成本是重读一遍 mapped SAM（§2.5）。
