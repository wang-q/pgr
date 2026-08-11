# paf 命令族代码审核记录（2026-08-05）

对 `pgr paf` 命令族（`index` / `query` / `to-bed` / `to-fas` / `to-maf` /
`to-vcf` / `to-gfa` / `graph` / `stat`，共 9 个子命令）及其库文件
（`libs/paf`：`parser`、`record`、`cigar`、`index`、`fasta`、`query`、
`msa_build`、`to_maf`、`to_fas`、`vcf`、`poa_compact`、`graph`、
`maf_import`）和全部测试/文档进行审核。以下仅保留有借鉴意义的结论；
验证过程已精简。

## 排除的疑点（安全不变量，经核验无需修复）

- `query` 的 subset / syntenic 过滤中 `idx.id_to_name(*qid).unwrap_or("")` 带
  默认值、`output_paf` 用 `unwrap_or("?")`，均非 panic。
- `index/query.rs::project` 空 CIGAR 分支：`len <= 0` 返回 `None`；带 CIGAR 分支仅
  在 `q_min < q_max && t_min < t_max` 才 `Some`，仅覆盖 I/D 的区间正确返回 `None`。
- `output_paf` 的 `matches`/`block_len`/`gi`/`bi`/`cg` 取完整源 CIGAR 而非投影子
  区间，属文档已声明的已知限制（`libs/paf/query.rs` 头注释），非 bug。
- `--merge-distance > 0` 时强制要求 `--fasta-tsv`（`anyhow::bail!`），与
  `docs/paf.md` 一致。

## 记录项（已改 / 未改）

- **`output_paf` 的 `query_length`/`target_length`（PAF 第 2/7 列）已修复**
  （2026-08-11）：`PafIndex` 新增 `seq_lens`（建索引时从 PAF 列 2/7 收集，
  first-wins），持久化格式 v4 → v5（旧 `.paf.idx` 报"请重建"）；输出列
  填充真实长度。测试：libs roundtrip + output_paf 长度断言；cli_paf 断言
  从 `A\t0\t` 更新为 `A\t100\t`。
- `to-vcf` 的部分删除样本回退为 best-effort 等位基因且不做完全左对齐
  （`docs/paf.md` 与 `after_help` 已声明）。
- `query --transitive --merge-distance` 依赖 `--fasta-tsv` 重算 CIGAR；若用户未
  提供 `-f` 会硬报错，属有意约束。

## 修复的缺陷（根因模式）

### 坐标投影正确性（Zero-Panic / 精度）

- **空 CIGAR 投影左越界时 over-count 长度**（`index/query.rs::project`）：查询起点
  落在记录 target 区间左侧时，原按查询区间原始宽度投影，投影到 query 上过长。
  修复：将查询跨度 clamp 到实际重叠 `[ts,te)` 与记录 target extent 的交集。
- **空 CIGAR 且 `-` 链的子区间投影坐标错误**：`-` 链 target 偏移应映射到 query
  的*末端*侧 `[query_end - off - len, query_end - off)`，原实现未镜像（与带 CIGAR
  的 `rc_to_forward` 不一致）。修复：空 CIGAR 分支对 `-` 链做末端镜像。
- **`maf to-paf` 反向链区间越界可能 usize 下溢**（`maf_import.rs`）：`-` 链
  `start + size > src_size` 时 `src_size - end` usize 下溢。修复：调用前
  `checked_add` + 越界 `bail!`。

### 数据安全（`-o` 覆盖保护）

- **`paf` 家族各输出命令此前未调用 `ensure_outfile_distinct`**，`-o` 指向输入会
  静默覆盖。修复：`index`（保护全部 infiles）、`query`/`to-bed`（保护 infile）、
  `to-fas`/`to-maf`/`to-vcf`/`to-gfa`（保护 infile + fasta_tsv）、`graph`/`stat`
  （保护 infile + fasta_tsv）就地补齐。
- **查询类子命令未保护 `--subset-sequence-list` / `--syntenic-filter` 辅助输入**
  （`query`/`to_bed`/`to_fas`/`to_maf`/`to_vcf`/`to_gfa`）：`-o` 指向过滤文件会把
  辅助输入覆盖成输出。修复：6 个命令输入列表补上 `subset_list`/`syntenic_filter`。

### 参数校验（CLI）

- **POA 评分参数接受非法值**（`args.rs`）：`--match` 可负数/0、罚分可正数，POA
  引擎静默产生逆向/无效打分。修复：`--match` 必须 > 0、`--mismatch`/`--gap-open`/
  `--gap-extend` 必须 ≤ 0，clap 解析期拒绝。
- **查询过滤参数接受负值被静默当作"关闭"**（`--min-dist`/`--min-output-len`/
  `--merge-distance`/`--min-chain-length`）。修复：`parse_non_negative_i32` 解析期
  拒绝。
- **`--min-identity` 越界未被校验**（此前 `value_parser!(f64)` 接受任意浮点）。
  修复：`parse_min_identity` 拒绝 `[0.0,1.0]` 之外的值。

### 变体正确性 / 对齐

- **VCF 中 DEL 在 INS 之后被静默丢弃**（`vcf.rs`）：连续 target gap（INS）后紧跟
  target 非 gap 但有 query gap（DEL）时，旧代码取 `msa[0][col_start-1]`（恰为 INS
  的 target gap）作锚点并丢弃该 DEL。修复：DEL 锚点取 INS 之前最后的 target 非
  gap 碱基（经 `t_aln_pos`）。
- **MSA 构建对相接片段不合并，导致 per-name 去重丢碱基**（`msa_build.rs`）：
  合并判据用 `qs < last.2`，target 删除把连续 query 切成恰好相接的两片
  （`qs == last.2`）时不会合并，去重静默丢弃第二片碱基。修复：判据改 `qs <=
  last.2`。

### 输出坐标 / 图节点序列 / 命令可用性

- **`to-fas` pairwise `-` 链头坐标用正向区间而非反向坐标**：头显示的应是反向
  坐标（`src_size - qe`），与 `to-maf`/MSA 输出一致。修复：改用 `q_start_maf`
  推导显示坐标。
- **`graph` 节点序列未按代表段正向取向存储**：DSU 代表段缺失、由 `-` 取向成员段
  填充时直接存正向区域，`-` 步长与存储序列不一致。修复：由 `-` 段填充时反向互补。
- **`to-bed` 未接受 `--fasta-tsv`，导致 `--merge-distance` 无法使用**：`to-bed`
  的查询参数来自共享 `add_query_args`（含 `--merge-distance`），但命令未注册
  `-f`，而 `run_query` 在 `--merge-distance > 0` 时强制要求 `-f`。修复：为 `to-bed`
  补 `add_optional_fasta_tsv_arg` 并纳入 `ensure_outfile_distinct`。

## 文档修复

`docs/paf.md` 引言"查询类子命令"列表补上遗漏的 `to-fas`。

## 结论

`paf` 命令族审核完成（累计修复 15 处缺陷、补回归测试与文档澄清），经多轮纵深复核
收敛，未再发现新问题。
