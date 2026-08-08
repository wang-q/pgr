# pgr fas 命令族代码审核记录（2026-08-05）

对 `pgr fas` 命令族（全部 20 个子命令）以及相关库文件（`libs/fmt/fas`、
`libs/alignment`、`libs/fas_multiz`、`libs/fas_xlsx`、`libs/nt`、`libs/io`、
`libs/ds/{intspan,crossover}`）和全部测试/文档进行审核。以下仅保留有借鉴意义
的结论；逐轮验证过程已精简。

审核范围：
- **信息**：`check` / `cover` / `link` / `name` / `stat`
- **子集**：`filter` / `slice` / `subset`
- **转换**：`concat` / `consensus` / `join` / `multiz` / `refine` / `replace`
- **文件**：`create` / `separate` / `split`
- **变异**：`to-vcf` / `to-xlsx` / `variation`

审核重点：数据安全（`-o` 不得覆盖输入，含 `.loc` 侧车索引与分块/逐物种输出）、
Zero Panic（畸形输入不 panic）、坐标/长度边界处理、算法正确性、文档一致性。

## 排除的疑点（安全不变量，经核验无需修复）

- 全部 `unwrap()`/`unreachable!` 均为 clap `required` 参数或 `value_parser` 约束
  枚举，运行期不可达。
- `-o` 覆盖保护覆盖全部单文件输出命令（`ensure_outfile_distinct`）；`separate`/
  `split` 输出为目录，采用逐输出路径 `same_path` 反向检查。
- **POA 回溯 `curr_j - 1` 下溢**（`poa/align.rs` E 态）：进入 E 态时恒有
  `curr_j >= 1`——所有模式下 `e[*][0]` 均为 `neg_inf`，任何转移进 E 态都要求
  `e[*][0]` 等于真实分值，恒不可能。经 3000 余随机图 × 3 模式 fuzz 无 panic。
  判定为不可达，不加守卫（避免推测性防御代码）。
- **`banded_align_refs` traceback `i-1`/`j-1` 下溢安全**：`node` 恒为被选为最优的
  有限前驱状态 ⇒ 满足其自身 `i>0`/`j>0` 守卫；row-0 仅插串链且 trace 恒为 I 标志，
  回退到 `(0,0)` 时循环条件已为假而退出。带约束单调保证 `cell(i,0)` 在带内。
- `best_crossover` 的 `debug_assert_eq!`：四个切片均以 `map_a.len()` 构建，长度恒
  相等，release 下为 no-op。
- `align_seqs_quick` 等长校验：quick 模式前提是已对齐输入（axt/maf 转换场景），
  对歧异长度输入报错而非 panic，属正当行为。
- `Range::from_str` 大坐标溢出：生产路径 `decode`→...→`parse_i32` 对溢出数字串
  返回 `None`，不 panic。
- `write_vcf_block` 的 `pos_idx`：`checked_sub(1)` 双重校验，无 panic。
- `fas_xlsx` 绘制逻辑：颜色索引 `fold % color_loop` 恒落在已注册区间（`--colors`
  限定 `[1,15]`），格式键必然存在。
- 共享 IO 辅助：`read_runlist`（拒绝空 run）、`read_names`（跳过空行与 `#`）、
  `same_path`（canonicalize + dev/inode）、`PgrWriter`（Drop 不 panic）。

## 记录项（未改，低风险 / 待决策）

- `separate`/`split` 的物种名/染色体名经 `sanitize_filename` 清洗后，两个不同名称
  可能碰撞到同一输出文件名（如 `a/b` 与 `a_b`），静默合并到同一文件。与文件名清洗
  方案固有行为一致。
- `run_pipeline` 串行分支对 `proc_block` 错误直接传播（`?`），并行分支收集后统一
  报错；对畸形 block 的 `next_fas_block` 错误则跳过并告警。行为一致。
- `to-xlsx` 非 outgroup 下 `sub.freq = min(freq, N-freq) <= 0.5` 恒成立，
  `--min-freq>0.5` 会过滤掉全部 SNP 输出空表。属参数语义边界。
- `merge.rs` 非参考物种回退仅在"完全未映射"时生效，对"映射到 gap"不回退到另一块；
  参考物种刻意只取 block A（保持 ungapped 不变量），非参考物种同样"第一块优先"，
  属设计取舍。
- 理论溢出（不可达）：`banded_align` 的 i32 乘法需块物种数 >8400；`poa` 的 u32
  权重需 >4e8 次累加；`topological_sort` 重叠非互斥 clique 仅手工构造畸形图可触发。
- `cover --name` 在"输入无块"（报错）与"有块但无该物种"（输出空 `{}`）行为不一致，
  均无数据损坏。
- `link --best` 函数 doc 措辞与实现略有出入，但 CLI 帮助文本已准确描述，行为与文档
  一致。

## 已知限制（有意保留）

- **`multiz` 同一输入文件在某窗口内含多个非连续 block 时，`merge_window` 只取第一
  个重叠 block，其余静默丢弃**。正确修复需按参考坐标拼接非连续 block，而现有 DP/
  `merge_conflicting_refs` 是为重叠 locus 的 re-align 设计，引入拼接逻辑风险高，
  超出当前合并架构。仅在 `--radius` 大到使两个非连续 block 的 `±radius` 扩展并入
  同一窗口时触发；默认 `--radius 30` 通常不会。
- **`multiz` `merge_window` 非 DP 回退在参考序列逐字符不等时丢弃整窗口**：DP 合并
  失败回退简单拼接时要求所有 block 参考序列 `entry_seq_equal`，否则对整窗口返回
  `None` 静默丢弃。参考不一致时本就无法单义拼接，属设计取舍，与"非连续 block 数据
  丢失"同源。
- **`refine --outgroup` 的 trim 删除"内群全 gap、外群有碱基"的列后，外群 range 坐标
  未随之收缩**，输出的 range 与序列长度不一致。与 kent 原版一致，修复需贯通 trim
  签名。
- **`merge_conflicting_refs` 合并后参考序列与其声明 range 的 ungapped 长度可能不一
  致**：仅在"同一窗口两个输入 block 参考序列真实不一致"（矛盾输入）时触发，属共享
  坐标系合并架构下的设计取舍。

## 修复的缺陷（共 35 处，根因模式）

### 数据安全（`-o` / 输出路径同输入保护）

- **流式命令允许 `-o` 覆盖输入**：先打开输出再读取输入，`-o` 指向输入会在读取前
  截断。修复：`stat`/`variation`/`filter`/`join`/`name`/`consensus`/`refine`/
  `multiz`/`replace`/`slice`/`subset`/`concat`/`to-vcf`/`to-xlsx`/`link`/`cover`
  统一加 `ensure_outfile_distinct`；涉及辅助输入的命令同时保护 `--required`/
  `--runlist`/`--replace-tsv`/`--sizes`。
- **`create` 的 `-o` 可覆盖基因组 `.loc` 侧车索引**：先开（截断）writer 后
  `open_indexed` 读索引；`-o` 名为 `{genome}.loc` 会先截断，随后 mtime 判定"新鲜"
  不重建，`.loc` 被永久损坏。修复：`{genome}.loc` 一并加入保护列表。`check` 同理。
- **`separate`/`split` 输出文件可能覆盖输入**（输出为目录、文件名动态生成）。修复：
  打开每个输出前用 `same_path` 与所有输入比对，命中即 `bail!`。

### Zero-Panic / 越界

- **不等长 block 越界 panic**：`get_subs`/`get_indels`/`alignment_stat`/
  `align_seqs_quick` 假设列数一致，不等长时复用 `seqs[0].len()` 索引越界。修复：
  函数开头校验所有序列等长。
- **外群序列短于内群越界 panic**：`polarize_subs`/`polarize_indels` 用 sub 坐标直接
  索引外群序列。修复：加 `ensure!(og_idx < og.len())` / `ensure!(end <= og.len())`。
- **`slice_block` 对全 gap 第二条序列 panic**（子切片落在 indel 岛内时 `ss` 为空，
  `ss.min()`/`ss.max()` 对空 `IntSpan` panic）。修复：`ss_ints.is_empty()` 时跳过。
- **`trim_head_tail` 对全 gap 比对 panic**（`--chop` 时 `remove(cur_len - 1)` 在
  `cur_len == 0` 下溢）。修复：头部移除量 `min(..., len)`，尾部序列为空时 `break`。
- **外部对齐器返回越界序列 id 时 panic**（`refine --engine` 非 builtin 解析外部 MSA
  输出，`>N` 的 `N` 直接索引 `out_seqs`）。修复：写入前校验 `idx >= out_seqs.len()`
  返回友好错误。
- **`multiz` 窗口边界 off-by-one**：`ref_overlaps_window` 原用 `start < window.end &&
  end > window.start`，落在窗口边界上的单碱基 block 被判不重叠。修复：改为含端点的
  `start <= window.end && end >= window.start`。
- **`IntSpan` runlist 解析空 token 静默注入坐标 0**（前导/连续逗号退化为 `(0,0)`）。
  修复：遇到逗号开头的空 run 返回 "empty run" 错误；尾随逗号仍作无害换行终止符。
- **`align_to_chr` 对空序列 intspan panic**（对空 `IntSpan` 调 `min()`/`max()`）。
  `fas` 现有调用方均不可达该路径，但 `align_to_chr` 是公开库函数、与对称的
  `chr_to_align`（返回错误）不一致。修复：入口对空 intspan 返回友好错误。
- **`multiz` 对倒置参考区间（`start > end`）u64 下溢 panic**（`width = e - s`）。
  修复：`derive_windows_from_blocks` 跳过 `start > end` 的倒置参考 entry。

### 算术溢出

- **`to-xlsx` 序列数 >32 时颜色索引溢出**：`u32::from_str_radix` 解析与内群序列等长
  的二进制串，>32 时超出 u32 位宽返回溢出错误。修复：改为按字节折叠模 `color_loop`
  计算索引（对短串结果一致、对长串无溢出）。
- **`to-xlsx` `--wrap` 大值导致 u16 溢出**（`wrap + 3` 回绕）。修复：
  `saturating_add(3)`；`paint_indel`/`paint_sub` 的游标推进同理防溢出。

### 算法正确性

- **`banded_align` I 态 gap 延伸成本用错计数**：`k` 被外层循环单元格索引遮蔽，误用
  为列内物种数。修复：变量改名 `k_a` 并用于全部 gap 延伸计算。
- **POA 共识忽略节点权重，多数碱基不敌首序列骨架**：`generate_consensus` 只用边
  权重忽略 `NodeData.weight`。修复：按 SPOA heaviest-bundle 语义把节点权重计入路径
  得分。
- **`to-xlsx` 列游标推进错误导致单元格重叠/数据损坏**（固定 `col_cursor += 1`，但
  多碱基 indel 实际占用 `length.min(3)` 列）。修复：`paint_sub`/`paint_indel` 返回
  实际占用列数，游标按其推进。

### 对齐语义

- **`refine_block` 误加等长校验，破坏 `refine` 的不等长重比对语义**：`refine` 用途
  恰是重比对不等长序列，等长校验使合法输入报错。修复：**移除**该等长校验（`refine`
  走 `align_seqs` POA/外部队列；`--chop` 由 `trim_head_tail` 空序列保护兜底）。
  `align_seqs_quick` 的等长校验保留（quick 模式前提是已对齐输入）。

### 坐标 / 共识

- **`slice` 对负链参考/物种产生空输出或反向范围**（`chr_to_align`/`align_to_chr` 使
  `ss_start > ss_end`）。修复：子切片与各物种输出时把 `start/end` 交换归一为
  `start <= end`（保留原链向）。
- **POA 共识对全 gap 首序列输出空**（首序列全 gap 时 gap 节点权重最小胜出，共识退化
  为全 gap）。修复：`consensus_block` 生成共识前剔除全 gap 序列，使只有真实碱基参与
  投票。

### 集合与去重

- **`concat`/`subset` 对 `--required` 中重复物种名输出重复条目**。修复：读取
  `needed` 后按首次出现顺序去重（`HashSet` 判重）。
- **`subset` 块内重复物种名静默丢序列**（`HashMap` collect 保留最后一次）。修复：
  `entry.entry(name).or_insert(e)` **首现优先**，与 `concat` 首匹配语义对齐。

### 多序列合并（multiz）

- **`multiz` 渐进合并失败时静默丢弃物种**。修复：`match` 捕获 `None` 分支，
  `log::warn!` 提示被丢弃的 block 及其物种。
- **`multiz` 覆盖率过滤误删单碱基参考块**（`DupeTree::add` 忽略零宽单碱基区间）。
  修复：移除该覆盖率过滤（每个窗口本就从至少一个输入的参考区间派生，覆盖恒成立）。
- **`merge_conflicting_refs` 对仅存在于单块的物种静默截断**（在 crossover 处按
  `pos < cut` 用 map_a、`pos >= cut` 用 map_b，单块物种半段全塌成 `-`）。修复：按
  物种**在哪些块中存在**选择映射——共享物种保持 `cut` 拼接，仅存在于 A 的用 `map_a`
  全程、仅存在于 B 的用 `map_b` 全程，不再截断。

### 边界 / 非 IUPAC

- **`create` 对超出参考基因组范围的坐标 abort**（与"chr 不存在跳过并告警"不一致）。
  修复：`create_from_links` 捕获 `"slice error"` 开头的错误，记录警告并 `continue`。
- **`separate --rc` 对非 IUPAC 字符反向互补后乱码**（`NT_COMP` 未知字符 → 255 哨兵
  渲染为 `ÿ`）。修复：`separate` 内联反向互补，`NT_COMP` 返回 255 时保留原始字节。
- **`slice` 对带 gap 的参考物种中止整个命令**（`chr_to_align` 边界小于 runlist 上界
  触发 `?` 传播）。修复：捕获 `chr_to_align` 两端点任一失败，`log::warn!` 并
  `continue` 跳过该子区间。

### to-xlsx 显示 / 样式

- **`to-xlsx` 外群含 IUPAC 歧义碱基时命令失败**（样式名 `sub_{obase}_unknown` 只对
  标准碱基注册）。修复：该分支加 `.or_else(|| format_of.get("sub_N_unknown"))` 兜底。
- **`paint_indel` 换行 off-by-one 与小节双重计数**（内部换行后调用方又后置换行，
  `--wrap 1/2/3` 时 `sec_cursor` 连续递增两次留空小节）。修复：换行预检上移到调用方
  统一处理（`col_cursor + width > wrap + 1` 才换行），移除内部换行。

### 文档一致性（一次性小修，已精简）

`consensus` 输出格式描述、`calling 变异` 中英混杂、`collect_subs`/`collect_indels`
外群语义 doc、`trim_head_tail` doc 与实现不符等均已修正。

## 结论

`fas` 命令族审核完成（累计修复 35 处缺陷并补回归测试与文档澄清），经多轮纵深复核
收敛，未再发现需要修复的新缺陷。
