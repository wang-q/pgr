# sd / rept / align 代码与文档审核记录（2026-08-03）

对新增命令族（sd / rept / align）的代码与文档进行七轮审核。范围：
sd（8 命令 + libs/sd）、rept（6 命令 + libs/pl）、align（pgi/lastz +
libs/pgi、libs/lastz），约 8000 行代码逐文件审完；文档
docs/{sd,rept,align-pgi,align-lastz}.md 全部核对。最终 956 测试全绿、
clippy 干净。

## 修复的缺陷（9 处）

### 命令层

1. **sd/run.rs 越界 panic**：解析 `elem.bed` 时短行直接取 `f[4]`，越界
   panic。加 `f.len() < 8` 检查（与 cover.rs 一致）。
2. **repeat.rs 两处吞 IO 错误**：`map_while(Result::ok)` 在读取中途出错时
   静默截断（e_align 的 PSL 过滤、run_profex_per_chr 的 Profex 输出）。
   改为 `let line = line?;` 传播错误。
3. **trf.rs 特殊字符文件名**：`fa split name` 生成 `sanitize(name).fa`，
   trf 却用原始名拼 `${chr}.fa`——染色体名含 `/\():` 或双下划线时找不到
   文件。改用 `sanitize_filename(chr)`。
4. **e_align 参数校验缺失**：`--min-identity` 无范围校验（>1 全拒、<0
   全过）；kmer/smer/window/parallel 无正值校验。已加 `(0,1]` 与 `>0`
   校验。
5. **sd run/search/cross 的 `--min-identity`**：同样无范围校验，三处统一
   加 `(0,1]` 校验。

### libs 层

6. **lastz 静默失败**：`run_lastz` 对 lastz 失败只打日志、返回 Ok——所有
   job 失败时调用方拿到空结果无提示。改为统计失败数并 `bail`（实测损坏
   输入报 `lastz failed for 1 of 1 jobs`）。
7. **pgi build u32 溢出**：`pos_start: positions.len() as u32` 在 >42 亿
   k-mer 记录时静默截断索引。加 `payloads.len() <= u32::MAX` 防御检查。

### 文档一致性

8. **s_align soft-mask 说明缺失**：命令有 soft-mask 警告行为但帮助文本与
   用户文档未提及（e_align 提了）。两处补齐。
9. **docs/rept.md 缺 e-align 章节**：命令已实现但用户文档无命令文档
   （用户指出后补齐，属文档一致性修复）。

## 定位记录（未改，待决策）

1. **tube 工作流对"库 vs 基因组"结构性失效**：对照实验（酵母 + repbase，
   1051 万种子一致）greedy 出 2220 个 PSL 块、tube 只有 4 个。根因：
   tube 的 merge 只在相邻对角桶对（宽 64 bp）间独立进行、每桶对单独累计
   覆盖（MIN_COV=85）；库比对种子稀疏，跨桶对的链被切断、cov 不足被丢。
   FastGA 的 tube 面向高密度全基因组自比对。e-align 默认 greedy 正常；
   修复需改 tube 跨桶 merge，属 pgi 算法改动，暂不处理。
2. **60,423 → 75,413 数据勘误**：e_align 对 MG1655+tncentral 历史记录
   60,423 bp，git archive 独立编译 c17a3d0 验证当前代码稳定输出
   75,413（1.63%）——差异来自 tncentral 库 16:30 更新（6073 → 6093 条，
   嵌入 header 拆分）及编译产物时序，非代码 bug。笔记
   [[repeat-masking.md]] §2.3.5 已加勘误。

## 鲁棒性验证（无 panic）

* 全 N 序列：e-kmer/s-kmer exit 1（合理），align pgi / e-align 正常退出；
* 截断 .pgi 索引 → 友好报错 "truncated index records"；
* `-p 0`（rayon 容错）、`--step 0`（fa window 自带校验）实测不 panic；
* 索引兼容校验、除零保护、负链坐标转换、95% 去重、重叠分母防护均确认。

## 文档审核结论

* docs 引用的 14 个 sd/rept/align 命令全部真实存在；
* 参数默认值与 CLI 一致（sd、align-pgi、align-lastz、rept 5 命令）；
* 示例语法与 CLI 一致；`--self`/`--syn` 引用正确；
* 补齐项：rept.md 的 e-align 章节、s-align soft-mask 说明、SD 搜索前勿用
  自比对提醒、TnCentral 示例路径修正。

## 低风险记录项（未改）

* `ctx.rs` 的 `tempdir.path().to_str().unwrap()`（非 UTF-8 路径罕见）；
* `decompose.rs` 负链投影 `gend - end` 依赖 header 与序列长度一致
  （cluster 内部保证）；
* cluster/cover 的 u32→i32 坐标转换（仅 >2.1 Gb 染色体才溢出）；
* rept 命令族不预检输入文件存在性（依赖子进程报错，可接受）。

## 回归保护

新增 5 个集成测试：trf 特殊字符名、e_align 非法 identity、s-align 端到端、
sd run 端到端、pgi 溢出检查路径。全套 956 通过。

## 复核轮（第八轮，2026-08-03）

对上述 9 处修复逐条对照当前代码复核，全部在位；另发现并修复 3 处新缺陷、
2 处文档错误。最终 958 测试全绿（956 + 新增 2），fmt 与
`cargo clippy -- -D warnings` 干净。

### 新修复（3 处）

1. **cluster.rs 重叠 union 漏连嵌套区间**：按 start 排序后仅用 `windows(2)`
   相邻对做 union-find。当 A 包含 B、C 且 B/C 不相交时（A 与 B、C 均重叠），
   C 只与 B 相邻比较，漏掉与 A 的直接重叠，被拆成独立簇。改为扫描时跟踪
   已见区间的最大右端：任何与先前区间重叠的当前区间必然与最大右端区间重叠。
   新增回归测试 `nested_overlapping_intervals_form_one_cluster`。
2. **align lastz 省略 query 未真正启用 self 模式**：`self_mode`（省略 query
   或 `--self`）与传给 `RunLastzOptions` 的 `is_self` 分离——省略 query 时
   `is_self` 仍为 false，目录输入会跑 n×n 全交叉、单文件不传 `--self`
   （产生序列 vs 自身的平凡比对），与文档"省略 query 即 --self"相悖。改为
   传 `self_mode`；`--self` 语义（丢弃平凡自比对）与 `align pgi` self 模式
   一致。新增集成测试 `command_align_lastz_omitted_query_is_self`；既有
   `test_align_lastz_single_input_self` 断言的是旧行为（"s {" 平凡 hit），
   改为断言 LAV 头部记录 `--self`。
3. **`--min-identity` 校验范围与文案不符**：sd run/search/cross 与 rept
   e_align 四处写 `(0.0..=1.0)`（放行 0.0）但报错文案称 `(0, 1]`。统一收窄
   为 `> 0.0 && <= 1.0`。

### 文档修正（2 处）

1. align-lastz.md：`-s, --preset` 的 `-s` 短选项实际不存在（clap 未注册），
   改为 `--preset`。
2. align.md：lastz 示例 `-o out.psl` 错误（lastz 输出为目录下的 LAV），改为
   `-o lastz_out` 并注明输出为 LAV、需 `pgr lav to-psl` 转换。

### 新增定位记录（未改）

1. **decompose.rs 静默丢弃无法解析的 FASTA 头**：非 cluster 格式的普通 FASTA
   输入逐条静默跳过（无警告、输出空），与 `parse_or_warn` 风格不一致；建议
   后续加警告。
2. **e-align PSL 过滤静默跳过畸形行**：`let Ok(psl) = ... else continue`
   无日志；sd 引擎用 `parse_or_warn`。建议统一。
3. **run_lastz self 模式仍构建 n×n job 列表**：内存 O(n²)，运行时才按
   basename 跳过；大目录 self 对齐可提前过滤。
4. **pgi u32 溢出防护无测试**：需 >42 亿 k-mer 记录，不现实；保持代码审查
   覆盖。
5. **e-align 的 identity 定义**：用 `Psl::ident()`（gap-compressed，不含
   insert 碱基），与 sd 文档化的 `(matches+repeats)/block_len` 不同，rept.md
   未说明；建议补注。
6. **wave.rs 的 panic**（`unreachable!` / `panic!`）均为 Myers 算法不变量，
   有穷举与随机测试兜底；未发现用户输入可触发的路径。

### 复核验证结果

* 原笔记 9 处修复全部对照确认在位（run.rs 越界检查、repeat.rs 两处吞 IO
  错误、trf 特殊字符文件名、e_align 与 sd 的参数校验、lastz 失败统计 bail、
  pgi u32 溢出防护、s_align soft-mask 文档、rept.md e-align 章节）。
* `cargo fmt --check` 干净；`cargo clippy -- -D warnings` 干净；
  `--all-targets` 有 3 个既有 `src/libs/plot/dot.rs` 测试告警（plot 模块，
  不在本次范围，未改）。
* 全量 958 测试通过（lib 463 + 各集成测试；新增 cluster 单测与 align-lastz
  集成测试各 1）。
