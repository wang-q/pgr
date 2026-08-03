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
