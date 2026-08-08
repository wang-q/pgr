# `pgr dist` 命令族代码与文档审核记录（2026-08-08）

对 `pgr dist` 命令族（`mini` / `mash` / `frac` / `hv` / `pgi`）及相关库
（`libs/hash`、`libs/hv`、`libs/syncmer`、`libs/pgi/dist`、`libs/pgi/to_hv`、
`libs/par`、`cmd_pgr/args` 的 dist 共享参数）及文档进行审核，修复缺陷并补
回归测试。以下仅保留有借鉴意义的结论；逐轮验证过程已精简。

## 修复的缺陷（根因模式）

### 数值语义（含方向/空输入/无关系列/边界）

- **`dist pgi` containment 用 `min(total1,total2)` 作分母（对称），与命令族约定不一致**：
  其余子命令（mini/mash/frac/hv）与 `docs/dist.md` 均定义为
  "Containment = 交集 ÷ 第一个集合大小"（有方向）。修复：分母改回第一个索引
  `total1`，与命令族对齐。回归 `command_pgi_dist_containment_directional`。
- **空输入 / 空==空 数值语义不一致**：`dist hv` 对空输入 0/0 输出 NaN；`dist pgi`
  空==空给出 jaccard=0/距离 1（应视为完全相同距离 0）。修复：对 `union==0`
  （空==空）统一给 jaccard=1、mash=0；空首集合 containment=0；`.hv` 路径双空
  视为 sim=1。与 sketch 家族语义对齐。回归 `command_dist_hv_empty_inputs_no_nan`、
  `command_pgi_dist_empty_indexes`。
- **`dist hv` FASTA 路径对无关系列输出负 Jaccard/Containment**：`inter` 是零均值
  点积，对无关系列可随机为负，`union` 随之偏大、jaccard/containment 为负（数学
  上无效）。`.hv` 文件路径因 `inter` 经 `usize` 转（饱和为 0）已天然 clamp，仅
  FASTA 采样路径受影响。修复：`inter = hv_dot(s1,s2).max(0.0).min(card1).min(card2)`。
  回归 `test_calc_distances_disjoint_no_negative`。
- **`mash_distance` 注释声称 "bounded to [0,1]" 但小 k 时可超 1**：实现仅对
  `jaccard==0` 返回 1，未 clamp。修复：clamp 到 `[0,1]`（与 Mash `min(1,dist)`
  一致）。回归 `test_mash_distance_bounded_to_one`。

### Zero-Panic / 参数校验

- **`--hasher mod` 配偶数窗口触发 minimizer_iter 断言 panic**：`seq_mins` 的
  `mod` 分支用 `MinimizerBuilder::new_mod().width(opt_window as u16)`，而该
  builder 在 `iter()` 时**断言窗口必须为奇数**。显式偶数 `-w`（如 `-w 4`）直接
  panic，违反 Zero-Panic。默认窗口均为奇数故默认不触发。修复：`mod` 分支前校验
  `opt_window % 2 == 1`，偶数返回友好错误。回归
  `command_dist_mini_mod_rejects_even_window`。

### 文档与实现一致性

- **`--hasher` 对 `--sampler syncmer` 被静默忽略**：syncmer DNA 恒用 2-bit 滚动
  哈希、蛋白恒用 RapidHash，完全不使用 `--hasher`；但帮助/注释暗示其生效。修复：
  明确 `--hasher` 仅对 minimizer/FracMinHash 生效，syncmer 忽略之。
- **`--dim` 帮助称"需为 32 的倍数"但实现不校验、任意正数可用**：`hash_hv_bit`/
  `hash_hv_sparse` 支持任意维。修复：改为"建议为 32 的倍数（对齐/性能），实现不
  强制"。

## 记录项（未改，低风险 / 待决策）

- `seq_mins` 的 `mod` hasher 路径用 `opt_window as u16` 传给
  `MinimizerBuilder::width`。若 `--window > 65535`（仅显式超大 `-w` 可达）会静默
  截断为 `w & 0xffff` 而不报错。属理论性边界，未改。

## 结论

`dist` 命令族审核收敛。第 4 轮复核（全部子命令、共享参数 `resolve_kmer_window`、
`--parallel` 1..=1024 校验、输出列格式、文档）确认无新缺陷；`cargo build`/
`clippy`/`fmt` 干净，`cli_dist`/`cli_dist_mini`/`cli_dist_mash`/`cli_dist_frac`/
`cli_pgi` 及 `--lib hash::`/`hv::`/`pgi::dist` 全部通过。
