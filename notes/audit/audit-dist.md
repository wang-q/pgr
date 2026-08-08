# `pgr dist` 命令族代码与文档审核记录

对 `pgr dist` 命令族（`mini` / `mash` / `frac` / `hv` / `pgi`）及相关库
（`libs/hash`、`libs/hv`、`libs/syncmer`、`libs/pgi/dist`、`libs/pgi/to_hv`、
`libs/par`、`cmd_pgr/args` 的 dist 共享参数）及文档（`docs/dist.md`、
`notes/design/hv.md`、`notes/benchmarks/*`）进行审核。缺陷按轮次记录；
关键结论均经实测验证（详见"验证"一节）。

## 第 1 轮（2026-08-08）

### 代码 / 行为缺陷

#### D1. `dist pgi` 的 containment 用 `min(total1, total2)` 作分母，与文档及命令族不一致
- 位置：[dist.rs](file:///home/wangq/Scripts/pgr/src/libs/pgi/dist.rs#L86-L90)
  `containment = inter / total1.min(total2)`。
- 其余子命令（mini/mash/frac/hv 的 FASTA 路径与 `.hv` 路径）及
  [docs/dist.md](file:///home/wangq/Scripts/pgr/docs/dist.md) 均定义为
  "Containment = 交集 ÷ 第一个集合大小"（有方向，以第一个输入为分母）。
- 实测（`s1.pgi` 小、`s2.pgi` 大，共享 5752）：
  - `dist pgi s1 s2` → containment 0.5024
  - `dist pgi s2 s1`（交换参数）→ containment 仍 0.5024
  - 即 pgi 的 containment 是**对称**的（min 分母），而命令族约定是**有方向**的。
- 影响：`dist pgi` 的 containment 语义偏离文档与命令族；文档未说明该偏差。
- 建议：与命令族对齐改为 `inter / total1`（第一参数为分母）；或在文档中
  显式说明 pgi containment 用较小集合作分母、为对称口径。

#### D2. `dist hv` 空输入输出 NaN（FASTA 路径与 `.hv` 路径都受影响）
- 位置：[hv.rs](file:///home/wangq/Scripts/pgr/src/libs/hv.rs#L372-L377)
  `calc_distances` 中 `union = card1+card2-inter` 可为 0，`jaccard =
  inter/union` 与 `containment = inter/card1` 出现 0/0 → NaN；
  [hv.rs](file:///home/wangq/Scripts/pgr/src/cmd_pgr/dist/hv.rs#L235-L239)
  `run_hv_files` 中 `union`/`jaccard` 同样 0/0 → NaN。
- 实测：`dist hv empty.fa empty.fa` → `0 0 0 0 NaN NaN NaN`；
  `dist hv empty.fa t.fa` → 同样 NaN。
- 对比：`dist mini/mash/frac` 对空==空正确处理为 jaccard=1、距离 0；
  违反"Zero-Panic / 不输出 NaN"的项目原则。
- 建议：`calc_distances` 与 `run_hv_files` 对 `union==0`（空==空）给
  jaccard=1、mash=0；对空首集合的 containment 给 0（与 sketch 家族一致）。

#### D3. `dist pgi` 空==空给出 jaccard=0 / 距离 1，与命令族不一致
- 位置：[dist.rs](file:///home/wangq/Scripts/pgr/src/libs/pgi/dist.rs#L81-L85)
  `union==0 → jaccard=0`。
- 实测：两个 0 unique k-mer 的 `.pgi`（如 4 nt 序列）`dist pgi a.pgi a.pgi`
  → `0 0 0 0 1.0000 0.0000 0.0000`（完全相同但距离 1）。
- 对比：sketch 家族与 `dist hv` 对空==空视为完全相同（距离 0）；
  pgi 自身的 `identical_indexes_zero_distance` 测试也断言相同索引距离 0
  （但该测试用非空索引，未覆盖空索引）。
- 建议：`union==0` 时应像 `set_distances`/`mash_sketch_distances` 一样
  返回 jaccard=1、mash=0。

#### D4. `mash_distance` 注释声称"bounded to [0,1]"，但小 k 时可超过 1
- 位置：[hash.rs](file:///home/wangq/Scripts/pgr/src/libs/hash.rs#L378)
  `SetDistances::mash` 注释写 "Mash distance (bounded to [0, 1])"；
  [hash.rs](file:///home/wangq/Scripts/pgr/src/libs/hash.rs#L129-L135)
  的实现仅对 `jaccard==0` 返回 1，未对 `>1` clamp。
- 实测：k=1、j=0.05 → 2.35；k=2、j=0.01 → 1.96。小 k 时 `dist mini -k 1`
  等会打印 `>1` 的距离，与"bounded [0,1]"注释不符。
- 建议：要么在 `mash_distance` 内 clamp 到 1.0（与 Mash 文档一致），
  要么修正注释。属次要问题。

### CLI / 帮助 / 文档缺陷

#### D5. `dist hv --sampler syncmer` 时 `--hasher` 被静默忽略，但帮助与注释声称其生效
- 位置：syncmer DNA 恒用 2-bit 滚动哈希、蛋白恒用 `RapidHash`（
  [syncmer.rs](file:///home/wangq/Scripts/pgr/src/libs/syncmer.rs#L248-L254)）；
  `dist hv` 的 syncmer 路径经
  [hv.rs](file:///home/wangq/Scripts/pgr/src/cmd_pgr/dist/hv.rs#L131-L148)
  走 `load_hv_from_fasta_syncmer`，完全不使用 `--hasher`。
- 但 [hv 帮助](file:///home/wangq/Scripts/pgr/src/cmd_pgr/dist/hv.rs#L15-L18)
  写 "Samplers, hash algorithms, ... are the same as the sketch-distance
  family"（暗示 hasher 对 syncmer 生效）；
  [args.rs](file:///home/wangq/Scripts/pgr/src/cmd_pgr/args.rs#L804-L807)
  `sampler_arg` 注释写 "protein syncmers hash s-mer bytes with `--hasher`"
  —— 均与实现不符（蛋白 syncmer 也硬编码 RapidHash）。
- 实测：`dist hv ... --sampler syncmer --hasher rapid/fx/mod` 输出完全相同。
- 建议：文档明确 `--hasher` 仅对 minimizer 采样生效、syncmer 忽略之；
  修正 `sampler_arg` 注释。

#### D6. `dist hv --dim` 帮助/文档称"需为 32 的倍数"，但实现不校验、任意正数可用
- 位置：[hv.rs](file:///home/wangq/Scripts/pgr/src/cmd_pgr/dist/hv.rs#L84)
  help "The dimension size should be a multiple of 32"；
  [docs/dist.md](file:///home/wangq/Scripts/pgr/docs/dist.md) "需为 32 的倍数"；
  但 [execute](file:///home/wangq/Scripts/pgr/src/cmd_pgr/dist/hv.rs#L102)
  仅校验 `opt_dim > 0`，`hash_hv_bit`/`hash_hv_sparse` 均支持任意维。
- 实测：`--dim 100`（非 32 倍数）正常运行。
- 建议：把"应为 32 的倍数"改为"建议为 32 的倍数（性能/对齐）"或在实现中
  校验，避免文档宣称硬性约束而实际未实施。

### 测试覆盖缺口

#### D7. 无 `dist pgi` 的 CLI 集成测试；`dist hv` 无空输入/NaN 覆盖
- 经核对：[cli_pgi.rs](file:///home/wangq/Scripts/pgr/tests/cli_pgi.rs) 已含
  `dist pgi` 的 CLI 测试（`command_pgi_dist_identical_and_disjoint`、
  `command_pgi_dist_param_mismatch_fails`），D7 关于"无 pgi CLI 测试"的表述
  不准确；真实缺口是**空索引**（D3）与 **hv 空输入**（D2）无回归覆盖。
- 建议：补 `dist pgi` 空索引回归、`dist hv` 空输入回归。

## 修复与验证（第 1 轮，2026-08-08）

全部缺陷均已修复并补回归测试：

- **D1**：`dist pgi` containment 改回**有方向**——分母用第一个索引
  `total1`（与命令族一致），见
  [dist.rs](file:///home/wangq/Scripts/pgr/src/libs/pgi/dist.rs#L81-L93)。
  新增回归 `command_pgi_dist_containment_directional`。
- **D2**：`dist hv` 空输入不再 NaN——`calc_distances` 与 `run_hv_files` 对
  `union==0`（空==空）给 jaccard=1、mash=0，空首集合 containment=0；`.hv`
  路径额外处理余弦 0/0（双空视为 sim=1）。新增回归
  `command_dist_hv_empty_inputs_no_nan`。
- **D3**：`dist pgi` 空==空不再 jaccard=0/距离 1，改为 jaccard=1、mash=0。
  新增回归 `command_pgi_dist_empty_indexes`。
- **D4**：`mash_distance` 现 clamp 到 [0,1]（与 Mash `min(1,dist)` 及自身
  注释一致）。新增单测 `test_mash_distance_bounded_to_one`。
- **D5**：修正 `args.rs` `sampler_arg` 注释与 `hv.rs` after_help——明确
  `--hasher` 仅对 minimizer/FracMinHash 生效，syncmer（DNA 与蛋白）忽略之
  （DNA 2-bit 滚动哈希、蛋白 RapidHash）。
- **D6**：`--dim` 帮助与 `docs/dist.md` 把"需为 32 的倍数"改为"建议为 32
  的倍数（对齐/性能），实现不强制"。
- **D7**：更正为真实缺口——新增 `dist pgi` 空索引与 `dist hv` 空输入回归
  （见上）。

验证：
- `cargo build` 通过；`cargo clippy --all-targets` 零新增告警；
  `cargo fmt --check` 干净。
- `cargo test --test cli_dist`（9）/ `cli_pgi`（10，含 2 新）/
  `cli_dist_mini`（5）/ `cli_dist_mash`（3）/ `cli_dist_frac`（3）全通过。
- `cargo test --lib hash::`（含新单测）/ `pgi::dist`（4）全通过。
- D4 clamp 仅影响极小 Jaccard，`dist mash` 字节级兼容测试（k=21/s=1000）
  不受影响（实测通过）。

## 第 2 轮（2026-08-08）

对 `dist hv` 的 FASTA 采样路径与 `.hv` 文件路径、`calc_distances` 数值语义、
`hash_hv_bit` 的 AVX2/标量尾部一致性、以及共享参数 `resolve_kmer_window` 做
第二轮深审，发现并修复如下问题。

### 代码 / 行为缺陷

#### D11. `dist hv` FASTA 路径对无关系列输出负的 Jaccard/Containment
- 位置：[hv.rs](file:///home/wangq/Scripts/pgr/src/libs/hv.rs#L394-L404)
  `calc_distances` 中 `inter` 直接取 `hv_dot(s1,s2)`，未 clamp 到 ≥0。
- `inter` 是点积，是**零均值**的共享 k-mer 数噪声估计；对完全无关系列，
  它可随机落在负值。`union = card1+card2-inter` 随之偏大，`jaccard =
  inter/union` 与 `containment = inter/card1` 为**负**（数学上无效，违反
  "Zero-Panic / 不输出无效值"原则）。
- 实测（15 对 100kb 随机无关系列，D=4096）：约 1/3 出现负值，如
  `jaccard -0.0089`、`containment -0.0174`、`inter 0`（负 f32 转 usize 饱和）
  、`mash 1.0000`（`mash_distance` 对负 jaccard 得 NaN，经 `min(1.0)` 兜底为
  1）。
- 对比：草图家族（mini/mash/frac）对无关系列正确输出 jaccard=0、
  containment=0、mash=1。`.hv` 文件路径因 `inter` 经 `usize` 转（饱和为 0）
  已天然 clamp，故只影响 FASTA 采样路径。
- 修复：`inter = hv_dot(s1,s2).max(0.0).min(card1).min(card2)`。新增单测
  `test_calc_distances_disjoint_no_negative`（50 次随机无关系列对，断言
  jaccard/containment ≥0 且非 NaN）。实测修复后无关系列输出
  jaccard=0、containment=0、mash=1。

#### D12. `hash_hv_bit` 测试参考 `hash_hv_bit_serial` 未随 D10 尾部修复更新
- 位置：[hv.rs](file:///home/wangq/Scripts/pgr/src/libs/hv.rs#L501-L532)
  （测试模块）。D10（第 1 轮）把生产路径 `hash_hv_bit` 与
  `hash_hv_bit_avx2` 的尾部维度改为接收每个种子的 ±1 贡献（仅末尾保留
  −N 常量已不再正确）；但作为对比基准的 `hash_hv_bit_serial` 仍只处理
  `hv_d/64` 个完整 chunk，把尾部留在 −N。
- 影响：`test_hash_hv_bit_avx2_serial_vs_simd`（及任何含尾部的维数）在
  非 64 倍数维（如 1056/1064）上失败——AVX2 输出尾部是正确值、参考基准是
  旧行为，二者不一致。该测试注释"1056 = 33×32 + 0 tail；tail keeps −N"
  亦过时（1056/64=16 余 32，实为 32 维尾部）。
- 修复：为 `hash_hv_bit_serial` 补齐与生产路径一致的尾部处理，并更新注释。
  修复后所有 `hv::` 单测通过（20 个），含 1056/1064 尾部维。

### 次要观察（未修复，符合"简洁优先"）

- **D13（观察）**：`seq_mins` 的 `mod` hasher 路径用 `opt_window as u16` 传给
  `minimizer_iter::MinimizerBuilder::width`（[hash.rs](file:///home/wangq/Scripts/pgr/src/libs/hash.rs#L110-L118)）。
  若 `--window > 65535`（仅能通过显式超大 `-w` 触发，正常使用不可达），会
  静默截断为 `w & 0xffff` 而不报错。属理论性边界；按"不为不可能发生的场景
  写错误处理"原则不做改动，仅记录。

### 验证（第 2 轮）
- `cargo build` 通过；`cargo clippy --all-targets` 对本次改动零新增告警；
  `rustfmt --check src/libs/hv.rs` 干净（`benches/hv_ann_recall.rs` 的 fmt diff
  为无关的既有差异，未改动）。
- `cargo test --lib hv::`：20 个全通过（含新 `test_calc_distances_disjoint_no_negative`
  与修复后的尾部 parity 测试）。
- `cargo test --test cli_dist`（11）/ `cli_dist_mini`（5）/ `cli_dist_mash`（3）/
  `cli_dist_frac`（3）全通过。
- 实测：`dist hv` 对 15 对随机无关系列不再输出负 jaccard/containment。

## 第 3 轮（2026-08-08）

对 `dist` 各子命令的 `--hasher mod` 路径与参数校验做第三轮深审，发现并修复
一个 Zero-Panic 违规。

#### D14. `--hasher mod` 配偶数窗口触发 minimizer_iter 断言 panic
- 位置：[hash.rs](file:///home/wangq/Scripts/pgr/src/libs/hash.rs#L110-L126)
  `seq_mins` 的 `mod` 分支用 `minimizer_iter::MinimizerBuilder::new_mod()`
  `.width(opt_window as u16)`。`minimizer_iter` 的 builder 在 `iter()` 时
  **断言窗口必须为奇数**（`builder.rs:183 assert width must be odd`）。
- 影响：`pgr dist mini --hasher mod -w <偶数>` 与
  `pgr dist hv --hasher mod -w <偶数>` 直接 **panic**（`thread 'main'
  panicked ... assertion left == right failed: width must be odd`），违反
  "Zero-Panic"原则。默认窗口（mini w=5、hv w=1）均为奇数所以默认不触发；
  仅显式传偶数 `-w`（如 `-w 4`，或大窗口截断为偶数）触发。
- 实测：`dist mini a.fa --hasher mod -w 4 -k 21` → panic；`-w 5` 正常。
  同样触发 `dist hv` 的 minimizer 路径。
- 修复：在 `seq_mins` 的 `mod` 分支前校验 `opt_window % 2 == 1`，奇数窗口
  才调用 minimizer_iter，偶数窗口返回友好错误
  `--hasher mod requires an odd window size`。同时覆盖 `dist mini` 与
  `dist hv` 两个消费者。新增回归
  `command_dist_mini_mod_rejects_even_window`。
- 说明：`seq_sketch`（同用 minimizer_iter + `w as u16`）同样有该约束，但仅被
  自身测试引用、无命令调用，不在 dist 可达路径内。

### 验证（第 3 轮）
- `cargo build` 通过；`rustfmt --check` 对改动文件干净。
- `cargo test --test cli_dist_mini`：6 个全通过（含新
  `command_dist_mini_mod_rejects_even_window`）。
- 实测：`dist mini/hv --hasher mod -w 4` 均返回友好错误而非 panic；
  `-w 5` 正常输出。

## 第 4 轮（2026-08-08，收敛）

对 `dist` 全部子命令（mini/mash/frac/hv/pgi）的 CLI、共享参数
（`resolve_kmer_window`、`--parallel` 1..=1024 校验）、输出列格式、`docs/dist.md`
文档、以及 `libs`（hash / hv / syncmer / pgi::dist / pgi::to_hv / par）做第四轮
复核，未发现新的代码 / 行为 / 文档缺陷：

- 所有 `unwrap()` 均作用于 clap 具默认值或 `required` 的参数（`get_one`），
  不构成用户输入 panic。
- `spawn_writer_and_pool` 的 writer 线程 `unwrap` 由 `join().map_err` 捕获并
  转为友好错误（进程不 abort）。
- `.hv` 文件路径、`dist pgi` 归并、空输入、无关系列（D11 修复后）等数值
  语义均与草图命令族一致。
- 文档与修复后的行为一致（pgi containment 有方向、hv 空输入 jaccard=1、
  `--hasher` 仅 minimizer 生效、`--dim` 非强制 32 倍数）。

验证（第 4 轮）：`cargo build` 通过；`cargo clippy --all-targets` 零新增告警；
`cargo test --test cli_dist`（11）/ `cli_dist_mini`（6）/ `cli_dist_mash`（3）/
`cli_dist_frac`（3）/ `cli_pgi`（10）与 `--lib hash::`（14）/ `hv::`（20）/
`pgi::dist`（4）全部通过。

**结论**：`dist` 命令族在第 1 轮（D1–D10）+ 第 2 轮（D11、D12，D13 观察）+
第 3 轮（D14）之后，第 4 轮已收敛、无新问题。

