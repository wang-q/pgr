# multiz 源码分析

> 初稿整理于 2026-02；2026-08-03 依据 `multiz-multiz/` 目录（v11.2）源码
> 逐文件重读补全（multiz.c / mz_preyama.c / mz_yama.c / mz_scores.c），
> 并记录 pgr 的直译进展。2026-08-11 再核对全部源码，修正 GAP 表第 4 条
> 映射、补充 maf.c/multi_util.c 分析、quasi-natural gap 逐对累加语义及两处
> 源码怪癖。设计/实现记录见 [[fas-multiz.md]]。

本文档分析 UCSC `multiz-tba` 软件包的核心组件 `multiz`。源码在仓库
`multiz-multiz/` 目录（2016 版 v11.2），环境已安装二进制
`/home/wangq/.cbp/bin/multiz`；自带测试数据 `multiz-multiz/test/`
（S288cvsRM11_1a.maf + S288cvsSpar.maf，与 pgr `tests/fas/*.slice.fas`
同源）。

## 1. multiz 概览

### 1.1 工具与输入输出

*   **功能**: 对两个多序列比对文件（MAF 格式）进行比对合并，生成新的多序列比对。
*   **核心假设**: 两个输入 MAF 文件的 Top Row 是同一参考序列；block 按参考
    坐标升序、参考上单覆盖（single-coverage）。

**命令行接口**：

```bash
multiz [R=?] [M=?] [L=?] [S=?] file1 file2 v [out1 out2] [nohead] [all]
```

*   `R=30`: DP 半径（香肠带扩展半径）。
*   `M=1`: 最小输出宽度（MIN_OUTPUT_WID）。
*   `L`/`S`: 大/小块断点宽度（LRG_BREAK_WID / SML_BREAK_WID）。
    *   **源码怪癖**：`main` 确实解析这两个参数（`multiz.c:189-196`），`multi_util.c:13-14`
        也定义了全局变量（`LRG_BREAK_WID=20`、`SML_BREAK_WID=2`），但**全仓库再无任何地方
        引用它们**——`multiz()` 与 `pre_yama()` 均不使用。即 `L=`/`S=` 是 v11.2 里的**死参数**，
        解析后即丢弃（`multi_util.c` 中 `OVERLAP_THRESHOLD`/`MIN_CHAIN`/`MIN_CLUSTER_CHAIN`/
        `OVERLAP_LEN_THREH`/`MIN_DISTANCE`/`MIN_SPB` 同理，均为未使用的遗留全局）。pgr 直译时
        不要为这两个参数留接口。
*   `v`: 参考固定模式。
    *   `0`: 两个参考都可微调（需要第二次 yama 对齐参考行）。
    *   `1`: 第一个文件的参考固定，第二个文件的参考可相对滑动。
*   `out1 out2`（可选）: 收集两个输入中未被使用的 block；缺省时合并结果与
    unused 都写 stdout。
*   `nohead`: 不输出 MAF 头；`all`: 输出单行块。
    *   **源码怪癖（已实测）**：是否输出单行块由全局 `row2` 控制，条件为
        `row2 == 0 || components->next != NULL`（`row2==0` 时连单行块也输出，
        见 `multiz.c:69,75,255,260`）。`multi_util.c:22` 中 `row2` **默认就是
        0**，而 `main` 对 `all` 的处理是 `row2 = 0`（`multiz.c:203-206`）
        ——故 `all` 在 v11.2 里**实为空操作**，默认即输出单行块。实测：同一
        输入在"无 `all`"与"有 `all`"两种调用下，单行 unused 块都照常输出。
        这与 usage 文本声称的 "null: not to output single-row blocks" 相矛盾，
        属源码文档/实现不一致（疑似 bug：若要"默认不输出单行块"，条件应反向
        为 `row2 != 0` 之类）。pgr 直译时应以"默认输出单行块"为准，不为 `all`
        单设语义。`fas-multiz.md` §2.6 沿用的"默认不输出"表述同样需按此更正。

### 1.2 核心算法一句话

"参考锚定下的 profile–profile 动态规划"：把两个 MAF 视作列剖面（各自
物种集合），以共同参考建立 DP 网格坐标与边界（LB/RB），用三状态
（C 替换 / D 删除 / I 插入）的 sum-of-pairs + 准自然 gap 成本做全局比对。

## 2. 源码逐文件分析

### 2.1 multiz.c —— 主流程

`main` 按参考染色体分组（`seperate_cp_wk`），逐染色体调用 `multiz()`：

1. 维护两个 block 流指针 a1/a2，按参考坐标滑动；
2. 未覆盖的前端 block 输出到 fpw1/fpw2（unused 收集，`MIN_OUTPUT_WID`
   过滤，`row2==0` 或非单行块才输出）；
3. 重叠区 `[beg,end]` 交给 `pre_yama` 合并；
4. `keep_ali` 切分被消耗的后半部分，`print_part_ali_col` 输出块前后
   未被合并的列；
5. 合并结果（`mafWrite(stdout, new_ali)`）与 unused 分离输出。

关键点：multiz 是**逐重叠区的块流合并**，输出多块；单覆盖假设下每个参考
位置只被合并一次。

**源码怪癖**：`multiz()` 内残留一处调试代码
`if (a1->components->start == 11305) test++;`（multiz.c:89-90），
`test` 变量从不被使用，纯死代码，可忽略。

### 2.2 mz_preyama.c —— 预处理与边界

`pre_yama(a1, a2, beg, end, radius, v)`：

1. **剖面**：K = a1 组件数（v=0 时减 1 排除参考），L = a2 的非参考组件数；
   A 剖面 = a1 的 K 行（v=1 含参考，参考"将被固定"），B 剖面 = a2 的非
   参考行。**B 剖面永远不含参考**——参考只作为坐标/边界锚，不重复参与
   打分。
2. **rmColDash**：去掉参考行以下全 dash 的列，返回旧列→新列映射（map1/map2）。
3. **LB/RB 建立**：两个参考去 dash 后逐位配对（第 k 个碱基 ↔ 第 k 个
   碱基），为每个 ref_a 碱基列建立点约束 `LB[i]=RB[i]=j`；参考行有 gap
   的列不设约束，由单调化继承。
4. **smooth**（见 §2.5）：单调化 + 半径香肠，端点归位（LB[0]=0, RB[M]=N）。
5. **yama**：带内 C/D/I DP。
6. **v=1**：`mafBuild(AL_new, K+L, ...)` 直接重建（参考行 = A 剖面的 a1
   参考，经 DP 列映射放置）。
7. **v=0**：先做一次不含参考的 yama（A=K 行）得 AL_new（K+L 行）；再做
   第二次 yama 精修参考行——A2 = a1 参考行（1 行）、B2 = AL_new（K+L 行）。
   **两套 LB/RB 并行建立**：一套由 A2↔AL_new（经 `mapping(A→AL_new)`，见
   `mz_preyama.c:276`），另一套由 B2 = a2 参考行↔AL_new（经
   `mapping(B→AL_new)`，`mz_preyama.c:308`），随后 `LB=MIN(LB,LB2)`、
   `RB=MAX(RB,RB2)` 合并（`M3 != N3` 则 fatal），最后
   `yama(A2,1,M3,AL_new,K+L,...)` → `mafBuild(K+L+1)`。即 a1 的参考被对齐进
   结果，a2 的参考只通过边界约束参与——这是 pgr 若想实现"参考也参与精修"可
   借鉴的模式（见 §3.3）。

    **源码怪癖**：v=0 分支的 `mapping(A,1,K,1,M,...)`（mz_preyama.c:276）以
    `a_row2=K` 传入，而此刻 A 在 `K--` 后仅含 0..K-1 共 K 行，读 `A[i][K]`
    属越界（应传 `K-1`）；同函数的 B 侧调用 `mapping(B,0,L-1,...)`
    （mz_preyama.c:308）传 `0..L-1` 则正确。属源码遗留问题，pgr 若直译 v=0
    分支应规避（pgr 目前未实现 v=0，见 §3.2）。

`mafBuild`：从 (K+L)×M_new 列矩阵重建 MAF 块——每行按参考坐标重算
start/size，`nc->size == 0`（全 gap 行）丢弃，`score = mafScoreRange`。

### 2.3 mz_yama.c —— C/D/I 三状态 DP

网格 M×N（M = A 剖面列数，N = B 剖面列数），每格三个状态 + 一个 trace
字节（`flag_c | flag_d<<2 | flag_i<<4`，FLAG_C=0 / FLAG_I=1 / FLAG_D=2）：

*   **C（替换）**：A 列与 B 列对齐。候选来自 (row-1, col-1) 的 C/D/I
    （diag 保存的上行值）；列分 = **K×L 全物种对** `Σ SS(A[row][i],
    B[col][j])`（含 A 侧的参考行）。
*   **D（删除）**：A 列插入、B 全 dash。候选来自 (row-1, col) 的 C/D/I；
    成本 `- n*L*gap_extend`（n = A[row] 非 dash 数）。
*   **I（插入）**：B 列插入、A 全 dash。候选来自 (row, col-1) 的 C/D/I；
    成本 `- n*K*gap_extend`（n = B[col] 非 dash 数）。

**准自然 gap 修正**：对每条候选路径（上一步为 C/D/I），按"最后两条边"
的 A/B dash 模式查 `GAP(s,t,u,v)` 表（见 §2.4），受 LB/RB 条件约束
（如 C 的 x 修正要求 `row>1 && col>LB[row-2]+1`）。注意此修正按物种对
(i,j) 双重循环**逐对累减**（见 §2.4 的重要语义），并非整列一次性查表。

**端部 gap 免费**：I 的修正只在 `row<M`（末行不收 open），C 只在
`col>1`（起点不收 open），D 只在 `0<col<N`（起点/终点不收 open）；
extend 成本恒收。

**初始化**：`dp[0].C=dp[0].D=dp[0].I=0`；行 0 只有 I（
`I[col]=I[col-1]-n*K*gap_extend`）；每行 `col=LB[row]` 处 I=MININT，
`col=LB[row-1]` 处 C=MININT。

**回溯**：从 (M,N) 取 max(C,D,I)，沿 trace 字节回放 edit script
（C/D/I 序列），再正放成输出矩阵：C→`new_col(A列+B列)`，I→`dashes+B列`，
D→`A列+dashes`。

### 2.4 mz_scores.c —— 打分与准自然 gap 表

*   **HOX70**（human-rodent，默认）：
    ```
    A   C    G    T
    A  91 -114  -31 -123
    C -114 100 -125  -31
    G  -31 -125 100 -114
    T -123  -31 -114  91
    ```
    gap open 400 / extend 30；`SS('-',x) = -gap_extend`，`SS('-','-') = 0`，
    其余未知字符 -100。
*   **GAP(s,t,u,v) 表**（s/t = 前一条边的 A/B 侧是否 dash，u/v = 当前边的
    A/B 侧是否 dash，1=dash）：
    16 种构型中 6 种收 gap_open：
    ```
    GAP(0,0,0,1) GAP(0,0,1,0) GAP(0,1,1,0)
    GAP(1,0,0,1) GAP(1,1,0,1) GAP(1,1,1,0)
    ```
    即边对（前边/后边）：(xx/x-)、(xx/-x)、(x-/ -x)、(-x/x-)、(--/x-)、
    (--/-x)。
    *   **重要语义**：准自然 gap 罚分在 yama 中对**每一对物种 (i,j) 逐对累加**
        ——内层双重循环 `for i in K: for j in L:` 里查 `GAP(s,t,u,v)` 并累减到
        候选值。因此 gap_open 不是对整列收一次，而是对"该处实际存在 gap 的
        每个物种对"各收一次；某列只有部分物种有 gap 时，罚分按有 gap 的物种数
        成比例。这正是 Altschul quasi-natural gap 与"固定 gap_open"的关键差异
        （后者对所有物种对一视同仁）。
*   **mafScoreRange**：MAF 块 SP 打分（列内所有物种对 SS + 相邻列 GAP2
    修正），供 `mafBuild` 写 `a score=`。

### 2.5 band（LB/RB）处理

multiz 的带不是固定带宽，而是**参考锚定的动态边界数组**：

1. **点约束**：两个参考去 dash 逐位配对，每行钉到对应列（LB=RB=配对列）。
2. **smooth()**：
    ```c
    // 单调化：LB 前向取 max（非递减），RB 后向取 min（非递增）
    for (i = j = 0; i <= M; ++i) LB[i] = j = MAX(j, LB[i]);
    for (i = M, j = N; i >= 0; --i) RB[i] = j = MIN(j, RB[i]);
    // 香肠：以 radi = MIN(M, radius) 扩展，端点强制归位
    for (i = M; i > radi; --i) LB[i] = MIN(MAX(LB[i] - radi, 0), LB[i - radi]);
    for (; i >= 0; --i) LB[i] = 0;
    for (i = 0; i < M - radi; ++i) RB[i] = MAX(MIN(RB[i] + radi, N), RB[i + radi]);
    for (; i <= M; ++i) RB[i] = N;
    ```
    带宽 = 参考锚定对角线 ± R。不变量：`LB[0]==0`、`RB[M]==N`、LB/RB 单调。
3. **yama 每行只在 [LB[row], RB[row]] 内算**（`col = LB[row]-1;
   while (++col <= RB[row])`），存储为逐行变长行（`tback_row[row] =
   tbp - LB[row]`），不是固定宽度矩阵；准自然 gap 修正条件也用 LB 保证
   带内有效。`yama` 开头断言带宽（`RB[row]-LB[row] >= MIN(N,10)`）与
   单调性，不满足直接 fatal。

### 2.6 maf.c / multi_util.c / util.c —— MAF I/O 与共享工具

**`maf.c`** 是迷你版 MAF 读写（header 注明 "version 12"，"Stolen from
Jim Kent & seriously abused"），依赖 zlib，`.gz` 输入由 `gzopen`/`gzgets`
透明支持：

*   **mafOpen**: 校验首行 `##maf version=%d`（失败 fatal），解析 `scoring=`
    字段。
*   **mafNext**: 跳过 `#`/空行，块必须以 `a` 行开始；逐 `s` 行解析
    `s src start size strand srcSize text`，做三重 sanity check：各组件
    text 长度一致、`start+size <= srcSize`、text 中非 dash 计数 == 声称的
    `size`；任一不符直接 fatal。注意 multiz 对坏数据是**硬失败（abort）**，
    pgr 的 "Zero Panic" 约定需转为友好报错。
*   **parseScoreLine**: 解析 `a` 行 `score=`、`amplifier=<row>`、`copy=<row>`
    ——把对应行 `paralog` 标为 `'a'`/`'c'`（默认 `'s'` 单例）。
*   **mafWrite**: 逐字段对齐列宽（src/start/size/srcSize 各取块内最长），
    并按各组件 paralog 输出 `amplifier=`/`copy=`。
*   **mafColDashRm**: 原位删除全 dash 列（供 make_part_ali_col 收尾）。

**`multi_util.c`** 提供 multiz 全程共享的工具：

*   **keep_ali(ali, beg)**: 从参考位置 beg 起保留 block 剩余部分（丢弃完全
    落在 beg 之前的组件），重算各组件 start/size 与块 score。
*   **make_part_ali_col / print_part_ali_col**: 按列区间 [cbeg,cend] 切出子块
    （全 gap 行剔除），供 multiz 输出块前后未被合并的列。
*   **mafPos2Col**: 参考位置 → 列号（遍历 text 跳过 dash）。
*   **seperate_cp_wk**: 把 cp_list 中 `src==chr` 的块移到 wk_list（按参考
    染色体分组，供 main 逐染色体调用 multiz）。
*   **retrieve_first**: 弹出链表头。
*   **parseSrcName / parseSrcName2**: 把 `src`（形如 `species.contig`）拆为
    name 与 contig。

**`util.c`**: `fatal`/`fatalf`（打印 argv0 前缀后 exit(1)）、`fatalfr`（追加
strerror）、`ckalloc`（0 长度按 1 分配、失败 fatal）、`copy_string`，全局
`argv0`。`util.h` 另定义了 `MAX`/`MIN` 宏与 `uchar = unsigned char`。

## 3. 对 pgr 的启示与落地

### 3.1 已落地：yama 引擎直译（2026-08-03）

pgr `libs::fas_multiz` 的 `banded_align.rs` 已按 §2.3–§2.5 直译：

*   C/D/I 三状态 + 准自然 GAP 查表 + 端部 gap 免费；
*   **全物种对评分（K×L）**：A 剖面含参考、B 剖面含参考（与 multiz 的
   "B 不含参考"略有差异——pgr 保留 ref-ref 对作为锚，LB/RB 与 ref-ref
   信号互为补充）；
*   LB/RB 逐位配对 + smooth 香肠 + 逐行变长存储；
*   边缘裁剪改为"全物种 gap 列"判定（multiz 的 rmColDash 语义在输出侧
   的对应）。

实证：无 LB/RB 时自由端 gap 会把列差整段堆到块端、删掉真实内容（Spar
4057→3724）；LB/RB 落地后 S288c 三输入合并各物种碱基数与输入逐碱基
一致（ref 3826 / RM 3834 / Spar 4057 / YJM 3822）。细节见
[[fas-multiz.md]] §4.1–§4.2。

### 3.2 已知差异（不追齐）

*   **打分矩阵**：pgr 的 `hoxd55` 与 multiz HOX70 数值完全相同（见
    §2.4），直接用 `hoxd55`，未引入 `hox70` 别名；
*   **v=0 模式**未实现（pgr 渐进合并天然等价 v=1，v=0 价值不明）；
*   **输入/分块语义**：multiz 是 MAF 块流、逐重叠区合并；pgr 是 fas
    窗口合并。多块输入中"仅共享参考且参考去 gap 不等"的块对在 pgr 的
    crossover 路径会被拒绝（需共享非参考物种打分），为已知限制；
*   **gap 参数**：multiz 固定 400/30，pgr 已硬编码同样的 400/30（不提供
    打分 CLI 参数），与 multiz 一致。

## 4. pgr MAF 实现现状

目前 `pgr` 已在 `src/libs/fmt/` 下实现了完整的 MAF 读写支持。

### 4.1 模块分布与功能

*   **读写统一**：位于 `src/libs/fmt/maf.rs`
    *   **读取 (Reader)**:
        *   **核心函数**: `next_maf_block`, `parse_maf_block`。
        *   **特性**: 支持 `a` (alignment) 和 `s` (sequence) 行解析，`a`
            行 `score=` 字段已解析到 `MafAli.score`。
        *   **坐标转换**: `MafComp::to_range()` 将 MAF 的 0-based 坐标
            转换为 1-based inclusive 格式（如 `chr:start-end`）。
        *   **负链处理**: 自动处理负链坐标，将其转换为相对于正链的坐标范围。
    *   **写入 (Writer)**:
        *   **核心结构**: `MafWriter`。
        *   **特性**: 支持输出标准 MAF 头信息 (`##maf`) 和对齐块，自动
            处理列宽对齐。

### 4.2 数据结构

读写使用同一套结构体：

| 结构体 | 用途 | 关键字段 |
| :--- | :--- | :--- |
| `MafComp` | `s` 行（一条序列组件） | `src`、`start`、`size`、`strand`、`src_size`、`text`（均为 `String`/`usize`） |
| `MafAli` | `a` 行 + block 内所有 `MafComp` | `score: Option<f64>`、`components: Vec<MafComp>` |

### 4.3 MAF 格式实现对比 (pgr vs multiz vs UCSC)

1.  **代码同源性**:
    *   `multiz` 的 `maf.c`（header 注明 "version 12"）标注了 "Stolen
        from Jim Kent & seriously abused"——精简版 mini-maf，移除了
        `linefile.h`/`common.h` 等依赖，直接使用标准 C 库函数。
    *   UCSC `chainnet` 中的 `maf.c` 是完整版，依赖 Kent Source 基础
        设施库。
2.  **功能差异**:
    *   UCSC 完整版支持 `i`/`q`/`e`/`r` 扩展行、复杂内存管理、大量辅助
        函数（`mafSubset`/`mafFlipStrand`/`mafScoreMultiz` 等）。
    *   multiz 精简版仅保留核心 `a`/`s` 行解析，足以支撑比对算法。
    *   pgr 解析专注于核心行并健壮跳过未知行；坐标系统与 UCSC/multiz
        一致（0-based start, 1-based size），内部提供 1-based inclusive
        转换；Rust 所有权模型替代 C 手动内存管理。

### 4.4 总结

`pgr` 的 MAF 模块已成熟（读、写、坐标标准化齐全），配合 `pgr maf
to-fas` 可与 multiz 测试数据（`multiz-multiz/test/*.maf`）互通；
`fas-multiz` 的 yama 直译则把 multiz 的算法核心完整落地。
