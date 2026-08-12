# anchr / ovlpr / App-Dazz 工具链盘点与 pgr 边界划分

> 2026-08-12。核对 `~/Scripts/anchr`、`~/Scripts/intspan`（ovlpr 前身）、
> `~/Scripts/App-Dazz`（Perl 前身）三者的功能，对照 pgr 现有命令，给出
> "anchr 留在外部、pgr 收通用原语"的边界划分方案。

## 1. 工具沿革

| 工具 | 语言 | 时期 | 关系 |
| :--- | :--- | :--- | :--- |
| `App::Dazz`（`dazz`） | Perl | 2017 | Daligner-based UniTig utils，最早形态 |
| `ovlpr` | Rust | 2020（intspan 时代） | 独立二进制；后整体移入 anchr |
| `anchr` | Rust | 现行 | "Assembler of N-free CHRomosomes"，ovlpr 命令并入其子命令 |

沿革链：`App-Dazz（Perl）→ ovlpr（intspan）→ anchr 子命令`。`ovlpr` 已无
独立二进制残留（`~/bin/ovlpr` 不存在），intspan 源码已无 ovlpr 代码。

## 2. anchr 定位与命令清单

anchr = **染色体级组装流程编排器**：0–9 编号 tera 模板（reads 处理 → 比对 →
anchors → 组装 → 统计）+ 底层命令。当前外部依赖：bbtools（bbduk/tadpole/
bbmap）、DALIGNER/App-Dazz（7 步 glue/fill）、hnsm、rgr、quorum 等。

命令（18 个，`src/cmd/`）：

| 命令 | 功能 | 类别 |
| :--- | :--- | :--- |
| `template` | 生成 0–9 流程 Bash 脚本 | 流程编排 |
| `dep` / `ena` | 依赖安装检查 / ENA 下载脚本 | 流程胶水 |
| `anchors` | 从 contigs 选 anchors（proper covered regions） | 特有策略 |
| `trim` / `mergeread` / `quorum` / `unitigs` | reads 处理（bbtools/quorum wrapper） | 原语（pgr 已覆盖） |
| `overlap` / `overlap2` | DALIGNER overlap 检测 | 原语（pgr 已覆盖） |
| `contained` / `orient` / `merge` | overlap 图：去包含 / 定向 / 合并 | 原语（pgr 已覆盖） |
| `covered` | 从 PAF/ovlp 算覆盖区间（--coverage/--base/--mean/--longest） | 原语（pgr 已覆盖） |
| `dazzname` / `show2ovlp` / `paf2ovlp` / `restrict` | DALIGNER 生态格式转换/过滤 | DALIGNER 特定 |
| `mergeread` | bbtools merge wrapper | 原语（pgr 已覆盖） |

App-Dazz 独有命令（anchr 无）：`cover`（第一文件被第二文件覆盖的可信
区间）、`group`（anchors 按长读分组）、`layout`（组内 anchors 布局 →
contig.fasta）。anchr 的 `7_glue_anchors` / `7_fill_anchors` 模板至今调用
`dazz group` / `dazz layout`（依赖 Perl App-Dazz + dazz_db 数据库）。

## 3. pgr 覆盖对照

> 2026-08-13 更新：`fq`/`asm` 业务已随迁移回到 anchr（anchr 自实现，
> 依赖 pgr 基础层）。下表为迁移后的分工。

| anchr/App-Dazz 命令 | 功能 | 迁移后归属 |
| :--- | :--- | :--- |
| `unitigs` / `overlap` / `overlap2` / `contained` / `orient` / `merge` | reads/unitig 组装与 overlap 处理 | anchr `asm` 命令组自实现 |
| `trim` / `mergeread` | reads 处理 | anchr `fq` 命令组自实现 |
| App-Dazz `cover` / `group` / `layout` | 覆盖区间 / 长读分组 / 组内布局 | anchr 自实现（分层策略，§5） |
| `covered` | 覆盖区间 | pgr `paf coverage`（无 cg:Z 已支持）+ `rg coverage` |
| `covered --mean/--longest` | 多文件均值/最长区间 | pgr 需组合（paf coverage TSV + awk/runlist） |
| `restrict` | 过滤到已知对 | pgr `paf query --subset-sequence-list`（语义不同） |
| `dazzname` | DALIGNER 式重命名 | 无兼容需求（DALIGNER 生态退役） |
| `paf2ovlp` / `show2ovlp` | DALIGNER ovlp 格式转换 | ❌ DALIGNER 生态特定 |

pgr 侧提供的基础：FASTA/FASTQ 读入（`libs/fmt`）、Phred 编码（`fq::qual`）、
配对 FASTQ 读取（`fq::pairs`）、k-mer、PAF、io/ds/loc/sys（anchr 依赖
pgr crate，见 [[design/fq-asm-migrate.md]]）。

## 4. 边界划分方案（2026-08-12 用户裁定方向：anchr 留在外部）

pgr 保持"通用基因组数据处理工具集"定位，不吸收 anchr 的流程编排逻辑。
划分原则：**anchr 留"流程 + 策略"，pgr 收"通用原语"，DALIGNER 生态退役**。

### 4.1 留在 anchr（约 6 项 + 模板体系）

- `template`（0–9 模板）：anchr 灵魂；pgr 明确不做 pl 大流程（原语路线，
  见 todo §3 "anchr 模板替换"）；
- `anchors`：anchor 选择策略（流程特有）；
- `quorum` / `mergeread` / `ena` / `dep`：外部工具集成 + 数据获取（流程胶水）；
- 0/1/3/8/9 步模板：reads 质控、bwa 比对、megahit/spades 组装、quast/busco
  统计——纯编排。

### 4.2 pgr 提供的基础层（fq/asm 已迁回 anchr）

- **2026-08-12/13**：`fq`/`asm` 业务逻辑迁回 anchr（reads 处理 + 组装
  归位组装器，双轨 golden 核对后从 pgr 删除，见 [[design/fq-asm-migrate.md]]）；
  pgr 保留基础层（FASTA/FASTQ 读入、Phred 编码、k-mer、PAF、io/ds/loc/sys），
  anchr 依赖 pgr crate。pgr 仍提供的通用原语：`paf coverage`、`paf graph`
  （分组/图）、`rg coverage`；anchr 自实现 `sam`（ihist/to-rg，处理
  `asm map` 的 SAM 输出）。
- 7_glue/7_fill 的 `dazz group/layout` → anchr 自实现（`paf graph` 分组 +
  anchr `asm layout`，见 §5）。

### 4.3 退役（DALIGNER 生态，不迁）

- `dazzname` / `show2ovlp` / `paf2ovlp` / `restrict`；
- `overlap`/`overlap2` 的 DALIGNER 依赖（换 anchr `asm ovlp`）；
- 7_glue/7_fill 的 App-Dazz（Perl）+ dazz_db 依赖（换 anchr 自实现
  group/layout，overlap 用 anchr `asm ovlp`）。

## 5. 待办与风险

- **7_glue/7_fill 替换评估（2026-08-12 已做，结论：分层保留在 anchr）**：
  `dazz group`/`layout` 的语义（读 `App-Dazz/lib/App/Dazz/Command/
  group.pm`/`layout.pm` 源码）：
  * **group**：只保留 anchor×long 的 overlap → multi-matched 长读剔除
    （同一长读多次匹配同一 anchor）→ anchor 图建边（边需 ≥coverage 个
    长读支持 + 距离判断 + 链方向一致）→ 连通分量 = 分组；
  * **layout**：用 overlap 端点建**有向图**（g_strand=0 才保留）→
    transitive reduction → 路径分解（linear/branched/cyclic）→ 按
    relation 距离拼接**锚序列**成 contig（长读只提供证据不参与拼接）。
  **与 pgr 对照**：`paf graph`（DSU 传递闭包）是"序列段平等"图，无
  anchor/long 分层、无证据计数/距离/方向一致性过滤、无 multi-matched
  剔除；`asm layout` 是 unitig×unitig 的 OLC 链化，不消费"锚经长读链接"
  的分层证据。**直接替换不可行**（语义差异大）。
  **结论**：group/layout 的**分层策略逻辑留在 anchr**（用 Rust 重写
  替换 Perl App-Dazz，overlap 检测用 anchr `asm ovlp`/`olc` 甩掉
  DALIGNER），pgr 提供覆盖 `paf coverage`、图构建 `paf graph`；
  ovlp.tsv 13 列格式可保留为 anchr 内部格式（或换 PAF）。**pgr 侧无需
  新增命令**，符合"anchr 留流程 + 策略"的边界原则。
- **`paf coverage` 无 cg:Z gap**：ovlpr `covered` 用 PAF start/end（不依赖
  CIGAR），pgr `paf coverage` 要求 `cg:Z`（无标签记录不贡献）——补丁方向：
  无 cg:Z 时退回用 target 区间算覆盖。**2026-08-12 已落地**（见 todo §5）。
- **anchr 模板替换整体**：见 anchr `notes/todo.md`（trim.era.sh 等
  bbtools 调用换 anchr fq/asm 命令链，用户自处理）。
