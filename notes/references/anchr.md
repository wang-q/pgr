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

| anchr/App-Dazz 命令 | 功能 | pgr 对应 | 状态 |
| :--- | :--- | :--- | :--- |
| `unitigs` | 从 reads 建 unitigs | `asm unitig` | ✅ 已覆盖 |
| `overlap` / `overlap2` | overlap 检测 | `asm ovlp` / `asm olc` | ✅ 已覆盖（自包含，不依赖 DALIGNER） |
| `contained` | 丢弃被包含 unitig | `asm olc` 的 contain 预过滤 | ✅ 已覆盖 |
| `orient` / `merge` | overlap 图定向/合并 | `asm layout` + `asm cns` | ✅ 已覆盖 |
| `covered` | 覆盖区间 | `paf coverage`（需 cg:Z）+ `rg coverage` | ✅ 基本覆盖（见 §5 gap） |
| App-Dazz `cover` | 被第二文件覆盖的可信区间 | `paf coverage -m N` | ✅ 已覆盖 |
| App-Dazz `group` | anchors 按长读分组 | `paf graph`（DSU 连通分量） | ✅ 已覆盖 |
| App-Dazz `layout` | 组内 anchors 布局 | `asm layout` + `asm cns` | ✅ 已覆盖 |
| `trim` / `mergeread` | reads 处理 | `fq` 系列命令链 | ✅ 已覆盖（见 todo §3 模板替换） |
| `covered --mean/--longest` | 多文件均值/最长区间 | 需组合（paf coverage TSV + awk/runlist） | ⚠️ 部分 |
| `restrict` | 过滤到已知对 | `paf query --subset-sequence-list` | ⚠️ 语义不同，实现极简 |
| `dazzname` | DALIGNER 式重命名 | `fa filter --simplify`（格式不同） | ⚠️ 无兼容需求 |
| `paf2ovlp` / `show2ovlp` | DALIGNER ovlp 格式转换 | 无 | ❌ DALIGNER 生态特定 |

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

### 4.2 迁移到 pgr（通用原语，大部分已就位）

- 已覆盖（anchr 模板改调用即可）：`asm unitig/ovlp/layout/cns/olc`、
  `paf coverage`、`paf graph`、`fq` 系列；
- 待补小 gap：`paf coverage` 支持无 `cg:Z` 的 PAF（见 §5）；
- 7_glue/7_fill 的 `dazz group/layout` → `paf graph` + `asm layout`
  （需先对照语义，见 §5）。

### 4.3 退役（DALIGNER 生态，不迁）

- `dazzname` / `show2ovlp` / `paf2ovlp` / `restrict`；
- `overlap`/`overlap2` 的 DALIGNER 依赖（换 pgr `asm ovlp`）；
- 7_glue/7_fill 的 App-Dazz（Perl）+ dazz_db 依赖（换 pgr）。

## 5. 待办与风险

- **7_glue/7_fill 替换评估（唯一有算法语义差异的迁移）**：对照 `dazz
  group`（anchors 按长读分组规则）与 `paf graph`、`dazz layout`（组内布局
  算法）与 `asm layout` 的具体语义差异，确认可替换后再动模板；
- **`paf coverage` 无 cg:Z gap**：ovlpr `covered` 用 PAF start/end（不依赖
  CIGAR），pgr `paf coverage` 要求 `cg:Z`（无标签记录不贡献）——补丁方向：
  无 cg:Z 时退回用 start/end 算覆盖；
- **anchr 模板替换整体**：见 todo §3 挂账条目（trim.era.sh 等 bbtools
  调用换 pgr 命令链，用户自处理）。
