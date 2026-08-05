# pgr pgi + LASTZ 混合比对（FastGA-gapfill 思路）

> 设计笔记（未实现）。背景：FASTGA 论文（*FastGA: fast genome alignment*，
> Bioinformatics Advances 5(1):vbaf238，DOI 10.1093/bioadv/vbaf238）提出混合方案
> FastGA-gapfill——以 FastGA 比对为锚点，对每对连续同向锚点之间的区间跑 LASTZ
> 填 gap，最后合并两者结果。pgr 的本地对应物：`pgr align pgi`（原生 FastGA 风格）
> 快速找片段，`pgr align lastz` 精修，两套 PSL 合并喂 `pgr pl chainnet`。
> 日期：2026-08-06。状态：方案已讨论，未实现（不着急）。
> 关联：[[pgi-align.md]]、[[sd.md]]、[[references/fastga.md]]、[[paf-pangenome.md]]。
> 命令命名（2026-08-06 与用户讨论）：不做 `pgr pl align --hybrid`，直接做
> `pgr align hybrid`——`pgr pl` 是"暂时没想好该放到哪边"的命令的临时存放处，
> hybrid 放 `pgr align` 下方便用户发现。

## 1. 背景与目标

- FASTGA 论文里 FastGA 找的是"几乎横跨整段序列、最大 gap 约 40 bp"的局部比对，
  不自己做跨大 gap 的 chaining；论文用 LASTZ 填锚点之间的 gap，把 FastGA 的速度
  和 LASTZ 的灵敏度（论文实测最高）结合起来。
- pgr 侧目标：`pgr align pgi` 做全基因组快速粗比对 → 对 pgi 没覆盖到（或锚点之间）
  的区间跑 LASTZ 补比对 → 两套 PSL 合并后进 chain/net 流程，得到更完整的比对结果。
- 意义：pgi 边界一般比真实同源区短 1–11 bp（链的种子锚定边界，见 `docs/sd.md`），
  且对接近 SD 身份下限（~90–93%）的拷贝可能漏块；LASTZ 补的正是这类缺口。

## 2. 现有资产盘点（全部已实现）

- `pgr align pgi`：FastGA 风格原生比对（syncmer 种子 → tube 链 → wave 扩展），
  输出 PSL。核心在 `src/libs/pgi/`（build/align/wave/mmap）+ `src/libs/ds/radix_sort.rs`
  （FastGA MSDsort 移植）。
- `pgr align lastz`：LASTZ 封装（`src/libs/lastz.rs`），输出 **LAV**，7 套预设
  set01..set07（set01 Human vs Chimp 最快，set07 Human vs Opossum 最远缘/灵敏），
  要求每个 FASTA 单序列，多 contig 需先 `pgr fa split name`。
- `pgr lav to-psl`：LAV → PSL 转换。
- `pgr pl chainnet`：native psl-chain-net-axt-maf 全链路，PSL 参数接受**文件或目录**
  （多个 PSL 直接合并链化）；`scripts/verify-pangenome.sh` 已验证
  FastGA PSL → chainnet 路径，pgi 的 PSL 与之格式一致。
- `pgr fa range` / `pgr fa split name`：按坐标提取区间、按记录切分序列。

## 3. 方案设计

### 3.1 流程

1. 粗比对：`pgr align pgi target.fa query.fa -o pgi.psl`
2. 确定补比对区间：
   - 首选（论文语义）：pgi 相邻同向锚点之间的 target/query 双侧 gap；
     可先 `pgr psl chain` 拿锚点骨架再取相邻 block 之间的空档。
   - 备选：全基因组无覆盖区；超长无锚点区间应跳过（多为真正的新序列）。
3. 提取区间序列并切单序列：`pgr fa range` / `pgr fa split name`
4. 精修：`pgr align lastz --preset <预设>`（预设由用户选择，见 §3.2/§3.6）→ LAV
5. 转换：`pgr lav to-psl` → `lastz.psl`
6. 合并去重：丢弃与 pgi 块重叠 >50% 的 lastz 块（或保长去短），避免 chain 分支碎片化
7. ChainNet：
   `pgr pl chainnet [--syn] target.fa query.fa psl_all/ -o out`

### 3.2 关键决策点

- **补比对区间**：锚点间 gap（命中率高、链得进去）vs 全未覆盖区（简单但白跑多）。
- **去重规则**：重叠阈值与保长/保 pgi 的选择，需小规模数据实测。
- **`--syn`**：syntenic 共线性比对加；重复/SD 分析必须不加（`pgr sd align` 明确
  规定 chain/net 精修非 `--syn`，否则重排同源丢失）。
- **预设**：做成用户选项（复用 `pgr align lastz --preset set01..set07`）。
  泛基因组主场景差异小，默认贴近的 set01/set02；远缘比较由用户自选
  set06/set07——**不默认远缘预设**（2026-08-06 与用户讨论后调整）。
- **方向一致性**：pgi 与 lastz 的 target/query 顺序必须一致，PSL 的 tName/qName
  前缀统一，否则 chain 链向混乱。

### 3.3 边界处理策略（定稿，2026-08-06 与用户讨论后）

pgi 块边界一般比真实同源区短 1–11 bp。处理方式**不是缩短 pgi 的 PSL 记录**
（那会丢失已找到的片段，且需 t/q 同步变换的新函数），而是：

1. pgi 的 PSL **原样保留**，一个区间都不动；
2. 仅当计算"补集区间"（交给 LASTZ 的范围）时，对 pgi 的 target 侧区间做
   `trim(n)`（n ≈ 25–50 bp，大于 pgi 边界误差 1–11 bp 即可）——补集因此向外
   大出一个缓冲带，真实边界落进 LASTZ 的搜索范围；
3. LASTZ 跑补集，覆盖真实边界，与 pgi 完整块边界产生**少量有意重叠**；
4. 合并时 pgi 完整块 + lastz 块，重叠 >50% 保长去重兜底（chain 分支防护）。

实现上只需要一维区间运算：`pgr psl to-range` → `pgr runlist span --op trim`
→ `pgr runlist span --op holes`，全部是现有命令，**不需要新增 PSL 记录变换
函数**。这块逻辑与 `rept`/`ir` 的 `run_repeat_runlist_pipeline`
（`libs/pl/repeat.rs` 的 Fill → Excise → Fill）同源，复用的正是
`IntSpan::trim/holes/fill/excise` 与 `pgr runlist span`。

### 3.4 适用场景（定稿）：仅共线性搜索

Hybrid（pgi 锚点 + LASTZ 补 gap）模式**只适合共线性（syntenic）搜索**，
即 `pgr pl chainnet --syn` 一路。原因：

- 非共线性（SD/重复）场景下 net 不做 syntenic 过滤，重叠冗余块会一起保留，
  输出碎片化、覆盖重复计数——去重只能缓解，不是正解；
- SD 场景的正确解法是给 `sd search` 补 BISER 式 `MAX_EXTEND` 边界扩展
  （见 [[sd.md]] §4.8），而不是引入混合比对；
- 本方案不计划覆盖非共线性用途。

### 3.5 论文 FastGA-gapfill 参数对照（2026-08-06 通读论文后补充）

论文 §5.2 的 FastGA-gapfill 与本方案的对应关系（详见
[[references/fastga.md]] §12.3）：

- **补 gap 的对象**：每对"顺序一致、方向一致、不重叠"、间隔 ≤1 Mb（默认）的
  锚点 → 双侧 bounding box。与我们"锚点间 gap / 未覆盖区补集"一致，且同样
  隐含共线性前提（§3.4）。
- **重叠缓冲**：论文默认 box 与锚点重叠 **1 kb** 以利 LastZ 播种；ALNfill
  `alngap -e` 默认也是 1 kb。我们的 §3.3 目前用 trim 25–50 bp——比论文小两个
  数量级。**建议实测时对比 50 bp vs 500 bp vs 1 kb**：重叠越大 LastZ 播种越稳、
  但合并时冗余/重叠越多（ALNfill 已把 `-e` 造成的重叠列为已知问题）。
- **box 去嵌套**：论文只保留最小 bounding box（无包含关系）——对应我们的
  `runlist holes` 天然产出非重叠区间，无需额外处理。
- **合并方式**：论文直接把两套输出 cat 合并（PAF），未做去重；我们计划保留
  >50% 重叠保长去重兜底（§3.2），仍值得保留。

> 论文的 FastGA-gapfill 灵敏度接近 LastZ、速度比 LastZ 快 19.3×–137.5×——
> 这是本方案可行性的直接证据，也是后续验证的对照目标。

### 3.6 ALNfill 实现对照（2026-08-06 源码通读后补充）

ALNfill（`alnfill-main/`，Chenxi Zhou）是论文 FastGA-gapfill 的工程化实现，
源码分析见 [[references/alnfill.md]]。对本方案有直接影响的实现细节：

- **gap 过滤**：`alngap` 只补双侧 gap 都在 [100, 1M] 的区间（`-l`/`-m`），
  且用哨兵覆盖染色体首尾端。我们的 holes 方案应加同样过滤：
  小于 100 bp 的洞让 pgi 自己处理，大于 1 Mb 的洞跳过（与 §4"超长无锚点区间
  应跳过"一致）。
- **去冗余**：`alngap` 默认对 PAF 做"双侧被覆盖 ≤50% 的贪心过滤"
  （reciprocal best，`-a` 关闭）。pgi PSL 在重复区可能有多映射块，
  算 holes 前值得先按链上锚点去冗余。
- **方向**：`alngap` 不读 PAF strand 列，混合方向的锚点对也会成 box；
  论文描述是"一致顺序方向"。我们定稿为仅共线性（§3.4），实现时可加
  方向过滤或直接交给 `chainnet --syn` 收尾。
- **LastZ 选项**：ALNfill 只用 `--format=PAF:wfmash --ambiguous=iupac`
  （lastz 默认打分）；pgr 复用 `pgr align lastz --preset`，**预设由用户选择**——
  泛基因组差异小，默认贴近的 set01/set02，远缘比较才选 set06/set07，
  不默认远缘（2026-08-06 与用户讨论后调整）。
- **坐标回移**：ALNfill 提取区间、跑 lastz、把区间坐标回移成全长坐标、输出完整
  PAF 一步完成；pgr 若做成 `pgr align hybrid`，需要 `fa range` 提取时记录
  offset，lastz 输出回移后再并 PSL。
- **内存**：ALNfill 把两个基因组整库读进内存（sdict strdup），超大基因组不现实；
  pgr 的 2bit/loc 区间提取无此问题。

## 4. 已知的坑

- LASTZ 输入单序列限制 → 提取/切分是主要工程量。
- `pgr align lastz` 只出 LAV，需 `pgr lav to-psl`。
- pgi 边界短 1–11 bp → lastz 补块与 pgi 块边界重叠 → 合并前去重。
  决策见 §3.3：pgi PSL 原样保留，仅 trim 补集计算区间，不缩 PSL 记录本身。
- 全未覆盖区含物种特异插入/着丝粒等无同源序列，LASTZ 白跑；应优先锚点间 gap。
- `pgr align lastz` 默认 query-depth 50 是"先到先得"式截断，补 gap 场景若覆盖深
  可能丢块，必要时调大 `--query-depth`。

## 5. 验证方案（后续执行）

- 用 MG1655 vs Sakai（已有测试数据）跑 pgi-only / hybrid / lastz-only 三路，
  对比 chainnet 输出的覆盖率与链完整性；hybrid 目标：覆盖接近 lastz-only，
  耗时接近 pgi-only（论文 gapfill 的结论形态）。
- 检查 hybrid 输出是否碎片化、有无重复覆盖块。
- **灵敏度口径（借鉴论文 §5.1）**：模拟含已知目标区的序列对（10 kb 块 = 相似区 +
  随机序列、目标区长度 100–5000 bp × 分歧度 1–40%、每组合 100 重复），
  统计"目标区被比对覆盖 ≥95%（两侧都算）"的恢复比例；同时统计目标区外的
  假阳性碱基比例。pgr 侧用 `pgr psl to-range` / `pgr runlist` 做覆盖统计即可。
- **真实数据口径（借鉴论文 §12.4，[[references/fastga.md]]）**：覆盖统计只按
  比对 start/end 计（比对内 gap 也算覆盖）；LastZ 侧可考虑喂 soft-mask 序列、
  按染色体对跑（与论文实验一致，作为对照时口径要对齐）。

## 6. 待办

- [ ] 手动脚本跑通小规模验证（临时，不入库）。
- [ ] 若效果达标，做成 `pgr align hybrid` 子命令（见头部命名说明），
      复杂逻辑放 `src/libs/`（如 `libs/align/`），`cmd_pgr` 保持薄壳。
- [ ] 按 AGENTS.md 要求补 `cargo fmt` / `cargo clippy` 与测试（如进入实现阶段）。
