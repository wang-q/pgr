# pbit PAF 驱动模式设计提案

> **状态**：设计提案，尚未实现。基于 [pbit.md](pbit.md) 的现有格式，探讨用 PAF 比对结果驱动
> LZ-diff 参考选择、以及直接存储 CIGAR 替代 LZ-diff 的可行性与方案。

## 1. 动机

### 当前局限

pbit 当前的参考选择是**按段位置索引匹配**（非 minimizer；[compressor.rs:280-336](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/compressor.rs#L280-L336)）：

- 按 contig 名匹配参考段（`contig_ref_groups[contig_name]`）
- 按段位置索引匹配（`ref_group_ids[seg_idx]`，clamped 到最后一段）
- 首段 k-mer 采样**仅用于方向检测**（`detect_rev_comp`，非参考选择）+ 逐段 delta 大小回退

**无法处理**：重排、转位、跨 contig 比对、大段 indel 导致的位置偏移。对标准群体基因组
（同源染色体对齐）足够，但对含结构变异的多样性基因组压缩率差。

### PAF 驱动的优势

用户用 minimap2/wfmash 等工具将样本比对到参考，生成 PAF（含精确 `=/X` CIGAR）。pbit 可利用
该比对结果：

1. **精确参考选择**：PAF 直接给出每段样本对应哪个参考段（含跨 contig/重排）
2. **精确方向**：PAF 的 strand 字段比 k-mer 采样检测更可靠
3. **精确差异**：CIGAR 的 `=/X/I/D` 操作已完整描述样本与参考的差异，可直接存储
4. **可复用比对信息**：解压后可还原比对关系，支持变异分析，无需重新跑比对工具

## 2. pgr 已有的 PAF + CIGAR 基础设施

pgr 已有完整的 PAF 处理栈，pbit 可直接复用：

| 组件 | 位置 | 说明 |
|------|------|------|
| `CigarOp` | [libs/paf/cigar.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/cigar.rs) | bit-packed u32（3位op + 29位length），已有 `from_raw`/`op()`/`len()`/`target_delta`/`query_delta` |
| `CigarStore` | [libs/paf/index/mod.rs:22](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/index/mod.rs#L22) | `Owned(Vec<CigarOp>)` / `Lazy(u64)` / `LazyReversed(u64)` |
| `PafIndex` | [libs/paf/index/](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/index) | coitrees 区间树索引，按坐标查询比对 |
| `build_pairwise_block` | [libs/paf/msa_build.rs:191](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/msa_build.rs#L191) | 从 CIGAR + FastaStore 重建比对序列（含正负链处理） |
| `extract_cigar` | libs/paf/cigar.rs | 从 PAF tags 解析 `cg:Z:` CIGAR 字符串为 `Vec<CigarOp>` |
| `reverse_cigar` | libs/paf/cigar.rs | 反转 CIGAR（用于负链） |
| `FastaStore` | libs/paf/fasta.rs | 序列存储，按名+区间提取 |

**关键洞察**：`CigarOp` 的 bit-packed u32 格式天然适合 pbit 的二进制存储。pgr 已有从 CIGAR +
参考序列重建样本序列的完整逻辑（`build_pairwise_block` → `build_maf_block`），pbit 解压时可复用。

## 3. 方案选项

### 方案 A：PAF 仅指导参考选择（LZ-diff 不变）

PAF 用于选择参考段和方向，delta 仍用 LZ-diff 编码。PAF 不存储在文件中。

```
压缩：PAF → 确定段级对应 → LZ-diff 编码差异 → flate2 压缩
解压：LZ-diff 解码 → 重建样本（与当前一致）
```

- **优点**：改动最小，文件格式不变，解压端无需改动
- **缺点**：不保留比对信息；LZ-diff 对含 indel 的段编码效率不如 CIGAR（indel 导致后续位置错位，
  LZ-diff 的 match 失效，退化为大量 literal）

### 方案 B：PAF 指导 + 存储比对元数据

PAF 指导参考选择，在 Collection 中存储比对摘要（ref_group_id, strand, query/target 坐标），
delta 仍用 LZ-diff 编码。

```
压缩：PAF → 确定段级对应 → LZ-diff 编码差异 → flate2 压缩；存储比对摘要
解压：LZ-diff 解码 → 重建样本；比对摘要可供查询
```

- **优点**：保留比对关系，解压后可查询；格式变化小
- **缺点**：仍用 LZ-diff 编码（indel 效率问题）；比对摘要与 LZ-diff delta 有冗余

### 方案 C：CIGAR 替代 LZ-diff 作为 delta 编码（推荐）

用 PAF 的 CIGAR 直接作为 delta 编码，替代 LZ-diff。CIGAR 的 `=/X/I/D` 操作完整描述差异，
解压时按 CIGAR 从参考序列重建样本序列。

```
压缩：解析 PAF → 提取每段 CIGAR + X/I 碱基 → bit-packed → flate2 压缩
解压：flate2 解压 → Vec<CigarOp> + X/I 碱基流 → 按 CIGAR 从参考段重建样本段
```

> `packed_data` 的完整二进制格式（CIGAR ops + X/I 碱基流）见 §8 "X/I 操作的碱基存储问题"。

- **优点**：
  - 复用 `CigarOp` bit-packed 编码（4 字节/op，高效）
  - 不需要重新计算 LZ-diff（省压缩时间）
  - 天然处理 indel（I/D 操作，不受位置错位影响）
  - 解压逻辑与 `build_pairwise_block` 一致
  - 保留比对信息（CIGAR 本身就是比对结果）
- **缺点**：
  - 依赖外部比对工具（用户需先跑 minimap2/wfmash）
  - 需处理未覆盖区域（PAF 没覆盖的样本序列）
  - 文件格式需扩展（delta 编码类型标志）

### 方案 D：混合编码（LZ-diff + CIGAR 可选）

delta 层支持两种编码，每段独立选择：
- 有 PAF 覆盖的段 → CIGAR 编码
- 无 PAF 覆盖的段 → LZ-diff 编码（回退到当前 minimizer 逻辑）

```
DeltaEntry {
    encoding: u8,  // 0 = LZ-diff, 1 = CIGAR
    packed_data: Vec<u8>,  // flate2(编码数据)
}
```

- **优点**：兼顾两种场景；未覆盖区域不丢数据
- **缺点**：解压端需支持两种解码路径；复杂度增加

## 4. 推荐方案：C + D 混合

**推荐方案 C 为主，D 为回退**：

- 有 PAF 输入时：PAF 覆盖的段用 CIGAR 编码，未覆盖段回退 LZ-diff
- 无 PAF 输入时：全部用 LZ-diff（当前行为，完全向后兼容）

理由：
1. CIGAR 对含 indel/重排的段编码效率优于 LZ-diff
2. `CigarOp` bit-packed 格式与 pbit 二进制风格一致
3. pgr 已有完整的 CIGAR 重建逻辑可复用
4. 混合模式保证未覆盖区域不丢数据，且与现有格式向后兼容

## 5. 格式扩展

### 5.1 DeltaEntry / DeltaMeta 扩展

pbit 现有**两个**结构体共享同一在盘头部：
- `DeltaMeta`（9 字节，仅头部）—— `Decompressor::new` 扫描所有 delta 头部构建
  `delta_meta: Vec<Vec<DeltaMeta>>`（[decompressor.rs:42](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/decompressor.rs#L42)），**不读 packed_data**
- `DeltaEntry`（完整，含 packed_data）—— 按需读取，通过 `meta()` 派生 `DeltaMeta`

`encoding` 字段**必须加入 `DeltaMeta`**：`get_contig` 决定走 LZ-diff 还是 CIGAR 解码路径时只看
`delta_meta`（不读 packed_data）。若 `encoding` 不在 `DeltaMeta`，则必须先读完整 delta 才能判断编码
——破坏随机访问。

```rust
pub struct DeltaMeta {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_size: u32,
    pub encoding: DeltaEncoding,  // 新增（在盘第 10 字节）
}

pub struct DeltaEntry {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_data: Vec<u8>,
    pub encoding: DeltaEncoding,  // 新增；meta() 自动带上
}

pub enum DeltaEncoding {
    LzDiff = 0,
    Cigar = 1,
}
```

**在盘格式**（Delta Data 区，每条 delta）：

```
offset  size  field
0       1     is_rev_comp
1       4     raw_length
5       4     packed_size
9       1     encoding        ← 新增（0 = LZ-diff, 1 = CIGAR）
10      packed_size packed_data
```

> encoding 字段是新增的。现有文件的 delta 头部是 9 字节（无 encoding），新文件为 10 字节。
> `DeltaMeta::read_header` / `write_header` 改为读写 10 字节。
> 本项目尚未正式发布，无需考虑向后兼容，直接升级 Header 版本号（1000 → 1001）。
>
> **raw_length 在 CIGAR 模式下的语义**：LZ-diff 的 `raw_length` 是样本段长度；
> CIGAR 模式下 `raw_length = Σ op.query_delta()`（query 轴长度，含 `X`/`I`，不含 `D`），
> 即样本段长度。`get_contig` 用 `raw_length` 做切片计算，CIGAR 段必须提供正确的 query 轴长度。

### 5.2 SegmentDesc 扩展（CIGAR 模式，固定大小存储）

CIGAR 模式下，每段需要记录比对坐标（样本段在参考段中的偏移）。采用**固定大小存储**：
SegmentDesc 新增 `ref_start` / `ref_end` 两个 u32 字段，LZ-diff 段填 0。

```rust
pub struct SegmentDesc {
    pub ref_group_id: u32,
    pub delta_id: u32,
    /// Segment-relative start offset within the reference 2bit record
    /// (= `target_start - seg_idx * segment_size`). 0 for LZ-diff segments.
    pub ref_start: u32,
    /// Segment-relative end offset (exclusive) within the reference 2bit
    /// record (= `target_end - seg_idx * segment_size`). 0 for LZ-diff segments.
    pub ref_end: u32,
    // is_rev_comp / raw_length / encoding 由 DeltaEntry 提供
}
```

**在盘格式**（Sample Index 区，每段，固定 16 字节）：

```
offset  size  field
0       4     ref_group_id
4       4     delta_id
8       4     ref_start         ← 新增（CIGAR 段：参考段内相对偏移；LZ-diff 段：0）
12      4     ref_end           ← 新增（CIGAR 段：参考段内相对结束；LZ-diff 段：0）
```

> **坐标语义（相对坐标）**：`ref_start`/`ref_end` 存储的是**相对于参考段（2bit 记录）起始的
> 偏移**，而非 PAF 的 contig 绝对坐标。计算方式：`seg_idx = target_start / segment_size`，
> `ref_start = target_start - seg_idx * segment_size`，
> `ref_end = target_end - seg_idx * segment_size`。解压时 `read_2bit_record` 读出的参考段
> 从 `seg_idx * segment_size` 开始，直接用 `ref_start..ref_end` 切片即可定位 CIGAR 对应的
> 参考区间，无需反查 `contig_ref_groups`。
>
> **为什么需要 ref_start/ref_end**：CIGAR 描述的是"样本段 vs 参考段某区间"的差异。
> 同一参考段可能被多个样本段引用（不同区间），需要 ref_start/ref_end 定位。LZ-diff 段
> 参考段是整条 2bit 记录，ref_start=0、ref_end=记录长度，但为保持固定大小统一填 0
> （解压时 LZ-diff 路径不读这两个字段）。
>
> **为什么不用变长存储**：变长方案（按 encoding 决定是否后跟 ref_start/ref_end）会破坏
> `SegmentDesc` 的 `Copy` trait，并让 serialize/deserialize 复杂化。固定大小方案虽然 LZ-diff
> 段多占 8 字节，但 Collection 整体经 flate2 压缩，连续零值压缩率极高，实际开销可忽略。
>
> **encoding 位置**：encoding 只放 `DeltaMeta`（§5.1），不放 SegmentDesc。一个 delta 的
> encoding 是其 packed_data 的固有属性，不会因引用它的 segment 不同而变化。`Decompressor::new`
> 扫描 delta 头时即可获知编码类型，`get_contig` 据
> `delta_meta[ref_group_id][delta_id].encoding` 分支解压路径，无需读 SegmentDesc。

### 5.3 CIGAR 的 bit-packed 存储

`Vec<CigarOp>` 直接序列化为 `Vec<u32>`（每个 CigarOp 是一个 u32），配合 X/I 碱基流，
然后 flate2 压缩。完整 `packed_data` 格式见 §8。

```
packed_data = flate2( u32 op_count + [CigarOp; op_count] + u32 base_count + [u8; base_count] )
```

解压时：`flate2_decompress(packed_data) → (Vec<CigarOp>, Vec<u8> X/I bases)`

> **空间效率**：一个 4096 bp 的纯 match 段 → `[=4096]` → 1 个 CigarOp = 4 字节 → flate2 后更小。
> 含 40 个 SNP 的段 → `=99 X1` × 40 + `=96` ≈ 82 个 CigarOp = 328 字节 → flate2 后约 100-200 字节。
>
> **压缩率 caveat**：上述估算需基准验证。CigarOp 是 `(op << 29) | len` 的 bit-packed u32，op 占
> 高位 3 bit、len 占低 29 bit。相邻 CigarOp 的 op 字段可能有局部重复（如 `=99 X1 =99 X1...`），
> 但 len 字段差异大。DEFLATE 对这种半结构化 u32 数组的压缩率通常不如对文本/LZ-diff 字节流。
> **SoA 优化候选**（Phase 8f）：将 ops 与 lengths 分离存储（两个独立数组），可提升 DEFLATE 对
> op 字段（高度重复）的压缩效率。作为可选优化，默认采用 interleaved（AoS）布局。

## 6. CLI 设计

### 6.1 新增 `--paf` 参数

```
pgr pbit create -r ref.fa -i sample.fa --paf sample.paf -o out.pbit
pgr pbit append in.pbit -i sample.fa --paf sample.paf
```

- `--paf` 可选，指定样本比对到参考的 PAF 文件
- 一个 `-i` 对应一个 `--paf`（或支持一个 PAF 含多个样本的比对）
- 省略 `--paf` 时使用当前 LZ-diff 模式（按段位置索引匹配参考，见 §1）

### 6.2 参数语义

```
pgr pbit create -r ref.fa \
    -i sample1.fa --paf sample1.paf \    # sample1 用 CIGAR 模式
    -i sample2.fa                        # sample2 用 LZ-diff 模式（无 PAF）
```

> 同一归档可混用两种模式（每样本独立选择）。delta 层的 encoding 字段区分。

### 6.3 `--name` TSV 模式与 `--paf` 的关系

现有 `--name` TSV 接受两列 `sample_name<TAB>fasta_path`（[create.rs:84-88](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/pbit/create.rs#L84-L88) /
[append.rs:48-52](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/pbit/append.rs#L48-L52)）。为支持 PAF，
扩展为**可选三列**：

```
# samples.tsv
sample1    sample1.fa    sample1.paf    # CIGAR 模式
sample2    sample2.fa                   # LZ-diff 模式（第三列缺失）
sample3    sample3.fa                   # LZ-diff 模式
```

**互斥规则**：
- `--paf` 仅与 `-i` 配套（`-i` 与 `--paf` 一一对应，按出现顺序配对）
- `--name` 与 `--paf` **互斥**（`--name` 模式下 PAF 由 TSV 第三列指定，不允许同时用 `--paf`）
- `--name` TSV 第三列缺失时该样本走 LZ-diff（向后兼容现有两列 TSV）

> **推荐**：多样本 + PAF 场景优先使用 `--name` TSV 三列模式。`-i`/`--paf` 按顺序配对
> 在样本数较多时容易出错（如 `-i s1.fa -i s2.fa --paf s1.paf --paf s2.paf` 与
> `-i s1.fa --paf s1.paf -i s2.fa --paf s2.paf` 顺序不同但意图相同，易混淆）。TSV 模式
> 将样本名、FASTA、PAF 显式绑定在同一行，无歧义。CLI 实现时对 `-i`/`--paf` 配对应做
> 计数校验（数量必须相等，否则报错）。

## 7. 压缩流程（CIGAR 模式）

```
append_sample_with_paf(sample_name, fasta_path, paf_path):
  1. 解析 PAF → 构建 query-side 区间树（按 query 坐标索引，见下方"PAF query-side
     索引"说明）
  2. 读取样本 FASTA → 按 segment_size 分段
  3. 对每个样本段 [seg_start, seg_end):
     a. 查 query-side 区间树：哪些 PAF alignment 覆盖此样本段区间？
     b. 选最佳 alignment（按 identity / 覆盖度）
     c. 提取该段对应的 CIGAR 片段（从 PAF CIGAR 中按坐标截取）
     d. 若 CIGAR 含 `M` 操作：对比样本碱基与参考碱基，拆分为 `=/X`
        （fallback：`pgr maf to-paf` 已在源头区分 `=/X`，此步骤仅处理
        minimap2/fastga 等直接输出的 `M` CIGAR）。参考碱基从
        `self.segments[ref_group_id].reference_dna()` 按 CIGAR target 坐标切片获取
        （target 坐标减去 `seg_idx * segment_size` 转为参考段内偏移）。
     d'. 提取 X/I 碱基：按 CIGAR 正向遍历顺序，收集 `X`/`I` 操作对应的样本碱基。
        **`-` 链记录**：CIGAR 描述的是 RC(query) vs forward(target) 比对，因此 X/I 碱基
        必须从 **RC(sample)** 提取（与 CIGAR 描述的比对方向一致），而非原始正向样本序列。
        解压时先正向应用 CIGAR（用存储的 RC 样本碱基）→ 得到 RC(sample) → 再 RC 得到 sample
        （见 §8 负链语义）
     e. 确定 ref_group_id（PAF target → 参考段映射）+ strand + ref_start/ref_end
     f. bit-pack CIGAR → flate2 压缩 → DeltaEntry(encoding=CIGAR)
     g. 记录 SegmentDesc(ref_group_id, delta_id, ref_start, ref_end)
     h. Delta 去重按 packed_data 字节比较（与 LZ-diff 一致），相同 CIGAR+XI 的段共享 delta_id
  4. 未覆盖段：回退 LZ-diff（对到同名参考段，当前逻辑）
```

### 关键问题：PAF 坐标→段映射

PAF alignment 的 query_start/end 是连续坐标，pbit 按 segment_size 切段。一条 PAF alignment
可能横跨多段。处理方式：

1. **保持固定段切分**：样本仍按 segment_size 切分
2. **CIGAR 截取**：对每段，从覆盖它的 PAF alignment 的 CIGAR 中截取对应区间
   - 用 `CigarOp::target_delta` / `query_delta` 做坐标投影
   - 一段可能被多条 alignment 覆盖 → 选最佳（identity 最高 + 覆盖度最大）
   - 一段可能部分被覆盖或跨多条 alignment 衔接点 → **整段回退 LZ-diff**（见 §11 决策 3，
     不拆段、不混用 CIGAR 与 LZ-diff）
3. **CIGAR target 投影跨参考段边界检查**：样本段按 query 坐标切分，但其 CIGAR 投影到
   target 轴时，因 indel 偏移，target 区间 `[target_start, target_end)` 可能跨越参考段边界
   （即 `target_start / segment_size != (target_end - 1) / segment_size`）。此时单个
   `ref_group_id` 无法覆盖完整 target 区间，`ref_end` 会超出参考段长度。处理方式：
   **整段回退 LZ-diff**（与决策 3 一致，不拆段、不跨段引用多 ref_group）。
   此情况在 indel 密集或段尾接近边界时可能发生，但 segment_size（默认 4096）远大于
   典型 indel，实际命中率影响小。

### 关键问题：PAF target → ref_group_id 映射

PAF 的 `target_name` + `target_start`/`target_end` 是 contig 内坐标，需映射到 pbit 的
`ref_group_id`。ref_groups 按 contig 分组、按 segment_size 切段顺序排列
（[compressor.rs:233-263](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/compressor.rs#L233-L263)），
`contig_ref_groups: IndexMap<String, Vec<u32>>` 提供 contig 名 → ref_group_id 列表。

映射规则：
```
ref_group_id = contig_ref_groups[target_name][target_start / segment_size]
```

参考段 i 覆盖 contig 内区间 `[i * segment_size, (i+1) * segment_size)`（最后一段可能短于
segment_size，需 clamped）。`ref_start`/`ref_end` 存储**相对参考段起始的偏移**（见 §5.2 坐标
语义），即 `ref_start = target_start - seg_idx * segment_size`、
`ref_end = target_end - seg_idx * segment_size`，其中 `seg_idx = target_start / segment_size`。
解压时 `read_2bit_record` 读出的参考段从 `seg_idx * segment_size` 开始，直接用
`ref_start..ref_end` 切片定位 CIGAR 对应的参考区间。

### 关键问题：未覆盖区域

PAF 没覆盖的样本序列区域：
- **回退 LZ-diff**：对到同名参考段（当前逻辑），encoding=LzDiff
- **存 raw**：直接存储原始序列（无参考），encoding=Raw（需第三种 encoding）
- **推荐**：回退 LZ-diff（复用现有逻辑，格式变化最小）

### 关键问题：PAF query-side 索引

pbit 压缩时需按**样本（query）段坐标**查询覆盖该段的 PAF alignment。现有
`PafIndex`（[libs/paf/index/mod.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/index/mod.rs)）
的 `reverse_trees` 虽按 query 坐标索引，但**不适合 pbit 直接使用**：

1. **仅覆盖 `+` 链记录**：`insert_record`（[mod.rs:179](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/index/mod.rs#L179)）
   的 mirror entry 只在 `rec.strand == '+'` 时插入。`-` 链记录在 `reverse_trees` 中无条目。
2. **元数据被交换**：mirror entry 的 query_id = 原 target_id，target_start/end = 原 query
   坐标，strand 强制为 `'+'`，CIGAR 被 `reverse_cigar`（I/D 交换）。这是为 BFS 设计的
   "角色互换"视图，与 pbit 需要的原始 PafMetadata 不一致。
3. **无公开查询方法**：`query()`（[query.rs:44](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/index/query.rs#L44)）
   只查 `self.trees`（target 侧）。`reverse_trees` 是 `pub(crate)`，无公开 query 接口。

**方案**：pbit 自建 query-side 区间树。读取 PAF 后，构建独立的
`BasicCOITree<PafMetadata, u32>`（key = query_id，interval = query 坐标），存储**原始未交换**
的 PafMetadata（含原始 strand、原始 target 坐标、原始 CIGAR）。不依赖 `reverse_trees`，
不改现有 PafIndex 的行为。

> **query_name → query_id 映射**：pbit 压缩端按 contig **名**（非 id）处理样本段，PAF 的
> `query_name` 是字符串。构建区间树时先建立 `IndexMap<String, u32>` 映射 query_name →
> query_id（按 PAF 出现顺序分配），区间树以 `query_id` 为 key。压缩端遍历样本段时按
> contig_name 查此映射获取 query_id，再查对应的区间树。

> 实现位置：`libs/pbit/` 内新增 query-side 索引构建逻辑（复用 `libs/paf/` 的 PAF 解析 +
> `coitrees` 区间树，但不依赖 `PafIndex` 结构体）。PAF 记录若无 `cg:Z:` CIGAR 标签
> （`extract_cigar` 返回空 Vec），该记录不插入索引（视为未覆盖，相关段回退 LZ-diff）。

### append 兼容性

`pgr pbit append` 复用 `Compressor::open_for_append`（[compressor.rs:198-267](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/compressor.rs#L198-L267)），
CIGAR 模式仅影响新追加的样本，已有数据不变。`open_for_append` 重建的 `Segment` 对象
（调用 `prepare` / `prepare_index`）供 LZ-diff 回退段使用，CIGAR 段不依赖它。新追加样本
若带 `--paf` 则走 CIGAR 路径，不带则走 LZ-diff 路径，与已有样本的编码类型互不影响。

## 8. 解压流程（CIGAR 模式）

```
get_contig / get_sample (CIGAR 段):
  1. 读取 DeltaEntry → 检查 encoding（从 DeltaMeta 获取，无需读 packed_data）
  2. 若 encoding == CIGAR:
     a. flate2 解压 packed_data → Vec<u32> → Vec<CigarOp> + X/I 碱基流
        （X/I 碱基按 CIGAR 正向遍历顺序连续存储，解压时按相同顺序消费）
     b. 读取参考段（ref_group_id → seek → read_2bit_record）
     c. 按 SegmentDesc.ref_start/ref_end 截取参考段区间（正向坐标）
     d. 按 CIGAR 重建样本段（正向应用 CIGAR，is_rev_comp 见下方语义说明）:
        - '=' : 从参考取 len 个碱基
        - 'X' : 跳过参考 len 个碱基，从 X/I 碱基流取 len 个样本碱基
        - 'I' : 从 X/I 碱基流取 len 个样本碱基（参考不前进）
        - 'D' : 跳过参考 len 个碱基
     e. 若 is_rev_comp: 对重建结果做反向互补
     f. 拼接所有段
  3. 若 encoding == LZ-diff: 当前逻辑（LZ-diff 解码）
```

> **`decode_delta` 签名扩展**：当前 `decode_delta(ref_group_id, delta_id)`
> （[decompressor.rs:260](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/decompressor.rs#L260)）
> 不接收 ref_start/ref_end。CIGAR 模式下步骤 2c 需要 SegmentDesc 的 ref_start/ref_end 来切片
> 参考段。实现时将 `decode_delta` 签名扩展为接收 `&SegmentDesc`（或额外传入 ref_start/ref_end），
> 由 `get_contig`/`get_sample` 在调用时从 SegmentDesc 提取。LZ-diff 路径忽略这两个字段。

### 负链 `is_rev_comp` 语义（CIGAR 模式）

`is_rev_comp` 来源于 PAF strand 字段，与 [build_pairwise_block](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/msa_build.rs#L191-L220)
的处理一致：

- **压缩时**：`is_rev_comp = (PAF strand == '-')`。CIGAR 直接取自 PAF（不重写为正向）。
  ref_start/ref_end 始终是**正向参考坐标**（PAF target 坐标，target 在 PAF 中总是正向）。
- **解压时**：先按 CIGAR 在正向参考区间上重建（得到 rev-comp 后的样本段），再对结果做
  `nt::rev_comp` 得到原始样本段。等价于 `build_pairwise_block` 中对 '-' 链先 rev-comp query
  再正向应用 CIGAR 的逆操作。

> **不变量**：ref_start/ref_end 永远是正向参考坐标；is_rev_comp 仅影响样本段方向。
>
> 此操作是 `build_pairwise_block` 的逆——后者对 '-' 链先 RC query 再正向应用 CIGAR，
> 前者先正向应用 CIGAR 再 RC 结果。

### X/I 操作的碱基存储问题

CIGAR 的 `X`（mismatch）和 `I`（insertion）操作需要存储样本特有的碱基（参考中没有）。
两种处理方式：

**方式 1**：CIGAR 之外额外存储 X/I 碱基
```
packed_data = flate2( Vec<CigarOp> + Vec<u8>(X/I 碱基) )
```
解压时按 CIGAR 遍历，遇到 X/I 从碱基流中取。

**方式 2**：用 CIGAR 的 `M` 操作替代 `=/X`
- `M`（match/mismatch）操作不区分 match 和 mismatch
- 样本碱基全部存储（match 区域也存），仅用 CIGAR 定位 indel
- 压缩率差（match 区域冗余存储），不推荐

**推荐方式 1**：CIGAR ops + X/I 碱基流。存储格式：
```
packed_data = flate2(
    u32 cigar_op_count,
    [CigarOp; cigar_op_count],
    u32 x_i_base_count,
    [u8; x_i_base_count]  // X/I 操作的碱基（编码见下）
)
```

> **N 碱基编码**：pbit 支持 ACGTN（[create.rs:21-22](file:///Volumes/ExtHome/Scripts/pgr/src/cmd_pgr/pbit/create.rs#L21-L22)），
> 2-bit 编码（A=0,C=1,G=2,T=3）无法表示 N。采用与 LZ-diff 一致的 5 状态编码
> （[pbit.md:322](file:///Volumes/ExtHome/Scripts/pgr/notes/design/pbit.md#L322) `A=0,C=1,G=2,T=3,N=4`，
> 需 3 bit）。实现上每碱基存 1 字节（ASCII，简单），或打包为 3-bit 流（紧凑但复杂）。
> 默认采用 1 字节/碱基（ASCII），3-bit 打包作为 Phase 8f 优化候选。

## 9. LZ-diff vs CIGAR 压缩率对比

| 场景 | LZ-diff | CIGAR (bit-packed) | 胜出 |
|------|---------|---------------------|------|
| 纯 match（高相似） | `!` back-ref + match-to-end，极紧凑 | `[=4096]`，1 op = 4B | LZ-diff 略优 |
| 稀疏 SNP（每 100bp 1 SNP） | literal + match，~每 SNP 2B | `=99 X1` ≈ 每 SNP 8B | LZ-diff 优 |
| 密集 SNP | 大量 literal | `Xn`，每段几个 op | CIGAR 优 |
| 含 indel | 位置错位→退化为 literal | `In Dn`，直接描述 | CIGAR 显著优 |
| 重排/转位 | 无法处理（参考选择错误） | PAF 坐标直接描述 | CIGAR 唯一可行 |

> **结论**：CIGAR 在含 indel/重排的场景下显著优于 LZ-diff；在纯 SNP 场景下 LZ-diff 略优。
> 混合模式（方案 D）可让用户按需选择。实际压缩率取决于数据特征，需基准测试验证。
>
> **CIGAR 命中率 caveat**：上述对比假设段被单条 PAF alignment 完整覆盖。实际中若一段跨多条
> alignment 的衔接点（如 alignment1 覆盖 [0,3000)、alignment2 覆盖 [3000,6000)，段为 [1000,5000)），
> 则该段被视为部分覆盖，按 §11 决策 3 整段回退 LZ-diff。高连续性 PAF（minimap2 `--paf-no-hit`、
> wfmash）少见；低连续性 PAF（mashmap 多 mapping、分段输出）会显著降低 CIGAR 命中率。
> Phase 8f 基准测试应包含"PAF 连续性 vs CIGAR 命中率"维度。

## 10. 实施阶段（草案）

### Phase 8a: CIGAR 编解码基础设施
- `libs/pbit/cigar_delta.rs`：CIGAR ↔ bit-packed 存储 + X/I 碱基流编码
- 复用 `libs/paf/cigar.rs` 的 `CigarOp`
- 公共 API（与 `segment.rs` 的 `Segment::add/get` 对应）：
  ```rust
  /// Pack CIGAR ops + X/I bases into a flate2-compressed byte buffer.
  pub fn pack_cigar(ops: &[CigarOp], xi_bases: &[u8]) -> Vec<u8>;
  /// Unpack a flate2-compressed buffer into (CIGAR ops, X/I bases).
  pub fn unpack_cigar(packed: &[u8]) -> anyhow::Result<(Vec<CigarOp>, Vec<u8>)>;
  /// Apply CIGAR to a reference slice, consuming X/I bases, producing sample seq.
  /// Used by Decompressor. Simplified variant of build_pairwise_block: produces
  /// raw sample sequence (no '-' gap insertion, no coordinate trimming), logic
  /// follows build_maf_block's =/X/M/I/D branches.
  pub fn apply_cigar(ref_seq: &[u8], ops: &[CigarOp], xi_bases: &[u8]) -> anyhow::Result<Vec<u8>>;
  ```
- 单元测试：CIGAR 往返、X/I 碱基流往返、空 CIGAR、纯 match、含 N 的 X/I 碱基

### Phase 8b: 格式扩展
- `format.rs`：`DeltaMeta`/`DeltaEntry` 新增 `encoding` 字段（10 字节头）；
  `SegmentDesc` 固定大小存储（新增 ref_start/ref_end，LZ-diff 段填 0，见 §5.2）
- Header 版本号升级（1000 → 1001），不保留旧格式读取路径

### Phase 8c: 压缩端（PAF 驱动）
- `compressor.rs`：新增 `append_sample_with_paf` 方法
  ```rust
  /// Append a sample using PAF-driven CIGAR encoding. Segments covered by PAF
  /// alignments are CIGAR-encoded; uncovered segments fall back to LZ-diff.
  pub fn append_sample_with_paf(
      &mut self,
      sample_name: &str,
      fasta_path: &str,
      paf_path: &str,
  ) -> anyhow::Result<()>;
  ```
- PAF → query-side 区间树构建（见 §7 "PAF query-side 索引"）→ 段级 CIGAR 提取 → CIGAR delta 存储
- 未覆盖段回退 LZ-diff

### Phase 8d: 解压端
- `decompressor.rs`：`get_contig` / `get_sample` 支持 CIGAR 解码
- `decode_delta` 签名扩展为接收 `&SegmentDesc`（CIGAR 段需 ref_start/ref_end，见 §8 说明）
- 复用 `libs/paf/msa_build.rs` 的 CIGAR 应用逻辑

### Phase 8e: CLI 集成
- `create.rs` / `append.rs`：新增 `--paf` 参数；`--name` TSV 扩展为可选三列（见 §6.3）
- 测试场景：
  - PAF 驱动往返（`=/X/I/D` CIGAR，`+` 链）
  - `-` 链 PAF 往返（验证 RC 语义：X/I 碱基从 RC(sample) 提取，解压后正向还原）
  - `M` 操作拆分（模拟 minimap2 未加 `--eqx` 的 CIGAR，验证 `M` → `=/X` 拆分正确性）
  - 混合模式（同归档内部分样本 CIGAR、部分 LZ-diff）
  - 未覆盖段回退（PAF 未覆盖的样本段走 LZ-diff）
  - CIGAR target 投影跨参考段边界 → 回退 LZ-diff（§7 决策 3c）
  - 重排场景（样本段比对到不同参考 contig）
  - PAF 无 CIGAR 标签 → 全部回退（决策 7）
  - PAF 文件为空 → 全部回退（决策 7）
  - PAF 单记录解析错误 → 跳过该记录 + warn（决策 8）
  - `--name` 三列 TSV（CIGAR + LZ-diff 混用）
  - `--paf` 与 `--name` 互斥校验
  - `-i`/`--paf` 数量不匹配 → 报错

### Phase 8f: 基准
- 压缩率对比：LZ-diff vs CIGAR vs 混合（不同数据特征）
- 压缩速度对比（CIGAR 省去 LZ-diff 计算，预期更快）
- 解压速度对比
- **PAF 连续性 vs CIGAR 命中率**：用 minimap2/wfmash/mashmap 分别生成 PAF，统计 CIGAR 模式
  命中段比例（部分覆盖回退 LZ-diff 的段比例）
- **优化候选验证**：SoA 布局（ops/lens 分离）对 flate2 压缩率的提升；3-bit X/I 碱基打包
  对空间的节省

## 11. 决策记录

1. **PAF 粒度**：采用 per-sample PAF，`--paf` 与 `-i` 一一对应。

2. **多 alignment 重叠**：选最佳（identity 最高 + 覆盖度最大），不合并。

3. **部分覆盖**：整段回退 LZ-diff，不拆段。触发回退的三种情况：
   (a) 最佳 alignment 未完整覆盖该段（段有未覆盖的 flank）；
   (b) 该段跨多条 alignment 的衔接点（需合并 CIGAR，不合并则覆盖不全）；
   (c) CIGAR 的 target 投影区间跨越参考段边界（`target_start / segment_size !=
       (target_end - 1) / segment_size`），单个 `ref_group_id` 无法覆盖完整 target 区间。
   判定标准：选最佳 alignment 后，检查其 query 覆盖区间是否完整包含 `[seg_start, seg_end)`，
   且 target 投影区间不跨参考段边界。任一条件不满足 → 整段回退 LZ-diff。

4. **`M` 操作处理**：`pgr maf to-paf` 已改为在源头区分 `=/X`（[cigar.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/cigar.rs)
   的 `cigar_from_alignment` 对比 MAF 对齐碱基，case-insensitive）。因此本项目通过
   UCSC pipeline → `maf to-paf` 路径生成的 PAF 已经是 `=/X`。pbit 压缩时仍保留 `M` → `=/X`
   拆分逻辑作为 fallback，处理 minimap2（未用 `--eqx`）、fastga 等直接输出的 `M` CIGAR。
   **注意**：`M` → `=/X` 拆分需读取样本序列（比对 M 区域的样本碱基与参考碱基）。压缩端有样本
   序列（从 FASTA 读取），可行；解压端无需此逻辑（CIGAR 已是 `=/X`）。`M` → `=/X` 拆分仅在
   CIGAR 含 `M` 操作时触发（如 minimap2 未加 `--eqx`）；纯 `=/X` CIGAR（来自 `maf to-paf`、
   minimap2 `--eqx`）无此开销。

5. **版本兼容**：无需考虑。本项目尚未正式发布，直接升级 Header 版本号（1000 → 1001），
   不保留旧格式读取路径。

6. **query-side 索引方案**：pbit 自建 query-side 区间树，不复用 `PafIndex.reverse_trees`。
   理由：(1) `reverse_trees` 仅覆盖 `+` 链记录，`-` 链完全缺失；(2) mirror entry 的元数据
   被交换（query↔target 角色互换、strand 强制为 `+`、CIGAR 被 reverse），不适合 pbit 需要
   原始 PafMetadata 的场景；(3) `reverse_trees` 无公开查询接口。pbit 在 `libs/pbit/` 内构建
   独立的 `BasicCOITree<PafMetadata, u32>`（key = query_id，interval = query 坐标），存储原始
   未交换的 PafMetadata。详见 §7 "PAF query-side 索引"。

7. **PAF 无 CIGAR 处理**：若 PAF 记录无 `cg:Z:` CIGAR 标签（`extract_cigar` 返回空 Vec），
   该记录不插入 query-side 索引（视为未覆盖，相关样本段回退 LZ-diff）。若整个 PAF 文件无
   任何 CIGAR，所有样本段回退 LZ-diff（等价于无 PAF 输入）。实现时不报错，仅 log 警告。

8. **PAF 解析错误处理**：单条 PAF 记录解析错误（格式错误、坐标非法、字段缺失等）时跳过
   该记录并 `log::warn`，相关样本段视为未覆盖，回退 LZ-diff（与决策 7 一致）。若整个 PAF
   文件无法打开或非 PAF 格式（如文件不存在、首行解析失败且无有效记录），返回错误终止压缩
   （`anyhow::bail!`），避免用户误以为 PAF 已生效但实际全部回退。
