# pbit PAF 驱动模式设计提案

> **状态**：设计提案，尚未实现。基于 [pbit.md](pbit.md) 的现有格式，探讨用 PAF 比对结果驱动
> LZ-diff 参考选择、以及直接存储 CIGAR 替代 LZ-diff 的可行性与方案。

## 1. 动机

### 当前局限

pbit 当前的参考选择是**简化版 minimizer**（[compressor.rs:280-336](file:///Volumes/ExtHome/Scripts/pgr/src/libs/pbit/compressor.rs#L280-L336)）：

- 按 contig 名匹配参考段（`contig_ref_groups[contig_name]`）
- 按段位置索引匹配（`ref_group_ids[seg_idx]`，clamped 到最后一段）
- 首段 k-mer 采样检测方向 + 逐段 delta 大小回退

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
压缩：解析 PAF → 提取每段 CIGAR → bit-packed Vec<CigarOp> → flate2 压缩
解压：flate2 解压 → Vec<CigarOp> → 按 CIGAR 从参考段重建样本段
```

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

### 5.1 DeltaEntry 扩展

```rust
pub struct DeltaEntry {
    pub is_rev_comp: bool,
    pub raw_length: u32,
    pub packed_data: Vec<u8>,
    // 新增：编码类型（0 = LZ-diff, 1 = CIGAR）
    // 存储在 9 字节头部的 reserved 区域或新增 1 字节
    pub encoding: DeltaEncoding,
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
> 本项目尚未正式发布，无需考虑向后兼容，直接升级 Header 版本号（1000 → 1001）。

### 5.2 SegmentDesc 扩展（CIGAR 模式）

CIGAR 模式下，每段需要记录比对坐标（样本段在参考段中的偏移）：

```rust
pub struct SegmentDesc {
    pub ref_group_id: u32,
    pub delta_id: u32,
    // 新增（CIGAR 模式）：
    pub ref_start: u32,   // 该段在参考段中的起始偏移（0-based）
    pub ref_end: u32,     // 该段在参考段中的结束偏移（exclusive）
    // is_rev_comp / raw_length 由 DeltaEntry 提供
}
```

> **为什么需要 ref_start/ref_end**：CIGAR 描述的是"样本段 vs 参考段某区间"的差异。
> 同一参考段可能被多个样本段引用（不同区间），需要 ref_start/ref_end 定位。

### 5.3 CIGAR 的 bit-packed 存储

`Vec<CigarOp>` 直接序列化为 `Vec<u32>`（每个 CigarOp 是一个 u32），然后 flate2 压缩：

```
packed_data = flate2(raw_u32_array)
```

解压时：`flate2_decompress(packed_data) → Vec<u32> → Vec<CigarOp>`

> **空间效率**：一个 4096 bp 的纯 match 段 → `[=4096]` → 1 个 CigarOp = 4 字节 → flate2 后更小。
> 含 40 个 SNP 的段 → `=99 X1` × 40 + `=96` ≈ 82 个 CigarOp = 328 字节 → flate2 后约 100-200 字节。

## 6. CLI 设计

### 6.1 新增 `--paf` 参数

```
pgr pbit create -r ref.fa -i sample.fa --paf sample.paf -o out.pbit
pgr pbit append in.pbit -i sample.fa --paf sample.paf
```

- `--paf` 可选，指定样本比对到参考的 PAF 文件
- 一个 `-i` 对应一个 `--paf`（或支持一个 PAF 含多个样本的比对）
- 省略 `--paf` 时使用当前 minimizer 逻辑（LZ-diff 模式）

### 6.2 参数语义

```
pgr pbit create -r ref.fa \
    -i sample1.fa --paf sample1.paf \    # sample1 用 CIGAR 模式
    -i sample2.fa                        # sample2 用 LZ-diff 模式（无 PAF）
```

> 同一归档可混用两种模式（每样本独立选择）。delta 层的 encoding 字段区分。

## 7. 压缩流程（CIGAR 模式）

```
append_sample_with_paf(sample_name, fasta_path, paf_path):
  1. 解析 PAF → 构建 PafIndex（coitrees 区间树，按 query 坐标索引）
  2. 读取样本 FASTA → 按 segment_size 分段
  3. 对每个样本段 [seg_start, seg_end):
     a. 查 PafIndex：哪些 PAF alignment 覆盖此区间？
     b. 选最佳 alignment（按 identity / 覆盖度）
     c. 提取该段对应的 CIGAR 片段（从 PAF CIGAR 中按坐标截取）
     d. 若 CIGAR 含 `M` 操作：对比样本碱基与参考碱基，拆分为 `=/X`
        （fallback：`pgr maf to-paf` 已在源头区分 `=/X`，此步骤仅处理
        minimap2/fastga 等直接输出的 `M` CIGAR）
     e. 确定 ref_group_id（PAF target → 参考段映射）+ strand + ref_start/ref_end
     f. bit-pack CIGAR → flate2 压缩 → DeltaEntry(encoding=CIGAR)
     g. 记录 SegmentDesc(ref_group_id, delta_id, ref_start, ref_end)
  4. 未覆盖段：回退 LZ-diff（对到同名参考段，当前逻辑）
```

### 关键问题：PAF 坐标→段映射

PAF alignment 的 query_start/end 是连续坐标，pbit 按 segment_size 切段。一条 PAF alignment
可能横跨多段。处理方式：

1. **保持固定段切分**：样本仍按 segment_size 切分
2. **CIGAR 截取**：对每段，从覆盖它的 PAF alignment 的 CIGAR 中截取对应区间
   - 用 `CigarOp::target_delta` / `query_delta` 做坐标投影
   - 一段可能被多条 alignment 覆盖 → 选最佳（identity 最高）
   - 一段可能部分被覆盖 → 覆盖部分用 CIGAR，未覆盖部分用 LZ-diff 或 raw

### 关键问题：未覆盖区域

PAF 没覆盖的样本序列区域：
- **回退 LZ-diff**：对到同名参考段（当前逻辑），encoding=LzDiff
- **存 raw**：直接存储原始序列（无参考），encoding=Raw（需第三种 encoding）
- **推荐**：回退 LZ-diff（复用现有逻辑，格式变化最小）

## 8. 解压流程（CIGAR 模式）

```
get_contig / get_sample (CIGAR 段):
  1. 读取 DeltaEntry → 检查 encoding
  2. 若 encoding == CIGAR:
     a. flate2 解压 packed_data → Vec<u32> → Vec<CigarOp>
     b. 读取参考段（ref_group_id → seek → read_2bit_record）
     c. 按 SegmentDesc.ref_start/ref_end 截取参考段区间
     d. 若 is_rev_comp: 反向互补参考区间
     e. 按 CIGAR 重建样本段:
        - '=' : 从参考取 len 个碱基
        - 'X' : 跳过参考 len 个碱基（样本碱基未知 → 需存储）
        - 'I' : 插入（样本碱基需存储）
        - 'D' : 跳过参考 len 个碱基
     f. 拼接所有段
  3. 若 encoding == LZ-diff: 当前逻辑（LZ-diff 解码）
```

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
    [u8; x_i_base_count]  // X/I 操作的碱基（2-bit 编码可进一步压缩）
)
```

> **优化**：X/I 碱基可用 2-bit 编码（A=0,C=1,G=2,T=3），压缩为原来的 1/4。

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

## 10. 实施阶段（草案）

### Phase 8a: CIGAR 编解码基础设施
- `libs/pbit/cigar_delta.rs`：CIGAR ↔ bit-packed 存储 + X/I 碱基流编码
- 复用 `libs/paf/cigar.rs` 的 `CigarOp`
- 单元测试：CIGAR 往返、X/I 碱基流往返、空 CIGAR、纯 match

### Phase 8b: 格式扩展
- `format.rs`：DeltaEntry 新增 encoding 字段；SegmentDesc 新增 ref_start/ref_end
- Header 版本号升级（1000 → 1001），不保留旧格式读取路径

### Phase 8c: 压缩端（PAF 驱动）
- `compressor.rs`：新增 `append_sample_with_paf` 方法
- PAF → PafIndex 构建 → 段级 CIGAR 提取 → CIGAR delta 存储
- 未覆盖段回退 LZ-diff

### Phase 8d: 解压端
- `decompressor.rs`：`get_contig` / `get_sample` 支持 CIGAR 解码
- 复用 `libs/paf/msa_build.rs` 的 CIGAR 应用逻辑

### Phase 8e: CLI 集成
- `create.rs` / `append.rs`：新增 `--paf` 参数
- 测试：PAF 驱动往返、混合模式、未覆盖回退、重排场景

### Phase 8f: 基准
- 压缩率对比：LZ-diff vs CIGAR vs 混合（不同数据特征）
- 压缩速度对比（CIGAR 省去 LZ-diff 计算，预期更快）
- 解压速度对比

## 11. 决策记录

1. **PAF 粒度**：采用 per-sample PAF，`--paf` 与 `-i` 一一对应。

2. **多 alignment 重叠**：选最佳（identity 最高 + 覆盖度最大），不合并。

3. **部分覆盖**：整段回退 LZ-diff，不拆段。

4. **`M` 操作处理**：`pgr maf to-paf` 已改为在源头区分 `=/X`（[cigar.rs](file:///Volumes/ExtHome/Scripts/pgr/src/libs/paf/cigar.rs)
   的 `cigar_from_alignment` 对比 MAF 对齐碱基，case-insensitive）。因此本项目通过
   UCSC pipeline → `maf to-paf` 路径生成的 PAF 已经是 `=/X`。pbit 压缩时仍保留 `M` → `=/X`
   拆分逻辑作为 fallback，处理 minimap2（未用 `--eqx`）、fastga 等直接输出的 `M` CIGAR。

5. **版本兼容**：无需考虑。本项目尚未正式发布，直接升级 Header 版本号（1000 → 1001），
   不保留旧格式读取路径。
