# chain

`pgr chain` 模块提供了一套用于操作和处理 **UCSC Chain** 格式比对文件的工具。这些工具是构建全基因组比对（Whole Genome Alignment）流水线（如 `pgr pl ucsc`）的核心组件。

## 核心定位

- **定位**：Chain 格式的高级处理与过滤工具。
- **输入**：Chain 格式文件（文本）。
- **输出**：处理后的 Chain 文件或 Net 文件。
- **互补**：
  - 上游：`lastz` -> `axtChain` (或 `pgr axt to-chain` 计划中) 生成 Chain。
  - 下游：`pgr net` (生成 Net), `pgr maf` (生成 MAF)。

## 子命令详解

### 1. `pgr chain sort`: 排序
对 Chain 文件按得分（Score）降序排列。

- **用途**: `chainPreNet` 和 `chainNet` 等下游工具通常要求输入是按得分排序的。
- **参数**:
  - `files`: 输入 Chain 文件列表。
  - `--input-list`: 包含文件路径列表的文件（用于批量处理）。
  - `--save-id`: 保留原始 Chain ID（默认会从 1 开始重新编号）。

### 2. `pgr chain split`: 拆分
将 Chain 文件按 Target 或 Query 序列名称拆分到不同文件中。

- **用途**: 并行处理或将大文件分割为按染色体组织的文件结构。
- **参数**:
  - `out_dir`: 输出目录。
  - `-q`: 按 Query 序列拆分（默认按 Target）。
  - `--lump <N>`: 将结果合并为 N 个文件（通过哈希分桶），避免产生过多小文件。

### 3. `pgr chain stitch`: 缝合
将具有相同 ID 的 Chain 片段缝合为一个完整的 Chain。

- **用途**: 修复由于并行处理或文件分割导致的同一 Chain 被打断的问题。
- **逻辑**: 自动检查 Target/Query 坐标和 Strand 的一致性，合并 Block 并更新 Header。

### 4. `pgr chain anti-repeat`: 去重复与低复杂度过滤
过滤掉主要由重复序列或低复杂度区域组成的 Chain。

- **用途**: 提高比对质量，去除无生物学意义的假阳性比对。
- **机制**:
  - **Degeneracy Filter**: 检查比对是否主要由低复杂度序列（如 `ATATAT...`）组成。
  - **Repeat Filter**: 检查比对是否落在软屏蔽（Soft-masked，小写字母）区域。
- **参数**:
  - `--target`/`--query`: 对应的 2bit 文件（必需，用于获取序列内容）。
  - `--min-score`: 过滤后的最低得分阈值。
  - `--no-check-score`: 高于此得分的 Chain 将跳过检查（视为可信）。

### 5. `pgr chain pre-net`: Net 前预处理
移除那些被更高得分 Chain 完全覆盖、没有机会形成 Net 的 Chain。

- **用途**: 显著减小 Chain 文件大小，加速后续的 `chainNet` 步骤。
- **机制**: 使用位图（BitMap）跟踪 Target 和 Query 基因组的覆盖情况。优先处理高分 Chain，若低分 Chain 的所有 Block 都已被高分 Chain 覆盖，则丢弃。
- **参数**:
  - `target_sizes`/`query_sizes`: 染色体大小文件。
  - `--pad`: 在 Block 周围添加 Padding，减少碎片。
  - `--incl-hap`: 是否包含 Haplotype 序列（`_hap`, `_alt`）。

### 6. `pgr chain net`: 生成 Net
将 Chain 文件转换为 Net 格式（Syntenic Nets）。

- **用途**: Net 格式表示了基因组之间的高级对应关系，能够区分 Orthologs（直系同源）和 Paralogs（旁系同源），并处理倒位（Inversions）和易位（Translocations）。
- **输出**: 同时生成 Target-referenced Net (`target_net`) 和 Query-referenced Net (`query_net`)。
- **参数**:
  - `target_sizes`/`query_sizes`: 染色体大小文件。
  - `--min-space`: 最小 Gap 大小（小于此值的 Gap 将被填充）。
  - `--min-score`: 最小 Chain 得分阈值。
  - `--incl-hap`: 是否包含 Haplotype 序列。

## 典型工作流 (UCSC Pipeline)

```bash
# 1. 排序 (Sort)
pgr chain sort raw.chain > sorted.chain

# 2. 预处理 (PreNet) - 去除被高分比对覆盖的冗余
pgr chain pre-net sorted.chain t.sizes q.sizes pre.chain

# 3. 生成 Net (Net)
pgr chain net pre.chain t.sizes q.sizes t.net q.net

# 4. 添加 Synteny 信息 (可选，通常配合 pgr net syntenic)
# ...
```
