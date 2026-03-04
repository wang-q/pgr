# pgr mat transform: 相似度到距离的转换方案

`pgr mat transform` 命令用于对矩阵中的数值进行数学变换。

这是将**相似度矩阵 (Similarity Matrix)** 转换为 **距离矩阵 (Distance Matrix)** 的核心工具，也支持归一化和其他数值调整。

## 背景 (Context)

聚类算法（如 UPGMA, NJ, Ward）和多维尺度分析（MDS）通常要求输入 **距离矩阵 (Distance Matrix)** 或 **相异度矩阵 (Dissimilarity Matrix)**，满足：
- $D(x, x) = 0$
- $D(x, y) \ge 0$
- $D(x, y)$ 越小表示越相似

然而，生物信息学上游工具（如 BLAST, MMseqs2, Diamond）或统计分析通常输出 **相似度 (Similarity)**，满足：
- $S(x, x) = Max$ (如 1.0 或 100)
- $S(x, y)$ 越大表示越相似

目前用户需要使用 `awk` 或外部脚本进行转换（例如 `100 - identity`），这不方便且容易出错（如未处理缺失值或自比对）。

## 常见转换模型 (Transformation Models)

我们支持以下几种常见的转换模式：

### 1. 线性反转 (Linear Inversion)
适用于有固定上限的相似度（如 Identity, Percent Similarity）。
$$D = Max - S$$
- **场景**: BLAST Identity (0-100) $\rightarrow$ $D = 100 - S$
- **场景**: Fraction (0-1) $\rightarrow$ $D = 1 - S$

### 2. 归一化线性反转 (Normalized Linear Inversion)
如果 $S$ 没有固定上限（如 Alignment Score），需先归一化。
$$D = 1 - \frac{S(x, y)}{\sqrt{S(x, x) \cdot S(y, y)}}$$
或者简单的:
$$D = 1 - \frac{S(x, y)}{Max(S)}$$

### 3. 对数转换 (Logarithmic)
适用于概率或乘性模型（类似于 Jukes-Cantor 校正）。
$$D = -\ln(S)$$
或者归一化后：
$$D = -\ln(\frac{S(x, y)}{\sqrt{S(x, x) \cdot S(y, y)}})$$
- **场景**: 序列一致性概率 $\rightarrow$ 进化距离

### 4. 倒数转换 (Reciprocal)
$$D = \frac{1}{S} - \frac{1}{Max}$$
- **场景**: 很少见，用于某些物理量的转换。

### 5. 特殊转换
- **Cosine Similarity**: $D = 1 - \cos(\theta)$
- **Correlation**: $D = \sqrt{2(1 - r)}$ 或 $D = 1 - r$

## 用法 (Usage)

```bash
pgr mat transform [OPTIONS] <infile>
```

### 参数 (Arguments)

- `<infile>`: 输入 PHYLIP 矩阵文件。

### 选项 (Options)

- `--op <METHOD>`: 变换操作 (默认: `linear`)。
  - `linear`: $val = val \times scale + offset$
  - `inv-linear`: $val = max - val$
  - `log`: $val = -\ln(val)$ (如果 $val \le 0$ 则设为 0 或 Inf)
  - `exp`: $val = \exp(-val)$
  - `square`: $val = val^2$
  - `sqrt`: $val = \sqrt{val}$
- `--max <FLOAT>`: 用于 `inv-linear` 的最大值 (默认: 1.0)。
- `--scale <FLOAT>`: 用于 `linear` 的缩放因子 (默认: 1.0)。
- `--offset <FLOAT>`: 用于 `linear` 的偏移量 (默认: 0.0)。
- `--normalize`: 是否在变换前基于对角线元素进行归一化 (需矩阵包含对角线数据)。
  - 归一化公式: $x_{norm}(i, j) = \frac{x(i, j)}{\sqrt{x(i, i) \times x(j, j)}}$
- `-o, --outfile <outfile>`: 输出文件名 (默认: stdout)。

## 常见场景 (Examples)

### 1. Identity (0-100) 转 Distance (0-1)

BLAST 等工具输出的 Identity 通常为 0 到 100。
目标公式: $D = (100 - Identity) / 100 = 1 - 0.01 \times Identity$。

使用 `linear` 操作：
```bash
pgr mat transform input.phy --op linear --scale -0.01 --offset 1.0 -o dist.phy
```

或者分两步（先反转再缩放）：
```bash
pgr mat transform input.phy --op inv-linear --max 100 | \
pgr mat transform stdin --op linear --scale 0.01 -o dist.phy
```

### 2. Identity (0-100) 转 Distance (0-100)

仅做反转：$D = 100 - Identity$。

```bash
pgr mat transform input.phy --op inv-linear --max 100 -o dist.phy
```

### 3. Similarity (0-1) 转 Distance (0-1)

标准的线性反转：$D = 1.0 - S$。

```bash
pgr mat transform input.phy --op inv-linear --max 1.0 -o dist.phy
```

### 4. 概率/乘性模型转换 (Log)

将序列一致性概率转换为进化距离（类似 Jukes-Cantor 校正的第一步）。
$D = -\ln(S)$。

```bash
# 假设输入矩阵是对角线为 1.0 的概率矩阵
pgr mat transform input.phy --op log -o dist.phy
```

### 5. 归一化并转换

如果输入的是未归一化的相似度得分（如 Alignment Score），且矩阵包含对角线（自比对得分）。
先归一化为 0-1，再转换为距离。

```bash
# 1. Normalize: S_norm = S_ij / sqrt(S_ii * S_jj)
# 2. Transform: D = 1.0 - S_norm
pgr mat transform raw_scores.phy --normalize --op inv-linear --max 1.0 -o dist.phy
```

## 注意事项

- **对角线处理**:
  - `pgr` 读取矩阵时通常会忽略对角线（设为 0），但 `transform` 命令会尝试保留对角线信息以支持 `--normalize`。
  - 如果输入文件没有对角线信息（如某些 PHYLIP 变体），`--normalize` 将无法正常工作（视为 0）。
- **数值稳定性**:
  - `log` 操作对 0 或负值敏感，程序会将其处理为极大值或 0。
  - 归一化时如果对角线为 0，结果将为 0。
