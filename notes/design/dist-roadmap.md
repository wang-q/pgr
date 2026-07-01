# dist 未来规划

> **实现状态注记**：本文档列出 `pgr dist` 尚未实现的距离度量与架构优化规划。当前已实现 seq/hv/vector 三个子命令。

## 1. Alignment-based Metrics (计划中)

*基于严格比对的距离计算，精度最高但速度较慢。*

- **Kimura 2-Parameter (K2P)**: 区分转换与颠换。
- **Jukes-Cantor (JC69)**: 基础核苷酸替换模型。
- **p-distance**: 简单的 Hamming 距离。

## 2. Scikit-learn/SciPy 兼容性与架构优化

*深度借鉴 `scikit-learn` 的距离计算架构，以提升可扩展性与性能。*

- **DistanceMetric 接口**: 统一所有距离计算（序列、向量、稀疏矩阵）的 API，便于未来扩展新度量。
- **分块计算 (Chunking)**: 引入类似 `sklearn.metrics.pairwise.gen_batches` 的机制，支持计算超大规模数据集的距离矩阵（流式处理，控制内存峰值）。
- **稀疏矩阵支持**: 借鉴 `sklearn` 的 CSR 优化策略，支持稀疏向量的高效距离计算。
- **SIMD 加速**: `pgr` 的基础线性代数库 (`libs::linalg`) 已经利用 `portable_simd` 实现了欧氏距离、点积等核心计算的向量化加速。计划进一步扩展到更多距离度量。

## 3. SciPy 兼容性扩展 (Metrics)

*计划支持更多 SciPy `spatial.distance` 标准度量：*

- **Bray-Curtis**: 生态学常用。
- **Canberra**: 对小数值敏感。
- **Chebyshev**: 切比雪夫距离 ($L_\infty$)。
- **Mahalanobis**: 马氏距离（考虑协方差）。
