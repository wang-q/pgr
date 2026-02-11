# KaKs_Calculator3.0 与 PAML 项目分析报告

本以此文档详细分析 `KaKs_Calculator3.0` 与 `paml-master` 的代码结构、核心功能及异同点，为 `pgr` 项目的潜在集成与参考提供依据。

## 1. KaKs_Calculator3.0 分析

### 1.1 概览
*   **定位**: 专注于计算非同义替换率 (Ka)、同义替换率 (Ks) 及其比值 (Ka/Ks) 的工具，支持多种经典算法及模型选择 (Model Selection) 和模型平均 (Model Averaging) 方法。
*   **语言**: C++
*   **架构**: 面向对象设计，核心逻辑封装在 `KAKS` 类中。

### 1.2 核心代码结构
*   **类继承体系**:
    *   `Base`: 基类，负责通用的序列处理和基础功能。
    *   `KAKS`: 核心派生类 (定义于 `src/KaKs.h`)，继承自 `Base`。它充当了算法的调度器。
*   **算法实现**:
    *   采用了策略模式的变体，不同的计算方法被封装为独立的模块或函数调用。
    *   **支持算法列表** (见 `src/KaKs.h`):
        *   `NG86`: Nei-Gojobori (1986) - *对应 `start_NG86()`*
        *   `LWL85`: Li-Wu-Luo (1985) - *对应 `start_LWL85()`*
        *   `LPB93`: Li-Pamilo-Bianchi (1993)
        *   `GY94`: Goldman-Yang (1994) - *最大似然法基础*
        *   `YN00`: Yang-Nielsen (2000) - *对应 `start_YN00()`*
        *   `MYN`: Modified YN00 - *修正的 YN00 方法*
        *   `MSMA`: Model Selection and Model Averaging - *通过 AICc 进行模型加权*
*   **辅助功能**:
    *   `AXTConvertor.cpp`: 提供了将多种格式 (Clustal, Msf, Nexus, Phylip, Pir) 转换为 AXT 格式的功能。
    *   `ConcatenatePairs.cpp`: 处理序列对的连接。

### 1.3 特点
*   **易用性**: 提供了 Windows GUI 支持 (`result4Win` 相关接口)。
*   **灵活性**: 允许用户指定特定方法或使用 MSMA 自动选择最佳模型。
*   **输出丰富**: 除了 Ka/Ks，还输出 AICc, Fisher P-Value, 详细的位点计数 (S-Sites, N-Sites) 等。

---

## 2. PAML (paml-master) 分析

### 2.1 概览
*   **定位**: **Phylogenetic Analysis by Maximum Likelihood**。这是分子进化领域的权威工具包，不仅限于 Ka/Ks 计算，更专注于复杂的系统发育分析、进化假设检验和祖先序列重建。
*   **语言**: ANSI C
*   **架构**: 过程式编程，包含多个独立的可执行程序，共享底层数学库。

### 2.2 核心组件
PAML 是一个工具箱，主要包含以下核心程序：
1.  **`codeml`** (`src/codeml.c`):
    *   **核心引擎**: PAML 中最复杂的程序。
    *   **功能**: 用于密码子 (Codon) 和氨基酸 (AA) 序列分析。
    *   **模型**: 支持 Branch models (检测特定分支的选择压力), Site models (检测特定位点的正选择), Branch-site models。
    *   **算法**: 实现了严格的最大似然估计 (MLE)，使用数值优化算法 (如 `ming2`, `Newton`) 求解参数。
2.  **`yn00`** (`src/yn00.c`):
    *   **功能**: 专门实现 Yang & Nielsen (2000) 方法的独立程序。
    *   **对比**: 相比 `codeml` 的迭代优化，`yn00` 是一种近似计数方法 (Counting method)，计算速度更快，但不如 ML 方法严谨。
3.  **`baseml`** (`src/baseml.c`):
    *   **功能**: 用于核苷酸序列分析 (GTR, HKY85 等模型)。
4.  **`mcmctree`**:
    *   **功能**: 贝叶斯分歧时间估算 (Bayesian dating)。

### 2.3 代码风格与实现
*   **全局状态**: 大量使用全局变量 (`com` 结构体, `noisy`, `GeneticCode` 等) 来管理配置和状态。
*   **数学库**: `tools.c` 包含了大量的统计分布函数、矩阵运算和数值优化例程。
*   **数据结构**: `paml.h` 定义了通用的 `DataType`, `ReadSeq` 等接口。
*   **控制文件**: 极其依赖 `.ctl` 文件进行参数配置，而非命令行参数。

---

## 3. 对比与总结

### 3.1 功能重叠点
*   **YN00 算法**:
    *   **PAML**: 原创实现 (`src/yn00.c`)。
    *   **KaKs_Calculator**: 复现/集成实现 (`src/YN00.cpp` & `start_YN00()`)。
*   **GY94 模型**:
    *   两者均支持，但 PAML 的 `codeml` 提供了基于树的完整 GY94 实现，而 KaKs_Calculator 主要关注成对比较。

### 3.2 关键差异
*   **核心目标**
    *   **KaKs_Calculator 3.0**: 快速计算成对 Ka/Ks，模型比较
    *   **PAML (codeml/yn00)**: 复杂的进化假设检验，树的构建与参数估算
*   **编程范式**
    *   **KaKs_Calculator 3.0**: C++ (OOP)，模块化较好
    *   **PAML (codeml/yn00)**: C (Procedural)，数学逻辑密集
*   **输入格式**
    *   **KaKs_Calculator 3.0**: AXT, Fasta (需转换)
    *   **PAML (codeml/yn00)**: PHYLIP, NEXUS (需严格格式化)
*   **配置方式**
    *   **KaKs_Calculator 3.0**: 命令行参数或 GUI
    *   **PAML (codeml/yn00)**: 详细的 `.ctl` 配置文件
*   **适用场景**
    *   **KaKs_Calculator 3.0**: 大规模成对序列扫描 (Pairwise Scan)
    *   **PAML (codeml/yn00)**: 精细的单基因家族分析，检测正选择位点
*   **拓展性**
    *   **KaKs_Calculator 3.0**: 易于集成新算法 (Strategy Pattern)
    *   **PAML (codeml/yn00)**: 难以修改核心逻辑，但功能极其强大

### 3.3 对 pgr 项目的启示与集成建议

1.  **算法移植**:
    *   若 `pgr` 需要实现内置的 Ka/Ks 计算，**KaKs_Calculator 的 C++ 源码**是更好的参考对象，其类结构 (`KaKs.h`) 更易于移植到 Rust (Trait/Struct)。
    *   PAML 的 `tools.c` 提供了极其宝贵的数值计算参考 (如 Gamma 分布、SVD 分解)，在实现复杂统计功能时可供借鉴。

2.  **流程集成**:
    *   **轻量级任务**: 对于 `pgr` 的常规 QC 或简单统计，可以直接参考 `KaKs_Calculator` 的逻辑实现快速计算。
    *   **重量级任务**: 对于需要发表级别的进化分析 (如 `pgr` 产生的多序列比对)，建议生成 PAML 兼容的输入文件 (PHYLIP + CTL)，调用外部 `codeml` 二进制文件，而不是尝试重写 PAML 的 ML 引擎。

3.  **格式支持**:
    *   KaKs_Calculator 对 AXT 的原生支持与 `pgr` (基于 AXT/MAF) 高度契合。`src/AXTConvertor.cpp` 中的逻辑可用于验证 `pgr` 的格式转换模块。

## 4. 迁移与集成计划 (Migration Plan)

### 阶段一：格式兼容与数据准备 (Format Compatibility)
*   **目标**: 确保 `pgr` 能生成 KaKs_Calculator 和 PAML 所需的高质量输入文件。
*   **任务**:
    1.  **增强 AXT 支持**: 完善 `pgr net to-axt` 和 `pgr axt` 模块，确保生成的 AXT 格式完全兼容 KaKs_Calculator。
    2.  **新增 PHYLIP 导出**: 实现 `pgr fas to-phylip`，支持 PAML 所需的严格交错/顺序格式 (Interleaved/Sequential)，并处理长序列名截断问题。
    3.  **CDS 提取与比对**: 开发或集成功能，将基因组比对 (MAF/AXT) 映射回 CDS 坐标，生成密码子比对 (Codon Alignment)。

### 阶段二：核心算法移植 (Native Implementation)
*   **目标**: 在 Rust 中原生实现轻量级算法，减少外部依赖，直接集成到 `pgr` 二进制中。
*   **任务**:
    1.  **基础统计**: 移植 `KaKs_Calculator` 中的 `NG86` 和 `LWL85` 算法。这些算法计算量小，易于并行化。
    2.  **YN00 移植**: 参考 `src/YN00.cpp` 和 `src/yn00.c`，在 Rust 中实现 YN00 方法。
    3.  **数学库集成**: 将 PAML `tools.c` 中的关键概率分布函数 (如 Gamma, Chi2) 移植到 Rust 数学模块中，或使用 Rust 生态中的现成库 (如 `statrs`)。

### 阶段三：高级功能 (Advanced Features)
*   **目标**: 支持基因组级扫描。
*   **任务**:
    1.  **滑窗分析**: 结合 `pgr fas window` 和内置 Ka/Ks 计算，实现全基因组滑窗选择压力扫描。
    2.  **并行化**: 利用 Rust 的 `rayon` 库，对大规模基因家族数据进行并行 Ka/Ks 计算。
