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

### 迁移与集成计划 (Migration Plan)

为了在 `pgr` 中集成 Ka/Ks 计算功能，我们将直接基于标准数据结构分阶段实现算法移植。

#### 阶段一：数据准备 (Data Preparation)
*   **目标**: 直接使用 `pgr` 内部的标准数据结构 (如 UCSC AXT/MAF) 为算法提供输入，不兼容 KaKs_Calculator 的私有格式。
*   **任务**:
    1.  **AXT/MAF 读取**: 利用 `pgr` 现有的解析器读取标准基因组比对文件。
    2.  **序列预处理**: 从比对记录中提取 CDS 序列，处理 gap 和反向互补，直接转换为算法所需的内存数据结构。

#### 阶段二：核心算法移植 (Native Implementation)
*   **目标**: 在 Rust 中原生实现轻量级算法，减少外部依赖，直接集成到 `pgr` 二进制中。
*   **库文件结构设计 (Library Architecture)**:
    为了保持 `pgr` 的模块化和高性能，我们将采用类似于 KaKs_Calculator 的策略模式，但使用 Rust 的 Trait 系统来实现。

    ```text
    src/libs/
    └── kaks/
        ├── mod.rs          # 模块入口，定义 Trait 和公共结构
        ├── result.rs       # 定义 KaKsResult 结构体
        ├── methods/        # 各种算法的具体实现
        │   ├── mod.rs
        │   ├── ng86.rs     # NG86 算法
        │   ├── lwl85.rs    # LWL85 算法
        │   └── yn00.rs     # YN00 算法 (最大似然法)
        └── stats/          # 统计学辅助函数
            ├── mod.rs
            └── gamma.rs    # Gamma 分布相关计算
    ```

    *   **核心 Trait 定义 (`src/libs/kaks/mod.rs`)**:
        ```rust
        pub trait KaKsEstimator {
            /// 计算一对序列的 Ka/Ks
            fn calculate(&self, seq1: &[u8], seq2: &[u8]) -> Result<KaKsResult, KaKsError>;
            
            /// 算法名称
            fn name(&self) -> &'static str;
        }
        ```

    *   **算法选择依据 (Algorithm Selection Rationale)**:
        初期仅移植以下三种代表性算法，覆盖了从“极速粗略”到“高精近似”的不同需求：
        1.  **NG86 (Nei-Gojobori, 1986)**:
            *   **定位**: 基准算法 (Baseline)。
            *   **理由**: 计算最简单、速度最快，假设所有替换发生概率相等。适合作为大规模数据的快速初筛工具，也是所有 Ka/Ks 分析的“标尺”。
        2.  **LWL85 (Li-Wu-Luo, 1985)**:
            *   **定位**: 进阶近似法。
            *   **理由**: 区分了转换 (Transition) 和颠换 (Transversion)，比 NG86 更符合生物学实际（转换通常发生频率更高）。在计算量增加不多的情况下提高了准确性。
        3.  **YN00 (Yang-Nielsen, 2000)**:
            *   **定位**: 近似法中的“黄金标准”。
            *   **理由**: 它是 PAML `yn00` 程序的核心算法。除了考虑转换/颠换，还通过迭代方法估算密码子频率，考虑了**密码子使用偏好 (Codon Usage Bias)**。在准确性上接近最大似然法 (ML)，但计算速度远快于 ML，是目前性价比最高的算法。
        
        *注意*: 暂不移植 GY94 等最大似然法 (Maximum Likelihood)，因为其涉及复杂的数值优化 (BFGS)，计算成本极高。对于需要 ML 级精度的用户，建议在后续版本中通过调用外部 PAML `codeml` 来实现。

        *   **关于 KaKs_Calculator 的模型选择 (Model Selection)**:
            `KaKs_Calculator` 的 `MS` 策略会遍历 **14 种核苷酸替换模型** (JC, F81, K2P, HKY, TNEF, TN, K3P, K3PUF, TIMEF, TIM, TVMEF, TVM, SYM, GTR)。
            *   **选择标准**: 使用 **AICc** (Corrected Akaike Information Criterion, 校正赤池信息量准则) 评估每个模型，选取 AICc 值最小的模型作为最佳模型。
            *   **实现代价**: 这需要基于 GY94 框架运行 14 次最大似然参数估计，计算量巨大，因此 `pgr` 初期不直接支持此功能。

    *   **结果结构体 (`src/libs/kaks/result.rs`)**:
        ```rust
        pub struct KaKsResult {
            pub ka: f64,
            pub ks: f64,
            pub ka_ks: f64,
            pub p_value: Option<f64>, // 用于 Fisher Exact Test (NG86/LWL85)
            pub aicc: Option<f64>,    // 用于模型选择 (AICc)
            pub ln_l: Option<f64>,    // 对数似然值 (ML 方法)
        }
        ```

*   **任务**:
    1.  **基础架构**: 建立上述目录结构，定义 Trait 和 Result。
    2.  **基础统计**: 移植 `KaKs_Calculator` 中的 `NG86` 和 `LWL85` 算法。这些算法计算量小，易于并行化。
    3.  **YN00 移植**: 参考 `src/YN00.cpp` 和 `src/yn00.c`，在 Rust 中实现 YN00 方法。
    4.  **数学库集成**: 将 PAML `tools.c` 中的关键概率分布函数 (如 Gamma, Chi2) 移植到 Rust 数学模块中，或使用 Rust 生态中的现成库 (如 `statrs`)。

### PAML (codeml) 分析

*   **定位**: 此时作为 "Gold Standard" (金标准)，主要用于**假设检验 (Hypothesis Testing)** 和复杂的进化分析，而非简单的 Ka/Ks 计算工具。
*   **模型体系**:
    *   **核苷酸模型 (baseml)**: 支持 JC69, K80, F81, F84, HKY85, T92, TN93, REV (GTR) 等几乎所有主流模型。
    *   **密码子模型 (codeml)**: 基于 GY94 (Goldman & Yang 1994) 和 MG94 框架。
        *   **频率模型**: F1x4, F3x4, F61 (Codon Table)。
        *   **位点模型 (NSsites)**: 用于检测正选择位点。
            *   **M0 (One-ratio)**: 所有位点 $\omega$ 相同。
            *   **M1a (NearlyNeutral)**: 只有保守 ($\omega<1$) 和中性 ($\omega=1$) 位点。
            *   **M2a (PositiveSelection)**: 增加正选择类别 ($\omega>1$)。通常与 M1a 比较。
            *   **M7 (Beta)**: $\omega$ 服从 Beta 分布 (0-1)。
            *   **M8 (Beta&w)**: Beta 分布 + 一个正选择类别。通常与 M7 比较。
*   **选择/评估标准**:
    *   **最大似然法 (Maximum Likelihood, ML)**: 计算给定模型下的参数，使得观测数据的似然值 ($lnL$) 最大。
    *   **似然比检验 (Likelihood Ratio Test, LRT)**:
        PAML **不自动选择最佳模型**。用户通常运行一对嵌套模型 (Nested Models，如 M1a vs M2a)，计算 $2\Delta\ell = 2(lnL_1 - lnL_0)$，然后查 $\chi^2$ 分布表来判断复杂模型是否显著优于简单模型（例如检测是否存在正选择）。
    *   **对比**: `KaKs_Calculator` 倾向于使用 AICc **自动挑选** 最佳拟合模型；而 PAML 倾向于让用户通过 LRT **验证科学假设**。

#### 阶段三：高级功能与优化
*   **目标**: 支持基因组级扫描。
*   **任务**:
    1.  **滑窗分析**: 结合 `pgr fas window` 和内置 Ka/Ks 计算，实现全基因组滑窗选择压力扫描。
    2.  **并行化**: 利用 Rust 的 `rayon` 库，对大规模基因家族数据进行并行 Ka/Ks 计算。
