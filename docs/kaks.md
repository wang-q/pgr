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

### 1.3 模型选择策略 (Model Selection)
*   `KaKs_Calculator` 的 `MS` 策略会遍历 **14 种核苷酸替换模型** (JC, F81, K2P, HKY, TNEF, TN, K3P, K3PUF, TIMEF, TIM, TVMEF, TVM, SYM, GTR)。
*   **选择标准**: 使用 **AICc** (Corrected Akaike Information Criterion, 校正赤池信息量准则) 评估每个模型，选取 AICc 值最小的模型作为最佳模型。

### 1.4 特点
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
    *   **算法**: 实现了严格的最大似然估计 (MLE)，使用数值优化算法 (如 `ming2`, `Newton`) 求解参数。
2.  **`yn00`** (`src/yn00.c`):
    *   **功能**: 专门实现 Yang & Nielsen (2000) 方法的独立程序。
    *   **对比**: 相比 `codeml` 的迭代优化，`yn00` 是一种近似计数方法 (Counting method)，计算速度更快，但不如 ML 方法严谨。
3.  **`baseml`** (`src/baseml.c`):
    *   **功能**: 用于核苷酸序列分析 (GTR, HKY85 等模型)。
4.  **`mcmctree`**:
    *   **功能**: 贝叶斯分歧时间估算 (Bayesian dating)。

### 2.3 模型体系与假设检验
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

### 2.4 代码风格与实现
*   **全局状态**: 大量使用全局变量 (`com` 结构体, `noisy`, `GeneticCode` 等) 来管理配置和状态。
*   **数学库**: `tools.c` 包含了大量的统计分布函数、矩阵运算和数值优化例程。
*   **数据结构**: `paml.h` 定义了通用的 `DataType`, `ReadSeq` 等接口。
*   **控制文件**: 极其依赖 `.ctl` 文件进行参数配置，而非命令行参数。

---

## 6. Bio-Tools-Phylo-PAML (BioPerl 模块) 分析

### 6.1 概览
*   **定位**: BioPerl 的一部分，专门用于解析 PAML (codeml, baseml, yn00) 的输出结果。
*   **语言**: Perl
*   **价值**: 提供了成熟的 PAML 输出解析逻辑，是 `pgr` 实现外部 PAML 调用后结果处理的重要参考。

### 6.2 核心功能
*   **解析器 (`Bio::Tools::Phylo::PAML`)**:
    *   支持解析 `codeml.mlc` 等主输出文件。
    *   提取 NG86 和 ML 方法生成的 dN/dS 矩阵。
    *   提取密码子使用表 (Codon Usage)。
    *   提取带有分支参数 (如 omega, dN, dS) 的系统发育树。
    *   提取模型参数 (Kappa, p0, p1, w0, w1 等)。

---

## 7. 对比与总结

### 7.1 功能重叠点
*   **YN00 算法**:
    *   **PAML**: 原创实现 (`src/yn00.c`)。
    *   **KaKs_Calculator**: 复现/集成实现 (`src/YN00.cpp` & `start_YN00()`)。
*   **GY94 模型**:
    *   两者均支持，但 PAML 的 `codeml` 提供了基于树的完整 GY94 实现，而 KaKs_Calculator 主要关注成对比较。

### 7.2 关键差异
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

### 7.3 对 pgr 的借鉴意义
*   **输出解析逻辑参考**: Bio-Tools-Phylo-PAML 提供了 PAML 输出解析的正则逻辑 (Regex)，可移植到 Rust 中用于 `pgr` 读取 `codeml` 结果（如矩阵提取、LRT 参数提取）。
*   **算法移植参考**: `KaKs_Calculator` 的 C++ 源码适合作为 Rust 原生实现的蓝本；PAML 的 `tools.c` 适合作为数值计算的参考。
*   **流程集成建议**:
    *   **轻量级任务**: 内置实现 NG86/LWL85/YN00。
    *   **重量级任务**: 生成 PHYLIP+CTL 调用外部 `codeml`，并使用移植的解析器读取结果。

---

## 5. 其他主流 Ka/Ks 计算软件

除了 KaKs_Calculator 和 PAML，目前生物信息学领域还有以下几款主流软件可用于 Ka/Ks 分析，各有其侧重：

### 5.1 HyPhy (Hypothesis Testing using Phylogenies)
*   **定位**: 现代化的系统发育假设检验工具，常被视为 PAML 的高性能替代者。
*   **特点**:
    *   **高性能**: 支持多线程 (OpenMP) 和 MPI 并行，甚至部分支持 GPU 加速，处理大数据集比 PAML 快得多。
    *   **方法丰富**: 独创了多种检测正选择的方法，如 **FEL** (Fixed Effects Likelihood), **SLAC** (Single-Likelihood Ancestor Counting), **FUBAR** (Fast, Unconstrained Bayesian AppRoximation), **MEME** (Mixed Effects Model of Evolution)。
    *   **Datamonkey**: 其 Web 版本非常流行。
*   **适用场景**: 大规模基因组数据的正选择扫描，以及需要利用复杂统计模型检测特定位点/分支选择压力的场景。

### 5.2 MEGA (Molecular Evolutionary Genetics Analysis)
*   **定位**: 综合性分子进化分析软件，以 GUI (图形界面) 著称。
*   **特点**:
    *   **易用性**: 集成了比对 (Muscle/Clustal)、建树 (NJ/ML/MP) 和选择压力分析 (Z-test for Selection)。
    *   **算法**: 主要支持基于 Nei-Gojobori 等计数法的 Z 检验，用于判断 Ka/Ks 是否显著大于/小于 1。
*   **适用场景**: 教学、初学者使用，或小规模数据的快速可视化分析。不适合集成到自动化流程中。

### 5.3 paPAML (Parallel PAML)
*   **定位**: PAML 的并行化封装工具。
*   **特点**: 使用 Perl 脚本将大的任务拆分，在集群上并行运行 `codeml` 或 `HyPhy`。
*   **适用场景**: 必须使用 PAML 及其模型结果，但受限于单线程速度的大型项目。

### 5.4 FastCodeML / Godon
*   **定位**: `codeml` 的高性能重写版本。
*   **特点**: 旨在复现 `codeml` 的分支-位点模型 (Branch-Site Model)，但通过优化数学库和并行化大幅提升速度。
*   **适用场景**: 专门针对分支-位点模型的全基因组扫描。

### 5.5 为什么 pgr 仍需自研/移植？
尽管有上述工具，`pgr` 选择原生移植仍有独特价值：
*   **零依赖 (Zero Dependency)**: 用户无需安装 PAML/HyPhy/Perl 等复杂的外部环境，`pgr` 单一二进制文件即可运行。
*   **格式无缝对接**: 直接支持 pgr/UCSC 的 AXT/MAF/2bit 格式，无需中间转换 (如转为 PHYLIP)。
*   **极速**: 针对 NG86/LWL85/YN00 等常用算法，Rust 实现可以比脚本调用外部工具快几个数量级，且比 PAML 的启动开销小得多。

---

## 8. 迁移与集成计划 (Migration Plan)

为了在 `pgr` 中集成 Ka/Ks 计算功能，我们将直接基于标准数据结构分阶段实现算法移植。

### 阶段一：数据准备 (Data Preparation)
*   **目标**: 直接使用 `pgr` 内部的标准数据结构 (如 UCSC AXT/MAF) 为算法提供输入，不兼容 KaKs_Calculator 的私有格式。
*   **任务**:
    1.  **AXT/MAF 读取**: 利用 `pgr` 现有的解析器读取标准基因组比对文件。
    2.  **序列预处理**: 从比对记录中提取 CDS 序列，处理 gap 和反向互补，直接转换为算法所需的内存数据结构。

### 阶段二：核心算法移植 (Native Implementation)
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

    *   **算法选择依据 (Algorithm Selection Rationale)**:
        初期仅移植以下三种代表性算法，覆盖了从“极速粗略”到“高精近似”的不同需求：
        1.  **NG86 (Nei-Gojobori, 1986)**:
            *   **定位**: 基准算法 (Baseline)。
            *   **理由**: 计算最简单、速度最快，适合大规模初筛，是所有 Ka/Ks 分析的“标尺”。
        2.  **LWL85 (Li-Wu-Luo, 1985)**:
            *   **定位**: 进阶近似法。
            *   **理由**: 区分了转换 (Transition) 和颠换 (Transversion)，比 NG86 更符合生物学实际。
        3.  **YN00 (Yang-Nielsen, 2000)**:
            *   **定位**: 近似法中的“黄金标准”。
            *   **理由**: 考虑了**密码子使用偏好 (Codon Usage Bias)**，准确性接近 ML，但速度远快于 ML，是目前性价比最高的算法。

        *注意*: 暂不移植 GY94 等最大似然法 (Maximum Likelihood)，因为其计算成本极高。对于需要 ML 级精度的用户，建议通过调用外部 PAML `codeml` 来实现。

    *   **任务**:
        1.  **基础架构**: 建立上述目录结构，定义 Trait 和 Result。
        2.  **基础统计**: 移植 `KaKs_Calculator` 中的 `NG86` 和 `LWL85` 算法。
        3.  **YN00 移植**: 参考 `src/YN00.cpp` 和 `src/yn00.c`，在 Rust 中实现 YN00 方法。
        4.  **数学库集成**: 移植 PAML `tools.c` 中的关键概率分布函数或使用 `statrs`。

### 阶段三：高级功能与优化
*   **目标**: 支持基因组级扫描。
*   **任务**:
    1.  **滑窗分析**: 结合 `pgr fas window` 和内置 Ka/Ks 计算，实现全基因组滑窗选择压力扫描。
    2.  **并行化**: 利用 Rust 的 `rayon` 库，对大规模基因家族数据进行并行 Ka/Ks 计算。
