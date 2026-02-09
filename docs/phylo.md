# 系统发育树内存数据结构对比：Newick Utilities vs Phylotree-rs

本文档对比了 `newick_utils` (C) 和 `phylotree-rs` (Rust) 两个项目在处理系统发育树时采取的不同内存数据结构设计。这两种设计分别代表了经典的指针链接式（Pointer-based）和现代的 Arena/池化（Arena/Vector-based）方法。

## 1. Newick Utilities (C)

`newick_utils` 采用经典的 C 语言指针方式，通过结构体指针将节点连接成树。

### 核心结构体

节点 (`struct rnode`):
定义在 `src/rnode.h`。
```c
struct rnode {
    char *label;                  // 节点名称
    char *edge_length_as_string;  //以此字符串形式存储（保留格式，如无长度）
    double edge_length;           // 数值形式存储
    void *data;                   // 泛型指针，用于挂载应用特定的数据

    struct rnode *parent;         // 父节点指针
    struct rnode *next_sibling;   // 下一个兄弟节点指针
    struct rnode *first_child;    // 第一个子节点指针
    struct rnode *last_child;     // 最后一个子节点指针（优化追加操作）
    
    int child_count;              // 子节点数量
    // ... 其他迭代和状态标记
};
```

树 (`struct rooted_tree`):
定义在 `src/tree.h`。
```c
struct rooted_tree {
    struct rnode *root;            // 根节点指针
    struct llist *nodes_in_order;  // 辅助链表，按后序遍历存储所有节点（便于线性遍历）
    enum tree_type type;           // 树类型（Cladogram/Phylogram）
};
```

### 设计特点
1.  左孩子-右兄弟 (Left-Child, Right-Sibling): 虽然它显式存储了 `first_child` 和 `last_child`，但兄弟节点之间通过 `next_sibling` 链表连接。这使得特定层级的遍历（如“所有子节点”）需要沿着链表跳转。
2.  分散内存 (Dispersed Memory): 每个 `rnode` 都是单独 `malloc` 分配的。这提供了极大的灵活性（树结构可以随意重组，无需移动内存），但可能导致较差的 CPU 缓存局部性（Cache Locality）。
3.  双重存储: `edge_length` 既存字符串又存 double。这是为了精确保留输入 Newick 文件的格式（例如区分 "0.0" 和 "0" 或无长度的情况），这是作为一个通用处理工具（Text Processing）的考量。
4.  侵入性扩展: 通过 `void *data` 指针，用户可以在运行时挂载任意数据结构到节点上，无需修改核心库代码。

---


## 2. Phylotree-rs (Rust)

`phylotree-rs` 采用了 Rust 中处理图和树结构的常见模式：Arena（竞技场） 或 Vector-backed 模式。所有节点存储在一个中心化的向量中，节点间的引用使用整数索引（ID）。

### 核心结构体

节点 (`struct Node`):
定义在 `src/tree/node.rs`。
```rust
pub struct Node {
    pub id: NodeId,               // 自身 ID (usize)
    pub name: Option<String>,     // 节点名称
    pub parent: Option<NodeId>,   // 父节点 ID (索引)
    pub children: Vec<NodeId>,    // 子节点 ID 列表 (动态数组)
    
    pub parent_edge: Option<EdgeLength>, // 父边长度
    // ...
    pub(crate) child_edges: Option<HashMap<NodeId, EdgeLength>>, // 子边长度映射
    pub(crate) subtree_distances: RefCell<Option<VecMap<EdgeLength>>>, // 缓存距离
}
```

树 (`struct Tree`)**:
定义在 `src/tree/tree_impl.rs`。
```rust
pub struct Tree {
    nodes: Vec<Node>,             // 核心存储：所有节点都在这个 Vector 中
    leaf_index: RefCell<Option<Vec<String>>>, // 缓存叶子节点索引
    // ...
}
```

### 设计特点
1.  Arena 内存布局: 所有 `Node` 结构体都紧凑地存储在 `Tree.nodes` 这个 `Vec` 中。这提供了极佳的缓存局部性，且内存分配只发生在 Vector 扩容时，而非每个节点创建时。
2.  索引引用 (Index-based Reference): 关系通过 `NodeId` (即 `usize` 索引) 维护。避免了 Rust 中自引用结构体的所有权（Ownership）和生命周期（Lifetime）噩梦。
3.  邻接表 (Adjacency List): 每个节点直接持有一个 `Vec<NodeId>` 存储所有子节点。访问所有子节点是极其快速的数组遍历，优于链表跳转。
4.  安全性: 避免了悬垂指针（Dangling Pointers）和内存泄漏。只要持有 `Tree` 对象，所有节点 ID 都是有效的（除非被标记删除）。

---

## 3. 对比总结

*   内存模型
    *   Newick Utilities: 分散式 (Malloc per node)
    *   Phylotree-rs: 集中式 (Vector Arena)
*   引用方式
    *   Newick Utilities: 原始指针 (`*rnode`)
    *   Phylotree-rs: 整数索引 (`usize`)
*   子节点存储
    *   Newick Utilities: 链表 (Next Sibling)
    *   Phylotree-rs: 动态数组 (`Vec<NodeId>`)
*   缓存友好性
    *   Newick Utilities: 低 (节点在堆上随机分布)
    *   Phylotree-rs: 高 (节点在内存中连续)
*   遍历性能
    *   Newick Utilities: 较慢 (指针跳转)
    *   Phylotree-rs: 极快 (数组索引)
*   重构代价
    *   Newick Utilities: 低 (仅需修改指针)
    *   Phylotree-rs: 中 (可能涉及 Vector 元素移动/Swap)
*   类型安全
    *   Newick Utilities: 低 (`void*`, 手动管理内存)
    *   Phylotree-rs: 高 (Rust 类型系统保障)
*   扩展性
    *   Newick Utilities: 运行时 (`void*` hook)
    *   Phylotree-rs: 编译时 (需修改结构体或用泛型)

---

## 4. API 设计对比

除了数据结构，两者在 API 设计风格上也大相径庭，反映了 C 和 Rust 语言特性的不同。

### Newick Utilities (C)

Newick Utilities 的 API 更加**过程式**，强调手动管理和灵活性。

*   **遍历 (Traversal)**:
    *   提供了 `struct rnode_iterator` 对象，用于深度优先遍历。
    *   提供了 `get_nodes_in_order(struct rnode *)` 函数，返回一个后序遍历的链表 (`struct llist`)。这是大多数操作（如计算）的首选方式，因为线性遍历链表比递归遍历树更快。
*   **拓扑修改**:
    *   `reroot_tree(...)`: 重新定根，直接修改指针结构。
    *   `prune`: 剪枝操作（虽然主要是 CLI 工具，但底层有对应逻辑）。
*   **克隆 (Cloning)**:
    *   提供了丰富的克隆函数：`clone_tree`, `clone_subtree`, `clone_tree_cond` (带谓词条件的克隆)。这在需要对树进行破坏性修改前非常有用。
*   **查询**:
    *   提供了 `nodes_from_labels` 和基于正则的 `nodes_from_regexp`，这非常强大，允许通过复杂的名称模式匹配节点。
    *   `is_cladogram`: 检查是否为有根树。

### Phylotree-rs (Rust)

Phylotree-rs 的 API 更加**面向对象**和**函数式**，利用 Rust 的特性提供安全性和易用性。

*   **遍历 (Traversal)**:
    *   提供了 `preorder`, `inorder`, `postorder`, `levelorder` 等方法，返回 `Vec<NodeId>`。这利用了 Rust 的迭代器模式，使得遍历非常符合直觉。
    *   `search_nodes(closure)`: 允许传入闭包进行灵活的节点查找。
*   **拓扑修改**:
    *   `add(Node)`: 添加节点，返回 ID。
    *   `add_child(Node, ParentId)`: 添加子节点。由于所有权在 Tree 内部，用户不需要手动 `malloc/free`，且 ID 引用保证了不会出现野指针。
*   **距离与比较**:
    *   内置了丰富的距离计算方法：`robinson_foulds`, `robinson_foulds_norm`, `weighted_robinson_foulds`, `kuhner_felsenstein`。这表明它更侧重于**系统发育分析**和**树的比较**，而不仅仅是操作。
    *   `get_partitions`: 获取树的二分（Bipartitions）集合，这是计算树距离的基础。
*   **安全性与错误处理**:
    *   大量使用 `Result<T, TreeError>`。例如，访问不存在的 ID 会返回 `Err(NodeNotFound)`，而不是段错误。
    *   `Option` 用于处理可选值（如名称、边长），避免了空指针异常。

---

## 5. 其他重要区别

除了核心数据结构和 API，两个项目在生态定位和底层实现上也有显著差异。

### 解析实现 (Parsing)

*   **Newick Utilities**: 采用**形式化语法**方法。
    *   使用 `Flex` (Lexer) 和 `Bison` (Parser) 进行词法和语法分析。
    *   这使得它对 Newick 格式的处理极其健壮，能够处理复杂的边缘情况，并且有严格的语法定义 (`newick_parser.y`)。
*   **Phylotree-rs**: 采用**手写状态机**方法。
    *   使用手写的字符遍历循环和状态枚举 (`Name`, `Length`, `Comment`) 实现解析。
    *   优点是零依赖，编译速度快；缺点是可能不如 Bison 生成的解析器严谨（例如源码中留有 `TODO: handle escaped quotes` 的注释），且维护复杂的语法变更更困难。

### 生态定位 (Ecosystem)

*   **Newick Utilities**: **Unix 工具集 (Suite)**。
    *   它不仅仅是一个库，更是一组功能单一、通过管道组合的命令行工具 (`nw_display`, `nw_stats`, `nw_reroot`, `nw_prune` 等)。
    *   遵循 Unix 哲学（Do one thing well），适合 Shell 脚本集成和管道处理。
    *   包含可视化功能（SVG/ASCII 绘图）。
*   **Phylotree-rs**: **分析库 (Library) + 多功能工具**。
    *   核心是一个 Rust Crate，旨在被其他 Rust 程序集成。
    *   提供了一个单一的“瑞士军刀”式 CLI 工具 (`phylotree`)，通过子命令 (`phylotree compare`, `phylotree generate` 等) 调用。
    *   功能上更侧重于**树的比较**（计算 Robinson-Foulds 距离等）和**随机树生成**，而不是图形化展示或复杂的文本处理。

### 总结建议

*   如果您需要**快速处理**、**绘图**或在 **Shell 脚本**中清洗 Newick 数据，**Newick Utilities** 是不二之选。
*   如果您需要**开发**高性能的系统发育分析软件，或者需要计算树之间的**拓扑距离**，**Phylotree-rs** 提供了更现代、安全的 Rust 接口。

---

## 6. pgr Newick 模块实施计划

本文档详细规划了在 `pgr` 项目中从零实现 Newick 格式处理模块组的路线图。该实现旨在结合 C 语言 `newick_utils` 的健壮性与 Rust `phylotree-rs` 的安全性。

### 理论基础：树比较核心概念 (Core Concepts of Tree Comparison)

在规划高级分析功能（如距离计算、拓扑比较）之前，有必要明确以下系统发育树比较的核心概念：

#### 1. 基础积木：划分 (Splits / Bipartitions)
**Splits** 是描述无根树拓扑结构的最基本单元。
*   **定义**: 树上的每一条**内部边**（Edge）都将所有叶子节点（Taxa）划分成两个互不相交的集合 $\{A, B\}$。这种划分被称为一个 Split，通常记作 $A|B$。
*   **性质**: 一棵树的拓扑结构可以**完全**由它所包含的所有 Splits 集合来定义。
*   **实现**: 在计算机中通常使用 **BitSet** 高效存储。例如叶子为 $\{A,B,C,D\}$，则 Split $\{A,B\}|\{C,D\}$ 可表示为 `1100`。

#### 2. Robinson-Foulds (RF) 距离
**RF 距离** 是基于 Splits 的最经典距离度量。
*   **计算公式**: $RF = |S_1 \setminus S_2| + |S_2 \setminus S_1|$。即两棵树中**不共享**的 Splits 总数（对称差）。
*   **特点**: **极度敏感**（"全有或全无"指标），但**计算极快**（线性复杂度 $O(n)$）。
*   **适用场景**: 快速检查树的完全一致性，或作为基础的差异度量。

#### 3. Triplet Distance (三元组距离)
**Triplet** 是**有根树 (Rooted Trees)** 的最小信息单元。
*   **定义**: 任意三个叶子 $\{x, y, z\}$ 在有根树中的拓扑关系。例如 $((x,y),z)$ 表示 $x,y$ 更亲缘。
*   **特点**: **仅限有根树**。比 RF 更**鲁棒**（High Resolution）。**tqDist** (Sand et al. 2014) 提供了针对一般多叉树 (General Trees) 的 $O(n \log n)$ 高效算法。

#### 4. Quartet Distance (四元组距离)
**Quartet** 是**无根树 (Unrooted Trees)** 的最小信息单元。
*   **定义**: 任意四个叶子 $\{a, b, c, d\}$ 在无根树中的拓扑关系。四者只有三种可能的二分拓扑：$ab|cd$，$ac|bd$，$ad|bc$。
*   **特点**: **最强鲁棒性**。目前公认衡量无根树相似度的最佳指标之一。但**计算复杂**（**tqDist** 算法实现了针对多叉树的 $O(d \cdot n \log n)$ 复杂度，解决了传统算法仅限二叉树或耗时 $O(n^2)$ 的问题）。

#### 5. Quartet Sampling (QS)
**QS** (Pease et al. 2018) 是一种现代的**分支支持度评估**方法，用于替代或补充传统的 Bootstrap。
*   **核心思想**: 传统 Bootstrap 只能给出“支持频率”，而 QS 利用 Quartet 的拓扑分布来区分**冲突 (Conflict)** 和 **信号缺失 (Lack of Signal)**。
*   **指标体系**:
    *   **QC (Quartet Concordance)**: 一致性 $[-1, 1]$。类似于 Bootstrap，越高越好。
    *   **QD (Quartet Differential)**: 偏向性 $[0, 1]$。衡量两种错误拓扑是否均匀分布。若 QD 接近 0，暗示存在特定方向的系统性冲突（如渐渗）。
    *   **QI (Quartet Informativeness)**: 信息量 $[0, 1]$。衡量有多少 Quartet 是有信息量的（非星状）。
*   **优势**: 能够深入剖析低支持率的成因（是数据太乱还是数据太少）。

#### 概念对比总结


| 概念 | 适用对象 | 核心逻辑 | 敏感度 | 计算复杂度 | 比喻 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Splits** | 无根/有根 | **边** (Edge) 的存在性 | - | - | 树的骨架 |
| **RF** | 无根/有根 | **数不同的边** | 极高 (脆弱) | 低 $O(n)$ | 严格考官：错一题扣大分 |
| **Triplets**| **有根树** | **数不同的三人组** | 中等 (鲁棒) | 中 $O(n \log n)$ | 民主投票：看三人小组意见 |
| **Quartets**| **无根树** | **数不同的四人组** | 中等 (鲁棒) | 高 $O(n \log n)$ | 民主投票：看四人小组意见 |

### 目标 (Goals)

1.  **原生支持**: 在 `pgr` 内部实现完整的 Newick 解析、存储和操作，不依赖外部二进制。
2.  **高性能**: 采用 Arena 内存模型，避免大量小对象的堆分配。
3.  **健壮性**: 能够处理带有注释 (`[&NHX...]`)、引号标签、无长度分支等复杂 Newick 变体。
4.  **功能丰富**: 提供统计、重定根、剪枝、格式化等常用功能。

### 架构设计 (Architecture)

#### 模块位置
实际实现采用了更紧凑的结构，将算法和遍历逻辑集成在 `tree.rs` 中：
```text
src/libs/phylo/
├── mod.rs          # 模块导出
├── node.rs         # 核心数据结构 (Node)
├── tree.rs         # 核心数据结构 (Tree) 及基础操作 (Traversal, Mod)
├── algo.rs         # 高级算法 (Sort, etc.)
├── parser.rs       # Newick 解析器 (基于 Nom)
└── error.rs        # 错误定义 (TreeError)
```

**设计决策说明 (Design Rationale)**:
*   **分层设计**: 核心数据结构和基础遍历逻辑保持在 `tree.rs` 中，而具体的应用层算法（如排序）分离至 `algo.rs`，保持代码结构清晰。

#### 核心数据结构 (Data Structures)
采用 **Arena (Vector-backed)** 模式，参考 `phylotree-rs` 但进行优化。

### 设计思路与对比 (Design Rationale & Comparison)

本模块的设计深受 `phylotree-rs` 启发，但在数据结构和功能侧重上做了以下权衡和改进：

1.  **内存模型 (Memory Model)**:
    *   **相同点**: 两者均采用 **Arena** 模式（即 `Vec<Node>` 存储所有节点，使用 `usize` 索引代替指针）。这是 Rust 中处理图/树结构的惯用做法，能有效避免引用循环，提高缓存局部性。
    *   **差异**: `pgr` 的 `Node` 结构更加精简。

2.  **节点设计 (Node Design)**:
    *   **`phylotree-rs`**:
        ```rust
        pub struct Node {
            // "EdgeLength" 代表进化距离（如碱基替换率或时间）
            pub parent_edge: Option<EdgeLength>,
            
            // 冗余存储：父节点同时也存了一份所有子节点的边长
            // 场景：在从上往下遍历计算时，无需访问子节点内存即可获取边长
            pub(crate) child_edges: Option<HashMap<NodeId, EdgeLength>>, 
            
            pub(crate) subtree_distances: RefCell<...>, // 内部可变性缓存
            
            // 软删除标记
            // 原因：在 Vec 存储中，真正删除元素会导致后续所有 Index 偏移，破坏树结构。
            // 做法：标记为 true，保留位置但视为不存在。
            pub(crate) deleted: bool, 
            
            // 原始注释字符串
            // 存储 Newick 中 [] 内的内容，未做解析。
            pub comment: Option<String>, 
            // ...
        }
        ```
        `phylotree-rs` 的设计包含大量**运行时缓存**（如距离矩阵、深度）和**冗余关系**（子节点边长表）。这使得它适合频繁查询的场景，但也增加了内存开销和维护状态一致性的复杂度。
    *   **`pgr::phylo`**:
        ```rust
        pub struct Node {
            // 核心设计：边长归属于子节点
            // 解释：在有根树中，每个节点（除根外）只有一条指向父节点的边。
            pub length: Option<f64>, 
            
            // 结构化属性替代原始 comment
            pub properties: Option<BTreeMap<String, String>>, 

            // 软删除标记 (响应用户建议)
            // 优点：删除节点时 O(1) 且保持 ID 稳定，避免牵一发而动全身。
            // 配合 Tree::compact() 方法在需要时进行垃圾回收。
            pub deleted: bool,
        }
        ```
        *   **去冗余**: 我们只存储 Newick 标准定义的最小信息量。计算属性（如深度、距离）将通过算法按需计算，不存储在节点中。
        *   **NHX 支持**: `phylotree-rs` 仅将注释视为字符串 (`Option<String>`)。我们将注释升级为 **`BTreeMap`**，原生支持 **NHX (New Hampshire X)** 格式的键值对。这使得 `pgr` 能更方便地处理元数据（如 `&&NHX:S=human`）。
        *   **确定性**: 使用 `BTreeMap` 而非 `HashMap`，确保在序列化输出时属性顺序固定（按键排序），保证 CLI 工具输出的**确定性 (Determinism)**，便于 diff 和测试。

3.  **解析策略 (Parsing Strategy)**:
    *   `phylotree-rs` 使用手写状态机，代码较难维护。
    *   `pgr` 引入 **`nom 8`**，利用 Parser Combinator 构建更健壮、易扩展且高性能的解析器。

```rust
// src/libs/phylo/node.rs
use std::collections::BTreeMap;

// 节点索引类型，轻量且安全
pub type NodeId = usize;

pub struct Node {
    pub id: NodeId,
    pub parent: Option<NodeId>,
    pub children: Vec<NodeId>,
    
    // Payload 数据
    pub name: Option<String>,
    pub length: Option<f64>,     // 分支长度
    pub properties: Option<BTreeMap<String, String>>, // 结构化 NHX 注释 (key=value)，按键排序
    pub deleted: bool,           // 软删除标记
}
```

```rust
// src/libs/phylo/tree.rs
use super::node::{Node, NodeId};

pub struct Tree {
    nodes: Vec<Node>,            // Arena 存储
    root: Option<NodeId>,        // 根节点 ID
}
```

#### 解析器选型 (Parsing Strategy)
已引入 **`nom`** (Parser Combinator) 库。
*   **理由**: `phylotree-rs` 的手写状态机难以维护且存在 TODO；`newick_utils` 的 Flex/Bison 难以在 Rust 中集成。`nom` 是 Rust 生态中最成熟的解析库，性能极高且易于处理嵌套结构和转义字符。
*   **实现**: 已在 `parser.rs` 中使用 `nom 8` 实现了完整的 Newick 语法解析。

#### 错误处理设计 (Error Handling Design)
为了提供良好的开发者体验和用户反馈，`phylo` 模块实现了专门的错误处理机制：
*   **TreeError 枚举**: 定义了 `ParseError` 和 `LogicError` 两类错误。
    *   `ParseError`: 携带 `line`, `column`, `snippet` 等上下文信息，当 Newick 格式错误时，能够精确指出出错位置（如 "Missing semicolon at line 1, column 15"）。
    *   `LogicError`: 处理运行时逻辑错误，如访问不存在的节点、在非连通图中查找 LCA 等。
*   **Doc Tests 覆盖**: 关键 API（如 `Tree::from_newick`, `get_path_from_root`）均包含 "Error handling" 的文档测试示例，明确展示了错误触发条件和处理方式。

### 实施步骤 (Implementation Roadmap)

#### Phase 1: 基础架构 (Infrastructure)
*   [x] 创建 `src/libs/phylo/` 目录结构。
*   [x] 实现 `Node` (`node.rs`) 和 `Tree` (`tree.rs`) 结构体。
*   [x] 实现基础方法：`add_node`, `add_child`, `get_node`, `remove_node` (soft), `compact` (gc)。

#### Phase 2: 解析与序列化 (Parsing & Serialization)
*   [x] 添加 `nom` 依赖。
*   [x] 定义 `TreeError` 枚举，提供详细的解析错误上下文（如出错位置）。
*   [x] 实现 Newick 语法定义 (BNF 转换)。
    *   [x] 处理 `Label` (支持引号和转义)。
    *   [x] 处理 `Length` (支持科学计数法)。
    *   [x] 处理 `Comment` (方括号内容)。
*   [x] 实现 `Tree::from_newick(str) -> Result<Tree>`.
*   [x] 实现序列化输出：
    *   [x] `to_newick()`: 紧凑格式。
    *   [x] `to_newick_with_format(indent)`: 支持缩进和换行的美化输出。
    *   [x] `to_dot()`: 导出 Graphviz DOT 格式用于可视化。
*   [x] 单元测试：覆盖各种 Newick 变体及序列化一致性。

#### Phase 3: 遍历与查询 (Traversal & Query) - 参考 phylotree-rs
*   [x] 实现迭代器：
    *   [x] `preorder`: 先序遍历。
    *   [x] `postorder`: 后序遍历 (适合计算，如 dp)。
    *   [x] `levelorder`: 层序遍历。
*   [x] 实现路径与距离查询：
    *   [x] `get_path_from_root(node_id)`: 获取从根到节点的路径。
    *   [x] `get_distance(node_a, node_b)`: 计算两个节点间的距离（边长总和 & 边数量）。
*   [x] 实现最近公共祖先 (LCA)：
    *   [x] `get_common_ancestor(node_a, node_b)`。
*   [x] 实现子树与查找：
    *   [x] `get_subtree(node_id)`: 获取子树所有节点。
    *   [x] `get_leaves()`: 获取所有叶子节点。
    *   [x] `find_nodes(predicate)`: 根据条件查找节点。
    *   [x] `get_node_by_name(name)`: 根据名称查找节点。

#### Phase 4: 高级算法 (Advanced Algorithms)
*   [x] **Reroot**: 实现 `reroot_at(node_id)`，涉及父子关系翻转和边长重新分配。
*   [x] **Prune**: 剪掉指定名称或正则匹配的节点。

#### Phase 5: CLI 集成 (Integration) - 模仿 Newick Utilities
我们的 CLI 功能将模仿 `newick_utils` 工具集。以下是目标工具列表及其功能映射：

*   **`nw_stats`** $\to$ **`pgr nwk stat`**
    *   **功能**: 树的统计信息 (节点数, 深度, 类型等)
    *   **状态**: **[x] 已实现** (支持多树处理, TSV/KV 输出, 统计二叉分枝)

*   **`nw_distance`** $\to$ **`pgr nwk distance`**
    *   **功能**: 计算节点间距离 / 树间距离
    *   **状态**: **[x] 已实现** (支持 root, parent, pairwise, lca, phylip)

*   **`nw_indent`** $\to$ **`pgr nwk indent`**
    *   **功能**: 格式化/缩进 Newick 树
    *   **状态**: **[x] 已实现** (支持紧凑/缩进输出)

*   **`nw_display`** $\to$ **`pgr nwk to-dot` / `to-forest` / `to-tex`**
    *   **功能**: 树的可视化 (Graphviz/LaTeX Forest)
    *   **状态**: **[x] 已实现** (支持 Graphviz DOT, LaTeX Forest 代码及完整文档导出)

*   **`nw_topology`** $\to$ **`pgr nwk topo`**
    *   **功能**: 仅保留拓扑结构 (去除分支长度)
    *   **状态**: **[x] 已实现** (支持 `--bl`, `--comment`, `-I`, `-L`)

*   **`nw_labels`** $\to$ **`pgr nwk label`**
    *   **功能**: 提取所有标签 (叶子/内部节点)
    *   **状态**: **[x] 已实现** (支持正则过滤, 内部/叶子筛选, 单行输出)

*   **`nw_reroot`** $\to$ **`pgr nwk reroot`**
    *   **功能**: 重定根 (Outgroup, Midpoint)
    *   **状态**: **[x] 已实现** (支持 Midpoint, Outgroup (LCA), Lax mode, Deroot)

*   **`nw_prune`** $\to$ **`pgr nwk prune`**
    *   **功能**: 剪枝 (移除指定节点)
    *   **状态**: **[x] 已实现** (支持正则/列表，自动清理，反选)

*   **`nw_clade`** $\to$ **`pgr nwk subtree`**
    *   **功能**: 提取子树 (Clade)
    *   **状态**: **[x] 已实现** (支持 context 扩展、单系群检查、正则匹配)

*   **`nw_order`** $\to$ **`pgr nwk order`**
    *   **功能**: 节点排序 (Ladderize)
    *   **状态**: **[x] 已实现** (支持 alphanumeric/descendants/list/deladderize)

*   **`nw_rename`** $\to$ **`pgr nwk rename` / `replace`**
    *   **功能**: 重命名节点 (Map file/Rule)
    *   **状态**: **[x] 已实现** (Split into `rename` & `replace`)

*   **`nw_condense`** $\to$ **`pgr nwk condense`**
    *   **功能**: 压缩树 (合并短枝/多叉化)
    *   **状态**: [ ] (`subtree` 命令支持 `-C` / `--condense` 选项来压缩子树)

*   **`nw_support`** $\to$ **`pgr nwk support`**
    *   **功能**: 计算/显示支持率 (Bootstrap)
    *   **状态**: [ ] (`reroot` 命令支持 `-s` / `--support-as-labels` 处理支持率)

*   **`nw_match`** $\to$ **`pgr nwk match`**
    *   **功能**: 匹配两棵树的节点
    *   **状态**: [ ]

*   **`nw_ed`** $\to$ **`pgr nwk ed`**
    *   **功能**: 编辑距离 / 树操作脚本
    *   **状态**: [ ]

*   **`nw_gen`** $\to$ **`pgr nwk gen`**
    *   **功能**: 生成随机树
    *   **状态**: [ ]

*   **`nw_duration`** $\to$ **`pgr nwk duration`**
    *   **功能**: (通常指时间树相关)
    *   **状态**: [ ]

---

## 7. 测试与验证 (Testing & Verification)

为了确保 `pgr nwk` 模块的正确性和鲁棒性，特别是作为 `newick_utils` 的替代品，我们采取了多层次的测试策略。

### 测试策略

1.  **单元测试 (Unit Tests)**:
    *   针对 `libs/phylo` 中的核心逻辑（解析、遍历、算法）。
    *   覆盖各种 Newick 格式变体（引号、注释、无长度、科学计数法等）。
    *   位于 `src/libs/phylo/` 源码文件中。

2.  **集成测试 (Integration Tests)**:
    *   针对 CLI 子命令的端到端测试。
    *   使用 `assert_cmd` 模拟命令行调用，验证标准输出（stdout）和退出代码。
    *   位于 `tests/` 目录下的 `cli_nwk_*.rs` 文件 (如 `cli_nwk_subtree.rs`, `cli_nwk_ops.rs` 等)。

3.  **兼容性测试 (Compatibility Tests)**:
    *   直接参考 `newick_utils` 的测试用例 (`tests/test_nw_*_args`)。
    *   确保 `pgr` 的输出与 `newick_utils` 在关键场景下一致（或语义一致）。

### 兼容性测试状态矩阵

以下列出了参考 `newick_utils` 原生测试套件实现的兼容性测试进度：

#### `nw_stats` -> `pgr nwk stat`

| Test Case (newick_utils) | 描述 | pgr 对应测试 | 状态 | 备注 |
| :--- | :--- | :--- | :--- | :--- |
| `def` | 默认输出 (Key-Value) | `command_stat` | **[x] Pass** | 验证节点数、叶子数、二叉分枝数等 |
| `fl` | Flat Line 输出 (TSV) | `command_stat_style_line` | **[x] Pass** | 对应 `--style line` 参数 |
| `many` | 多树处理 (Multi-tree) | `command_stat_multi_tree` | **[x] Pass** | 验证单文件多棵树的连续处理 |

#### `nw_clade` -> `pgr nwk subtree`

| Test Case (newick_utils) | 描述 | pgr 对应测试 | 状态 | 备注 |
| :--- | :--- | :--- | :--- | :--- |
| `def` | 默认提取 | `command_subtree_default` | **[x] Pass** | 基础子树提取 |
| `regex` | 正则匹配提取 | `command_subtree_regex` | **[x] Pass** | 支持 `-r` 正则匹配 |
| `context` | 上下文扩展 | `command_subtree_context` | **[x] Pass** | 支持 `-c` 扩展 N 层 |
| `monop` | 单系群检查 | `command_subtree_monophyly` | **[x] Pass** | 支持 `-M` 单系性验证 |

| `nw_indent` | 缩进/格式化 Newick 字符串 | `pgr nwk indent` / `fmt` | **[x] Pass** | |
| `nw_trim` | 修剪树 (Trim) | `pgr nwk trim` | [ ] | (部分功能被 `prune` 覆盖) |

*注：我们将实现其中大部分功能，具体命令名称可能会根据 `pgr` 的整体风格进行微调（例如使用动词作为子命令）。*

### 依赖管理 (Dependencies)
已在 `Cargo.toml` 中添加：
```toml
[dependencies]
nom = "8"  # 用于高性能解析
```

### 测试计划 (Testing)
1.  **Unit Tests**: 针对 parser 的每个组件编写测试 (已完成)。
2.  **Doc Tests**: 为关键 API 提供文档测试，覆盖正常用法和错误处理 (已完成)。
3.  **Integration Tests**: 使用 `newick_utils` 生成的标准文件作为输入，验证解析结果的一致性。
4.  **Property Tests**: (可选) 生成随机树并验证 `parse(to_newick(tree)) == tree`。

---

## 8. API 对比与差距分析 (pgr vs phylotree-rs)

以下列表对比了 `phylotree-rs` (参考版本) 提供的公开 API 与 `pgr::phylo` 当前实现的覆盖情况。

### 已实现 (Implemented)

这些核心功能已经移植或重构，能够满足基础操作需求。

*   **Tree Structure**:
    *   `Tree::new()`: 创建空树。
    *   `Tree::add_node()`, `Tree::add_child()`: 构建树结构。
    *   `Tree::get_node()`, `Tree::get_root()`: 访问节点。
    *   `Tree::len()`: 节点总数 (对应 `size()`)。
*   **Parsing**:
    *   `Tree::from_newick()`: 解析 Newick 字符串 (支持引号、注释、科学计数法)。
*   **Serialization**:
    *   `to_newick()`: 紧凑格式输出。
        *   *注*: `phylotree-rs` 返回 `Result<String>`, `pgr` 返回 `String` (Infallible)。
    *   `to_newick_with_format()`: 支持缩进的格式化输出。
    *   `to_dot()` (Graphviz): 输出 DOT 格式，可用于可视化 (**pgr 特有**)。
*   **Traversal**:
    *   `preorder`, `postorder`: 深度优先遍历 (迭代器风格)。
    *   `levelorder`: 广度优先遍历。
*   **Query**:
    *   `get_leaves()`: 获取所有叶子节点。
    *   `get_path_from_root()`: 获取根到节点的路径。
    *   `get_common_ancestor()` (LCA): 最近公共祖先。
    *   `get_distance()`: 计算节点间距离 (加权/拓扑)。
    *   `get_subtree()`: 获取子树节点集合。
    *   `find_nodes()`, `get_node_by_name()`: 查找节点。
    *   `get_height()`: 计算节点高度 (到最远叶子的距离)。
    *   `is_monophyletic()`: 判断是否为单系群 (**pgr 特有**)。
*   **Modification**:
    *   `reroot_at()`: 重新定根 (支持边长重分配)。
    *   `prune_where()`: 剪枝 (删除匹配节点及其子孙)。
    *   `remove_node()`: 软删除单个节点。
    *   `collapse_node()`: 压缩节点 (合并边长)。
    *   `compact()`: 物理删除软删除节点并重构树 (**pgr 特有**)。

### 统计与计算 (Statistics & Calculation)

*   `is_binary()`: 检查是否为二叉树。
*   `get_leaf_names()`: 获取所有叶子节点的名称列表。
*   `get_partitions()`: 获取树的二分 (Bipartitions) 集合 (返回 `HashSet<BTreeSet<String>>`)，是计算 RF 距离的基础。
*   `diameter()`: 树的直径 (最远叶子间距离)。
*   `robinson_foulds()`: 计算两棵树的 Robinson-Foulds 距离 (拓扑差异)。

### 未实现 (Missing / Differences)

以下功能在 `phylotree-rs` 中存在，但在 `pgr` 中尚未实现或有不同设计。

*   **Tree Statistics (统计指标)**:
    *   `cherries()`: 计算 Cherry 数量。
    *   `colless()`, `colless_yule()`, `colless_pda()`: Colless 平衡指数。
    *   `sackin()`, `sackin_yule()`, `sackin_pda()`: Sackin 平衡指数。
*   **Tree Comparison (树比较)**:
    *   `compare_topologies()`: 拓扑结构比较。
*   **Internal Caching**:
    *   `get_node_depth()`: 缓存节点深度 (目前需实时计算)。
*   **Traversal**:
    *   `inorder`: 中序遍历 (仅适用于二叉树，`pgr` 支持多叉树故未直接实现)。

### 计划中 (Planned)

*   `robinson_foulds()`:
    *   `weighted_robinson_foulds()`: 加权 RF 距离。
*   **Visualization**:
    *   **Graphviz DOT**: `to_dot()` 已实现。
    *   **LaTeX Forest**: `to_forest` (Raw code) 和 `to_tex` (Full document) 已实现。
    *   `print_entity()` (或类似): 在终端打印 ASCII 树状图，用于快速调试和展示。
*   **Tree Generation**:
    *   `generate_random_tree()` (Yule/Coalescent 模型): 主要用于模拟研究。`pgr` 侧重于处理真实数据，除非用于测试生成，否则优先级较低。
*   **Complex I/O**:
    *   `from_file()`: `pgr` 通常通过 CLI 处理文件读取，核心库只需处理字符串或 Buffer。
*   **Internal Caching**:
    *   `update_depths()`, `matrix`: `phylotree-rs` 缓存了大量中间状态。`pgr` 倾向于按需计算 (On-demand) 以保持轻量化。

---

## 7. Visualization Details (LaTeX Forest)

`pgr` 的可视化功能深度集成了 LaTeX Forest 包，配合精心设计的模板 (`docs/template.tex`)，能够生成出版级质量的进化树。

### 7.1 核心命令

*   **`pgr nwk to-forest`**: 生成原始 Forest 代码。适合嵌入现有 LaTeX 文档。
*   **`pgr nwk to-tex`**: 生成完整 `.tex` 文档。自动合并模板，可直接用 `xelatex` 编译。

### 7.2 样式系统 (Styles)

模板定义了四种核心样式，可以通过 Newick 文件中的 NHX 注释直接调用：

1.  **`dot` (节点圆点)**
    *   **效果**: 在节点处绘制实心圆点。
    *   **用法**: `[&&NHX:dot=red]` (指定颜色) 或自动应用于带名称的内部节点。

2.  **`bar` (垂直短杠)**
    *   **效果**: 在父节点与子节点的连线上绘制垂直短杠，常用于标记性状演化或事件。
    *   **用法**: `[&&NHX:bar=blue]`。

3.  **`rec` (背景矩形)**
    *   **效果**: 为整个子树（Clade）绘制背景矩形框。利用 `fit to=tree` 实现。
    *   **用法**: `[&&NHX:rec=LemonChiffon]`。常配合模板中定义的柔和色系使用。

4.  **`tri` (三角形)**
    *   **效果**: 在节点右侧绘制三角形，常用于表示折叠的子树 (Collapsed Clade) 或强调叶节点。
    *   **用法**: `[&&NHX:tri=green]`。

### 7.3 颜色与全局设置

*   **配色方案**: 模板内置了一组柔和的莫兰迪色系（如 `ChampagnePink`, `TeaRose`, `Celadon` 等）。
*   **自动对齐**: 默认启用 `tier=word`，强制所有叶节点对齐（Cladogram 风格）。
*   **字体支持**:
    *   **默认**: 使用 `Noto Sans` 系列（需安装），兼容性好。
    *   **高级 (`--style`)**: 保留模板中预设的 `Fira Sans` (英) 和 `Source Han Sans SC` (中) 设置，适合需要特定设计感的场景。

### 7.4 高级特性 (Advanced Features)

*   **Phylogram 模式 (`--bl`)**:
    *   绘制带分支长度的系统发育树。
    *   **自动比例尺**: 程序会根据树高自动计算合适的比例尺（如 0.01, 0.05, 1.0 等），并绘制在右下角。
*   **Forest 直通车 (`--forest`)**:
    *   允许将外部生成的 Forest 代码文件（非 Newick）直接嵌入模板生成 PDF。
*   **特殊字符处理**:
    *   Newick 名称中的下划线 `_` 会被自动转换为空格，避免 LaTeX 编译错误。

### 7.5 工作流示例

1.  **准备数据**: 使用 `pgr nwk comment` 命令或是手动为节点添加样式注释。
    ```bash
    # 为节点 A 和 B 的最近公共祖先 (LCA) 添加背景矩形和标签
    pgr nwk comment input.nwk --lca A,B --rec TeaRose --label Group1 > annotated.nwk
    ```
2.  **转换**:
    ```bash
    pgr nwk to-tex annotated.nwk > output.tex
    ```
3.  **编译**:
    ```bash
    tectonic output.tex
    ```

---

## 9. 使用示例 (Usage Examples)

### pgr nwk stat

统计 Newick 文件的基本信息（节点数、叶子数、二分歧节点数等）。

```bash
# 默认输出 (Key-Value 格式)
pgr nwk stat data.nwk

# 表格输出 (TSV 格式，适合后续处理)
pgr nwk stat data.nwk --style line
```

### pgr nwk indent

格式化 Newick 树，使其更易读，或压缩为单行。

```bash
# 默认缩进 (2个空格)
pgr nwk indent data.nwk

# 自定义缩进字符 (例如使用4个空格)
pgr nwk indent data.nwk --text "    "

# 压缩为单行 (Compact)
pgr nwk indent data.nwk --compact
```

```text
整合 indent.rs，测试在 `tests/cli_nwk_ops.rs` 中。
你看看 nw_indent 的代码，有没有什么值得参考的
把 tests/ 里 nw_stats 的相关测试迁移到本项目，需要的测试材料 可以拷贝到 `tests/newick/` 目录下
改进帮助文本。 参考 nw_indent 的帮助文本，按我们自己的样式 进行调整。
phylo.md 根据现在的代码状态，更新文档
```

---

## 10. 附录：Workflow 参考 (Appendix: Workflow Reference)

### Bootscan Workflow (`bootscan.sh`)

`newick_utils` 源码中的 `bootscan.sh` 展示了一个结合多种生物信息学工具进行重组检测的完整流程。这为 `pgr` 的 CLI 设计提供了实际应用场景参考。

**流程步骤详解：**

1.  **序列比对 (Alignment)**
    *   **工具**: `mafft`
    *   **操作**: 对输入的未比对序列文件（FASTA）进行多序列比对。
    *   **命令**: `mafft --quiet "$INFILE" > "$MUSCLE_OUT"`

2.  **切片 (Slicing)**
    *   **工具**: `infoalign`, `seqret` (EMBOSS 工具集)
    *   **操作**:
        *   获取比对长度 (`infoalign`).
        *   按指定步长 (`SLICE_STEP`) 和窗口大小 (`SLICE_WIDTH`) 遍历比对结果。
        *   将每个窗口切片并转换为 PHYLIP 格式 (`seqret`).
    *   **目的**: 准备滑动窗口数据以构建局部树。

3.  **构建系统发育树 (Tree Building)**
    *   **工具**: `phyml`
    *   **操作**: 对每个 PHYLIP 切片文件构建最大似然树 (Maximum Likelihood Tree)。
    *   **参数**: `-b 0` (无 bootstrap，求速度), `-o n` (不优化拓扑/分支/率参数)。
    *   **输出**: 生成一系列无根树文件。

4.  **重定根 (Rerooting)**
    *   **工具**: `nw_reroot` (Newick Utilities)
    *   **操作**: 使用指定的外群 (`OUTGROUP`) 对每棵树进行定根。
    *   **命令**: `nw_reroot $unrooted_tree $OUTGROUP > ${unrooted_tree/.txt/.rr.nw}`
    *   **对应 pgr**: `pgr nwk reroot`

5.  **距离提取 (Distance Extraction)**
    *   **工具**: `nw_distance`, `nw_clade`, `nw_labels`
    *   **操作**:
        *   计算参考序列 (`REFERENCE`) 到树中其他所有节点的距离 (`nw_distance -n`).
        *   提取相关标签 (`nw_clade`, `nw_labels`).
        *   生成包含 `(Position, Distance1, Distance2, ...)` 的表格数据。
    *   **对应 pgr**: `pgr nwk distance`, `pgr nwk subtree`, `pgr nwk label`

6.  **可视化 (Visualization)**
    *   **工具**: `gnuplot`
    *   **操作**: 将距离表格绘制成折线图。横轴为序列位置，纵轴为参考序列到其他序列的遗传距离。
    *   **原理**: 如果参考序列在某个区域与其他序列的距离显著变化（例如最近邻改变），提示可能发生了重组事件。
