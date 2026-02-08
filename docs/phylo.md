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
├── tree.rs         # 核心数据结构 (Tree) 及所有算法 (Traversal, Algo)
├── parser.rs       # Newick 解析器 (基于 Nom)
└── error.rs        # 错误定义 (TreeError)
```

**设计决策说明 (Design Rationale)**:
*   **整合原因**: 将遍历 (`iter`) 和算法 (`algo`) 整合进 `tree.rs` 是为了保持**高内聚**。这些操作紧密依赖 `Tree` 的内部实现（如 Arena 索引），集中管理可以简化可见性控制，符合 Rust 将方法定义在类型附近的惯例。当前 `tree.rs` 约 800 行，体积适中，无需过早拆分。

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

| Newick Utilities | 功能描述 | pgr 对应子命令 (暂定) | 状态 |
| :--- | :--- | :--- | :--- |
| `nw_stats` | 树的统计信息 (节点数, 深度, 类型等) | `pgr nwk stat` | **[x] 已实现** (支持多树处理, TSV/KV 输出, 统计二叉分枝) |
| `nw_display` | 树的可视化 (ASCII/SVG/Map) | `pgr nwk display` / `view` | [ ] |
| `nw_topology` | 仅保留拓扑结构 (去除分支长度) | `pgr nwk topology` | [ ] |
| `nw_labels` | 提取所有标签 (叶子/内部节点) | `pgr nwk labels` | [ ] |
| `nw_reroot` | 重定根 (Outgroup, Midpoint) | `pgr nwk reroot` | [ ] |
| `nw_prune` | 剪枝 (移除指定节点) | `pgr nwk prune` | [ ] |
| `nw_clade` | 提取子树 (Clade) | `pgr nwk clade` / `subtree` | [ ] |
| `nw_order` | 节点排序 (Ladderize) | `pgr nwk order` | [ ] |
| `nw_rename` | 重命名节点 (Map file/Rule) | `pgr nwk rename` | [ ] |
| `nw_condense` | 压缩树 (合并短枝/多叉化) | `pgr nwk condense` | [ ] |
| `nw_distance` | 计算节点间距离 / 树间距离 | `pgr nwk dist` | [ ] |
| `nw_support` | 计算/显示支持率 (Bootstrap) | `pgr nwk support` | [ ] |
| `nw_match` | 匹配两棵树的节点 | `pgr nwk match` | [ ] |
| `nw_ed` | 编辑距离 / 树操作脚本 | `pgr nwk ed` | [ ] |
| `nw_gen` | 生成随机树 | `pgr nwk gen` | [ ] |
| `nw_duration` | (通常指时间树相关) | `pgr nwk duration` | [ ] |
| `nw_indent` | 缩进/格式化 Newick 字符串 | `pgr nwk indent` / `fmt` | [ ] |
| `nw_trim` | 修剪树 (Trim) | `pgr nwk trim` | [ ] |

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

## 7. API 对比与差距分析 (pgr vs phylotree-rs)

以下列表对比了 `phylotree-rs` (参考版本) 提供的公开 API 与 `pgr::phylo` 当前实现的覆盖情况。

### 已实现 (Implemented)

这些核心功能已经移植或重构，能够满足基础操作需求。

*   **Tree Structure**:
    *   `Tree::new()`: 创建空树。
    *   `Tree::add_node()`, `Tree::add_child()`: 构建树结构。
    *   `Tree::get_node()`, `Tree::get_root()`: 访问节点。
*   **Parsing**:
    *   `Tree::from_newick()`: 解析 Newick 字符串 (支持引号、注释、科学计数法)。
*   **Serialization**:
    *   `to_newick()`: 紧凑格式输出。
        *   *注*: `phylotree-rs` 返回 `Result<String>`, `pgr` 返回 `String` (Infallible)。
        *   *注*: 实现代码位于 `src/libs/phylo/writer.rs`，但通过 `Tree` 方法暴露以保持兼容。
    *   `to_newick_with_format()`: 支持缩进的格式化输出。
    *   `to_dot()` (Graphviz): 输出 DOT 格式，可用于可视化。
        *   *注*: 这是 `pgr` 特有的功能，`phylotree-rs` 未直接提供。
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
*   **Modification**:
    *   `reroot_at()`: 重新定根 (支持边长重分配)。
    *   `prune_where()`: 剪枝 (删除匹配节点及其子孙)。
    *   `remove_node()`: 软删除单个节点。
    *   `compact()`: 物理删除软删除节点并重构树。

### 未实现但需要 (Missing & Planned)

*   **Internal Caching**:
    *   `get_node_depth()`: 缓存节点深度。

### 不太需要/低优先级 (Low Priority / Not Needed)

这些 API 要么用途有限，要么与 `pgr` 的设计哲学不符，或者根据当前需求被认为不重要。

*   **Comparison**:
    *   `robinson_foulds_distance()`: 计算两棵树的拓扑差异 (RF Distance)。
    *   `weighted_robinson_foulds()`: 加权 RF 距离。
    *   `get_partitions()`: 获取树的二分 (Bipartitions) 集合。
*   **Visualization**:
    *   `print_entity()` (或类似): 在终端打印 ASCII 树状图，用于快速调试和展示。
*   **Tree Generation**:
    *   `generate_random_tree()` (Yule/Coalescent 模型): 主要用于模拟研究。`pgr` 侧重于处理真实数据，除非用于测试生成，否则优先级较低。
*   **Complex I/O**:
    *   `from_file()`: `pgr` 通常通过 CLI 处理文件读取，核心库只需处理字符串或 Buffer。
*   **Internal Caching**:
    *   `update_depths()`, `matrix`: `phylotree-rs` 缓存了大量中间状态。`pgr` 倾向于按需计算 (On-demand) 以保持轻量化。

