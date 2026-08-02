# pgr plot dot：PAF 共线性图（dot plot）设计

> 设计稿，日期：2026-08-03。参考 FastGA 套件的 ALNplot（`FASTGA-main/ALNplot.c`，
> 手写 EPS 的静态共线性图）；FastGA 分析见 [[../references/fastga.md]]。

## 1. 动机与选型

pgr 已有 `pgr plot hh / venn / nrps`（LaTeX/TikZ 输出），但缺"直接看比对"的可视化。
本命令把 PAF 比对画成双轴共线性图（dot plot 风格），填补这一空白。

**选型：SVG，不用 TikZ。**

*   dot plot 是数据密集图（每条比对一条线段，全基因组可达数十万条），TikZ 的
    pdflatex 编译在十万级元素下会非常慢/吃内存（ALNplot 因此默认限 10 万条）；
    SVG 是流式文本，生成即字符串拼接，浏览器直接渲染。
*   与现有 TikZ 命令不冲突：插图类（venn/hh/nrps）走 TikZ，数据密集类（dot）走 SVG。
*   SVG 无外部依赖（手写 XML，与 ALNplot 手写 EPS 同思路），发表时可用
    rsvg-convert / inkscape / cairosvg 转 PDF/PNG。

## 2. 命令接口

```text
pgr plot dot [OPTIONS] <infile>
```

*   `infile`：PAF 文件（支持 stdin / `.paf.gz`）
*   `-o`：SVG 输出（默认 stdout）
*   `--min-len`（默认 100）：最小比对长度（block_length）
*   `--min-identity`（默认 0.7）：最小 identity（matches / block_length）
*   `--max-align`（默认 100000，0=全部）：最多绘制条数（按长度取 top）
*   `--width`（默认 2000）：SVG 宽度（像素），高度按两轴总长比例自动计算

## 3. 布局与绘制

### 3.1 轴布局

*   x 轴 = target 侧，y 轴 = query 侧（PAF 列 6-9 / 列 2-4）；
*   contig 按 PAF 中首次出现的顺序排列，各自长度取该名字记录中的长度；
*   偏移 = 前面所有 contig 的累计长度（contig 之间无 gap）；
*   同一条记录中同一名字的长度一致，不一致时取最大值。

### 3.2 坐标与颜色

*   scale = width / total_target；height = scale × total_query（两轴天然同比例）；
*   线段 `(t_start,q_start) → (t_end,q_end)`，起点终点各加边距；
*   颜色 = identity 的蓝色阶（对齐 wgatools dotplot 的 blues 色阶）：
    `--min-identity`（最浅）→ `--max-identity`（最深），低于 `--min-identity`
    不画；右下角有色阶图例（wgatools 的 `labelFontSize: 20` 思路：字号随图宽缩放）；
*   刻度与网格：两轴带 bp 刻度（nice step 1/2/5×10^n，目标间隔 ~120px），
    显示真实基因组坐标（zoom 时含偏移）；淡灰网格线 + contig 分隔线；
*   轴标签：contig 名（像素宽度不足 3×字号时省略）；
*   所有线宽与字号由 `--width` 推导（线宽 = width/300，字号 = width/60，
    tick 字号 = 0.8×标签字号），任意尺寸下整体等比。
*   默认两轴同比例（1 bp = 相同像素），框高随 query/target 长度比；
    `--square` 使方框区强制正方形（两轴独立缩放，共线角失真），整体 SVG
    仍含文字边距（不等宽高）。
*   x 轴名字水平（段中心下方）、y 轴名字旋转 -90°（与序列对齐）；相邻段
    用**黑色/深灰色交替**区分；y 轴旋转名字垂直重叠时左右分列（贪心按列
    分配），左边距按列数动态扩展，顶部边距按顶部短段名字长度扩展；图例
    位于 plot 上方右侧（半透明背景，不遮任何比对线段）。
*   刻度每段从 0 开始（染色体/质粒各自计长，不累加偏移）；`--range` 放大时
    单段显示真实基因组坐标。

### 3.5 美学细节

*   y 轴翻转：query 坐标 0 在底部，正链比对从左下到右上（传统 dot plot 方向）；
*   线段 opacity 0.5，密集重叠区自然加深，避免糊成一片；
*   左边距 200px（容纳 contig 名 + 刻度），底边距 110px（刻度 + contig 名 +
    图例），图例在右下角。

### 3.4 局部放大（--range）

*   单个 `--range chr:start-end`（1-based，intspan::Range 解析），作用于 target 侧；
    query 轴**自动对准**，用户无需指定第二个区间（比对是一对一映射）。
*   只保留与区间相交的比对并把 target 坐标裁剪/相对化，轴的该 chr 段长度 = 区间长；
*   query 自动聚焦到**显著匹配簇**：保留记录按 query 位置贪心合并（gap ≤ 100kb
    并入），保留覆盖 bp ≥ 最大簇 1% 的所有簇（**同一染色体远端匹配、跨染色体
    匹配都可见**，只丢弃微小噪声碎片）；每个簇改名 `chr#k` 作为独立轴段拼接，
    刻度显示各簇真实基因组坐标，标签显示原染色体名；簇内坐标保持绝对值，
    渲染时减簇起点（负偏移由 SVG clipPath 视觉裁剪），线段方向不失真。

### 3.3 过滤

*   长度 = `block_length`，identity = `matches / block_length`；
*   超出 `--max-align` 时按长度降序取前 N（记录顺序稳定，不排序输出）。

## 4. 输出示例

```xml
<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="..." height="...">
  <g id="red" stroke="#FF0000"><line .../></g>
  <g id="blue" stroke="#0080FF"><line .../></g>
</svg>
```

## 5. 测试计划

*   单元测试（`libs/plot/dot.rs`）：小 PAF → 线段条数/颜色/坐标、过滤（长度/identity/
    max-align）生效、空输入报错；
*   集成测试（`tests/cli_plot.rs`）：`tests/plot/dot.paf` → SVG 头 + 红/蓝分组存在。
