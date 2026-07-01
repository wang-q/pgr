# ms2dna 迁移到 Rust 设计与规划

## 目标与范围
- 目标：将 ms2dna（将 Hudson ms 输出的 0/1 单倍型转换为 DNA 序列）完整迁移到 Rust，保持 CLI 与输出格式兼容
- 范围：参数解析、随机数与种子管理、输入解析（ms 输出）、GC 设定、DNA 序列生成、FASTA 输出、多群体标签与批次命名

## 源码结构参考
- 主程序：读取参数、设定种子、迭代输入样本
- 参数/界面：getArgs、printUsage、printSplash
- 样本处理：读取 ms 输出、生成祖先序列与突变序列、FASTA 打印
- 结构体定义：样本结构与字段
- 随机数/工具：ran.c、stringUtil.c、sequence_data.c、eprintf.c

## 输入与输出
- 输入：ms 输出文本流或文件，格式包含
  - 首行命令行（解析 nsam, howmany, -r nsite, -I 群体与样本数）
  - 每个样本块的 segregating sites 与 positions 与 haplotypes
- 输出：FASTA 序列，带批次/群体/个体标签，单行序列（不换行）

## 行为概述
- 祖先序列生成：按 GC 比例（-g）随机生成祖先序列（A/T 与 G/C 各 0.5）
- 位点映射：将 segsites 的 positions 映射到 [0, nsite) 的整数位置，避免重复命中
- 变异生成：在映射位置按 GC/AT 决定派生碱基，遵循互补/随机二选一规则
- 标签规则：Lx（批次）/ Px（群体）/ Sx（样本序号）
- 位置精度增强：将 ms 的 4 位 positions 通过小扰动提升精度

## Rust 架构设计
- 模块划分
  - params：Args 的 Rust 映射与 clap 参数解析（-g、-s、-h、-v、文件列表）
  - rng：SeedManager 与 RNG 提供者（兼容 ms2dna 的 seed 行为）；提供 uniform_f64/choice
  - parser：解析 ms 输出流（首行命令行、样本块、positions、haplotypes）
  - model：Sample 结构体（nsam/howmany/nsite/maxlen/segsites/npop/sampleSizes/map/positions/haplotypes/seq/line）
  - generator：祖先序列与变异生成、positions 精度扰动、映射与避免重复命中
  - output：FASTA 打印与标签规则（单行输出）
  - cli：pgr 子命令 ms2dna 或独立可执行（建议子命令：pgr ms2dna）
- 数据结构映射
  - struct Args { h: bool, v: bool, s: u64, g: f64, files: Vec<PathBuf> }
  - struct Sample { howmany, nsam, nsite, maxlen, segsites, npop, sample_sizes, map, positions, haplotypes, seq, line }
  - 解析状态与流读取使用 BufRead，避免一次性读入大文件

## 随机数与复现性
- 采用轻量 SimpleRng；支持固定种子与系统种子
- 种子策略：命令行 -s 优先，否则使用系统时间与进程号组合；在每个输入流中一次性初始化并复用 RNG
- 固定种子下端到端输出可复现（允许微小浮点差异）

## CLI 设计
- 用法：pgr ms2dna [-g GC] [-s seed] [-o outfile] [--no-perturb] [-v] [inputFiles]
- 选项：
  - -g/--gc：GC 含量（0..1，默认 0.5）
  - -s/--seed：随机数种子（u64）；默认使用系统时间与PID
  - -o/--outfile：输出文件名；默认 stdout
  - --no-perturb：禁用位点位置微扰（保留原始 ms 位置精度）
  - -v/--verbose：打印运行路径、输入文件列表与使用的种子
  - -h/--help：打印帮助信息
- 行为：
  - 无文件时读取 stdin；多文件依次处理
  - 输出单行 FASTA 格式（不换行）


## 测试与验证
- 单元测试：
  - positions 解析与精度扰动、映射去重
  - 祖先序列生成与 GC 比例校验（统计）
  - 变异生成规则（AT/GC 派生）与边界条件
- FASTA 标签与单行输出
- 端到端：
  - 使用仓库 data/test.ms、simpleTest.dat 比对输出
  - 固定种子下比对输出（容差微小）
- 基准：
  - 大样本与大量 segsites 的解析与生成时间/内存曲线

## 迁移阶段计划
- 阶段 1：参数解析与最小骨架（stdin→FASTA 空输出占位），建立模块与结构体
- 阶段 2：解析首行命令行与样本块，构建 Sample
- 阶段 3：positions 精度提升与映射、祖先序列生成、变异写入
- 阶段 4：FASTA 输出与标签规则（单行输出）
- 阶段 5：多文件输入、错误处理与帮助/版本输出
- 阶段 6：端到端与性能优化

## 风险与规避
- 输入多样性：ms 输出格式变体（macs 提示）；对 segsites>nsite 情况输出警告并容错
- 浮点扰动：positions 微扰需控制边界与重复映射概率
- GC 比例：统计波动与随机性需通过批量测试验证
- 大文件：采用流式解析与按需分配，避免 OOM
  - macs 兼容：Rust 版仅提示“确保 positions/nsite 兼容”，不提供 -a 选项

## 集成方案
- 优先实现为 pgr 的子命令：pgr ms2dna
- 独立二进制备选：ms2dna（与 C 版名字一致），后续视需求可提供

## 参考文件
- 主程序：ms2dna.c
- 参数接口：interface.c
- 样本处理：sample.c
- 结构体定义：sample.h
