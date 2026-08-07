# RepeatMasker 安装与自定义库（TnCentral）使用记录

## 结论

RepeatMasker 本身**不绑定任何数据库**。开箱即可用自定义 FASTA 库：

```bash
RepeatMasker genome.fa -lib my_library.fa
```

TnCentral 的 FASTA（`~/data/repeats/tncentral.fa.gz`）完全走这条路，**不需要也不能"装进"程序**。
只有按物种自动选库的 `-species` 模式才需要 FamDB（Dfam/RepBase 的 H5 格式），TnCentral 不是这种格式。

## 本地安装状态（2026-08-07）

- 源码：`/home/wangq/Scripts/pgr/RepeatMasker/`，版本 4.2.4（官方最新）。
- 搜索引擎：RMBlast 2.14.1（CBP 编译版，`~/.cbp/bin`，要求 glibc ≤2.16，
  CentOS 7 可跑），经合并目录
  `/home/wangq/Scripts/pgr/RepeatMasker/rmblast-cbp-bin` 配置为默认引擎。
- TRF：`/home/wangq/.cbp/bin/trf`（系统已有）。
- FamDB：未配置（因此只有 `-lib` 模式可用；`-species` 不可用）。
- 已通过 `perl ./configure` 完成配置（RMBLAST_DIR 指向 `rmblast-cbp-bin`），
  `./RepeatMasker -h` 可运行。RMBlast 的 tar 包仍在 `~/Downloads/`，/tmp 里
  的解压副本可删。

> **为什么是合并目录**：configure 校验 RMBLAST_DIR 需同时存在
> rmblastn / dustmasker / makeblastdb / blastdbcmd / blastdb_aliastool /
> blastn 六个可执行文件（`RepeatMaskerConfig::validateParam` 逐个检查
> `-x`），而 CBP 只装了前两个。合并目录里 rmblastn+makeblastdb 软链 CBP 版
> （真正使用的引擎），其余四个软链官方包、仅为 configure 校验占位——
> RepeatMasker 4.2.4 运行时只调 rmblastn 和 makeblastdb（源码核实：
> dustmasker/blastdbcmd/blastdb_aliastool/blastn 除校验名单外无任何引用）。
> **注意**：四个占位软链指向官方预编译包（glibc ≥2.29），在 CentOS 7 上
> 本身跑不起来，只是不会被调用；若想彻底干净，可把这四个换成同年代老构建
> （如 blast+ 2.2.28）或 bioconda 包里的对应二进制。搬到 CentOS 7 时软链
> 目标路径需保持一致或改硬拷贝。

### 重新 configure（若以后目录再移动）

```bash
cd /home/wangq/Scripts/pgr/RepeatMasker
perl ./configure -perlbin "$(which perl)" \
  -trf_prgm /home/wangq/.cbp/bin/trf \
  -rmblast_dir /home/wangq/Scripts/pgr/RepeatMasker/rmblast-cbp-bin \
  -default_search_engine rmblast
```

configure 期间回答 "Configure FamDB now?" 为 `n`（当前不需要物种库）。

## 冒烟测试（TnCentral 库）

MG1655 前 100 kb（`tests/genome/mg1655.fa.gz`），`-pa 8`：

```bash
zcat ~/data/repeats/tncentral.fa.gz > /tmp/rmtest/tncentral.fa
RepeatMasker mg1655_chunk.fa -lib /tmp/rmtest/tncentral.fa \
  -pa 8 -e rmblast -dir /tmp/rmtest/out2
```

结果：3.12% 被遮蔽；检出 IS621、IS186B、IS1A、ISPpu12、ISEc39、
ISSoEn2、Tn7243 等真实 IS 序列，与 MG1655 已知内容吻合。

## 两个实测坑

1. **`-lib` 不接受 gzip**：直接把 `tncentral.fa.gz` 传给 `-lib` 会失败，
   makeblastdb 报 "Input doesn't start with a defline"。
   必须先 `zcat` 解压成 `.fa`。
2. **TnCentral 源库有 24/6093 条记录格式瑕疵**：部分序列行开头黏着
   accession 前缀（如 `In1223` 的序列以 `NX784502...` 开头，少数以
   `_PAJ...` 开头）。RepeatMasker 能容忍（按非法字符处理），
   但正式做金标准比对前建议把这些前缀清掉。
3. **RM 会把 .fa.gz 输入解压到文件旁**：RepeatMasker 对 `.fa.gz` 输入
   自动 `gunzip -c file.gz > file`（写到输入同目录，RepeatMasker:748-755）。
   **不要在仓库内直接对 .gz 跑 RM**，否则会在源码目录留下一份解压副本
  （2026-08-07 实测：`tests/genome/mg1655.fa` 因此反复出现）。做法：
  先 zcat 解压到 /tmp 再喂 RM，或用软链；若文件旁已存在同名解压版，
  RM 会直接报错退出。

## CentOS 7（glibc 2.17）部署兼容性（2026-08-07）

**问题**：NCBI 官方 `rmblast-2.14.1+-x64-linux-GLIBC_2.31.tar.gz`
预编译包要求的最高 glibc 符号为 **GLIBC_2.29**（本机 readelf 实测），
CentOS 7 只有 glibc 2.17，直接跑不起来。

**我们代码的真实依赖**（比之前判断宽松）：
- `parse_tab_row` 只消费 rmblastn tab 输出的 qseqid/qstart/qend 三列；
  18 列 outfmt 里 2.13+ 才新增的 kdiv/cpg 等列**没有被使用**。
- RepeatMasker 4.2.4 自身对 rmblastn <2.13 也只是退回 legacy 解析
  （NCBIBlastSearchEngine.pm setPathToEngine 里的特性开关），并非硬性拒绝。

**可选方案**（按推荐顺序）：
1. **CBP 安装的 rmblast 2.14.1（`~/.cbp/bin/rmblastn` + `makeblastdb`）**：
   2026-08-07 实测，版本 2.14.1、要求的最高 glibc 符号仅 **GLIBC_2.16**
   （官方预编译包是 2.29），CentOS 7（2.17）可直接运行。用它跑 60 kb MG1655
   片段 × TnCentral：63 条原始命中 / 25 条去重区间，与官方 2.14.1 结果
   完全一致；18 列 outfmt 与 v5 库格式均为正式行为，代码零改动。
2. **bioconda rmblast 2.14.1**：linux-64 conda 包按 glibc 2.17 兼容构建，
   CentOS 7 可直接运行，版本与本地验证完全一致、零结果差异。
   服务器装 micromamba（单静态二进制，无需 root）：
   `micromamba create -p ~/rmblast-env -c conda-forge -c bioconda rmblast=2.14.1`，
   然后 `pgr rept masker ... --rmblast-dir ~/rmblast-env/bin`。
3. **CentOS 7 源码编译 2.14.1**：需要 devtoolset-8+（C++14）与 boost、
   zlib/bzip2 开发包，编译耗时长，只作无 conda 时的备选。

**pgr 本体**：glibc 兼容性由项目发布流程保证（`.github/workflows/publish.yml`
用 `cargo zigbuild` 交叉到 glibc 2.17），此处不赘述。TRF 4.09 为静态二进制。

## 相关笔记

* [design/masker.md](../design/masker.md)：`pgr rept masker`
  实现设计（参数表、TRF 两阶段、验证结果）
* [design/repeat-masking.md](../design/repeat-masking.md)：重复标记总体方案与
  RepeatMasker 源码梳理（附录 A）
* [ecoli-repeats.md](../ecoli-repeats.md) §2.7/§2.8：RepeatMasker 金标准核对与
  masker 复刻对拍

## 参考

- 官方安装页：https://www.repeatmasker.org/RepeatMasker/（依赖、configure 流程、`-lib` 声明）
- RMBlast 下载：https://www.repeatmasker.org/rmblast/
- GitHub README："You can use it immediately with a custom library (`-lib mylib.fa`)"
- GitHub issue #289：额外库（Dfam partition 等）均通过 `-lib` 传入
