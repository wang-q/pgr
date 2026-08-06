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
- 搜索引擎：RMBlast 2.14.1，`/home/wangq/share/rmblast/bin`（默认引擎）。
- TRF：`/home/wangq/.cbp/bin/trf`（系统已有）。
- FamDB：未配置（因此只有 `-lib` 模式可用；`-species` 不可用）。
- 已通过 `perl ./configure` 完成配置（RMBLAST_DIR 指向 `~/share/rmblast/bin`），
  `./RepeatMasker -h` 可运行。RMBlast 的 tar 包仍在 `~/Downloads/`，/tmp 里
  的解压副本可删。

### 重新 configure（若以后目录再移动）

```bash
cd /home/wangq/Scripts/pgr/RepeatMasker
perl ./configure -perlbin "$(which perl)" \
  -trf_prgm /home/wangq/.cbp/bin/trf \
  -rmblast_dir /home/wangq/share/rmblast/bin \
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
- 因此若把 outfmt 收敛到 `6 qseqid qstart qend`（不影响命中集合），
  理论上任意版本 rmblastn 都能配；但老引擎（2.2.28 等）与 2.14.1
  在同样参数下的命中集合可能有差异，金标准对齐会打折扣。

**可选方案**（按推荐顺序）：
1. **bioconda rmblast 2.14.1**：linux-64 conda 包按 glibc 2.17 兼容构建，
   CentOS 7 可直接运行，版本与本地验证完全一致、零结果差异。
   服务器装 micromamba（单静态二进制，无需 root）：
   `micromamba create -p ~/rmblast-env -c conda-forge -c bioconda rmblast=2.14.1`，
   然后 `pgr rept masker ... --rmblast-dir ~/rmblast-env/bin`。
2. **CentOS 7 源码编译 2.14.1**：需要 devtoolset-8+（C++14）与 boost、
   zlib/bzip2 开发包，编译耗时长，只作无 conda 时的备选。
3. **老预编译 `ncbi-rmblastn-2.2.28-x64-linux.tar.gz`**（NCBI FTP
   `blast/executables/rmblast/2.2.28/`）：2015 年构建，面向当时主流系统
   （glibc 2.12/2.14），CentOS 7 大概率可跑、体积小；但需把 outfmt
   收到旧版支持的列，且命中集合与 2.14.1 有差异。

**别忘了 pgr 本体**：Rust 二进制默认链接构建机的 glibc。若在 Ubuntu
22.04 上编译，pgr 拿到 CentOS 7 同样会报 GLIBC_2.34 not found。
要整套上 CentOS 7 需：在服务器上 rustup 编译，或本机
`cargo zigbuild --release --target x86_64-unknown-linux-gnu.2.17`
交叉编译，或直接 musl 静态。TRF 4.09 为静态二进制，无此问题。

## 参考

- 官方安装页：https://www.repeatmasker.org/RepeatMasker/（依赖、configure 流程、`-lib` 声明）
- RMBlast 下载：https://www.repeatmasker.org/rmblast/
- GitHub README："You can use it immediately with a custom library (`-lib mylib.fa`)"
- GitHub issue #289：额外库（Dfam partition 等）均通过 `-lib` 传入
