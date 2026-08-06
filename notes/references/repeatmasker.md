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

## 参考

- 官方安装页：https://www.repeatmasker.org/RepeatMasker/（依赖、configure 流程、`-lib` 声明）
- RMBlast 下载：https://www.repeatmasker.org/rmblast/
- GitHub README："You can use it immediately with a custom library (`-lib mylib.fa`)"
- GitHub issue #289：额外库（Dfam partition 等）均通过 `-lib` 传入
