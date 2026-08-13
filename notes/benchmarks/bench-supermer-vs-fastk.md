# supermer 两段计数 vs FastK：端到端对照（2026-08-14）

> 目的：回答「pgr 的 FastK 式 super-mer 两段计数原型（`libs/kmer/supermer.rs`）
> 与参考实现 FastK 本体到底差多少、差在哪」。数据与命令见文末「复现」。

## 1. 环境与输入

* 机器：32 核（`nproc=32`），pgr `--release`（LTO），FastK 用
  `/home/wangq/.cbp/bin/FastK`（CBP 安装版，= 上游 ddea6cf，与仓库
  `FASTK-master` 快照计数核心语义一致，见 `design/kmer.md` §2.3）。
* 输入：mg1655 全基因组（4,641,652 bp）按 stride 7 切 150 bp reads
  → **663,072 条 reads，99,460,800 bp**（起点唯一、~21× 覆盖、无重复读，
  接近真实 Illumina 数据的冗余形态）。
* 端到端命令：FastK `-T32 -t1 -k<k> -N<root> reads.fa`（含读文件、建表、
  写 `.ktab`）vs `pgr kmer table -k<k> -o out.pkt reads.fa`（含读文件、
  建表、写 `.pkt`）。两者 unique k-mer 条目数**逐一对上**（k=31:
  4,554,264；k=100: 4,575,350），语义一致。

## 2. 端到端结果（wall / 峰值 RSS，`/usr/bin/time -v`）

| 工具 | k=31 | k=100 |
| :--- | ---: | ---: |
| FastK | **0.74 s / 411 MB** | **0.95 s / 907 MB** |
| pgr kmer table（直接路径） | 1.74 s / 1.46 GB | 1.60 s / 1.95 GB |
| pgr supermer（lib 计时，不含读文件） | 0.87 s（m=12） | 1.85 s（m=12） |

量级与 `unitig-bucket.md` §3.1 的 G37 实测一致（FastK 快 ~1.7–2.5×、
内存省 2–3.6×）。

## 3. FastK 阶段统计（`-v`）解读

k=31：79,568,640 窗口 → **7,042,193 个 super-mer（平均 11.3 窗口）** →
stage-2 加权 k-mer **17,677,673（savings 4.5×）**；阶段耗时
分发 0.573 s + 排序计数 0.144 s + 合并 0.024 s。

k=100：33,816,672 窗口 → **1,428,537 个 super-mer（平均 23.7 窗口）** →
stage-2 加权 k-mer **31,327,679（savings 仅 1.1×）**；阶段耗时
分发 0.580 s + 排序计数 0.311 s + 合并 0.062 s。
FastK 自行告警：`Too much of the data is in reads on the order of the
k-mer size`。

要点：
* FastK 实际用 **5-mer minimizer**（1024 core prefixes，本例 PAD 未触发），
  span 平均只有 11–24 窗口，比 pgr 原型的默认 m=12 更短；
* **k=100 短读的折叠失效（1.1×）是 FastK 同样存在的事实**——span≈整条读、
  无跨读折叠；它靠窗口总量少（150 bp 读只有 51 窗口）+ 工程效率维持 0.95 s。

## 4. 算法一致性对照（同一输入，pgr supermer 内部统计）

| 指标 | FastK k=31 | pgr m=5 k=31 | FastK k=100 | pgr m=5 k=100 |
| :--- | ---: | ---: | ---: | ---: |
| super-mer 实例 | 7,042,193 | 7,088,908 | 1,428,537 | 1,482,539 |
| stage-2 加权 k-mer | 17,677,673 | 17,959,008 | 31,327,679 | 31,337,222 |

pgr 原型的折叠行为与 FastK **几乎完全一致**（同一输入下 n_records /
n_entries 数值对上），说明算法层面已同构，差距不在 super-mer 划分本身。

## 5. 结论

1. **算法**：pgr 的 super-mer/minimizer 两段计数与 FastK 同构，折叠行为
   一致（k=31 均 ~4.5×，k=100 均 ~1.1×）；
2. **k=100 短读**：折叠失效是 FastK 也承认的（savings 1.1× + 告警）；
   此时两段式对 pgr 是纯开销（1.85 s vs 直接 0.92 s lib 计时），
   **super-mer 不是长 k 短读的解**；
3. **FastK 领先来源**：工程效率（C + 位打包 + 内存 907 MB vs pgr
   1.95 GB + 32 线程并行），不是算法；
4. **pgr 现状**：k=31 的 supermer lib 0.87 s 已接近 FastK 0.74 s；
   直接路径 k=100 lib 0.92 s 已优于 FastK 的折叠后路径在 pgr 侧的表现
   （1.85 s），端到端差距主要在读文件与内存模型。

## 6. 复现

```bash
mkdir -p /tmp/fastk_e2e/sort
python3 gen_reads.py   # 见下；产出 /tmp/fastk_e2e/reads.fa
/usr/bin/time -v /home/wangq/.cbp/bin/FastK -T32 -t1 -k31 -v \
  -N/tmp/fastk_e2e/fk31 -P/tmp/fastk_e2e/sort /tmp/fastk_e2e/reads.fa
/usr/bin/time -v ./target/release/pgr kmer table -k31 \
  -o /tmp/fastk_e2e/pgr31.pkt /tmp/fastk_e2e/reads.fa
# k=100 同理；条目数核对：
#   FastK: python 解析 .ktab stub 的 pindex 末项
#   pgr:   (pkt 大小 - 24) / (key_bytes + 4)
```

`gen_reads.py`：gzip 读 `tests/genome/mg1655.fa.gz`，对首条 ≥4 M 的 contig
按 `range(0, len-150+1, 7)` 写 150-mer FASTA。
