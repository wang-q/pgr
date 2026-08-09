# pgr kmer 命令吞吐 sanity（2026-08-09）

> AGENTS.md 要求：功能性新功能只需正确性测试 + 吞吐 sanity。本页记录
> `pgr kmer` 七个子命令（table/profile/hist/gc/qhist/qcheck/gsize）在
> release 构建下的吞吐与内存。机器：AMD Ryzen 9 7945HX（x86_64），
> stable 1.97，release。随机输入：5 Mb 基因组 FASTA（50×100 kb 随机
> DNA，perl `substr("ACGT", int(rand(4)), 1)` 生成）、1 万条 100 bp
> 随机 FASTQ reads（~2.1 MB，质量全 I，阈值 38 自动检测；reads 从
> /dev/urandom → base64 → tr 过滤的 1 Mb DNA 切块）。
>
> 注意：早期版本用 perl `("ACGT")[int(rand(4))]` 生成——那是**列表索引**
> 不是字符索引（下标 >0 返回空），序列退化为 ACGT 块拼接（唯一 k-mer
> 极少、命令更快），本页已用正确生成方式重测。

## 基因组型（5 Mb，k=17）

| 命令 | wall | max RSS | 说明 |
|---|---:|---:|---|
| `kmer table` | 0.29 s | 191 MB | ~5.0 M 唯一 k-mer，.pkt 紧凑落盘 |
| `kmer hist` | 0.18 s | 191 MB | 建表 + 聚合（.hist 固定 256 KB） |
| `kmer gc` | 0.30 s | 191 MB | 矩阵计算同 hist；随机低覆盖数据报
  "no maximal peak" 属预期（count≥2 无峰，KatGC 同行为），`-X` 可显式
  覆盖 |
| `kmer profile` | 0.41 s | 358 MB | self profile（排序合并，表 + profile
  双份内存） |

## Reads 型（1 万 × 100 bp，k=17，阈值 38）

| 命令 | wall | max RSS | 说明 |
|---|---:|---:|---|
| `kmer qhist` | 0.07 s | 67 MB | 质量偏置表 + 直方图 |
| `kmer qcheck` | 0.07 s | 65 MB | 建表 + 逐 read 判定（rayon 并行，read 间
  独立）；随机 reads 约 3.1%（308/10000）被判定有错（随机 k-mer 低覆盖
  预期）。并行化前单线程 0.50 s → ~7× |
| `kmer gsize` | 0.03 s | 44 MB | 随机 reads 无覆盖结构，peak=1 属预期
  （gsize 的峰值估计需要真实覆盖度数据，合成 30× 1 kb 测试见单测） |

## 备注

* 内存大头是 `KmerTable`（u128 key + u32 count，~20 B/唯一 k-mer）与
  profile（u16/位置）；5 Mb 细菌规模无压力，50 Mb 真菌 ~1 GB 量级
  （kmer.md §3.2 设计值）。
* qcheck 逐 read 判定为 O(read × k) 滚动查询（partition_point）；read 间
  独立，判定阶段已 rayon 并行（表只读），输出按输入序。
