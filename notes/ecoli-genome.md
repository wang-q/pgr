# E. coli genome

使用 MG1655、Sakai 和 SE11 三株大肠杆菌基因组，演示 `pgr` 的 UCSC pipeline 和比对可视化。
MG1655 是 K-12 实验室菌株，Sakai 是 O157:H7 致病菌株，两者基因组大小约 4.6 Mb 和 5.5 Mb；
SE11 是 O152:H28 共生株，除 4.8 Mb 染色体外还携带 6 个质粒，用于多 replicon（多染色体）测试。

## Download

```bash
# MG1655 (K-12)
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz |
    gzip -dc |
    pgr fa filter stdin --simplify |
    pgr fa gz stdin -o tests/genome/mg1655.fa.gz

# Sakai (O157:H7)
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz |
    gzip -dc |
    pgr fa filter stdin --simplify |
    pgr fa gz stdin -o tests/genome/sakai.fa.gz

# SE11 (O152:H28, 共生株; RefSeq GCF_000010385.1, ASM1038v1)
# 1 染色体 + 6 质粒 (pSE11-1 ~ pSE11-6)，共 7 个 replicon
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/010/385/GCF_000010385.1_ASM1038v1/GCF_000010385.1_ASM1038v1_genomic.fna.gz |
    gzip -dc |
    pgr fa filter stdin --simplify |
    pgr fa gz stdin -o tests/genome/se11.fa.gz
```

## Multi-replicon test data (SE11)

SE11 是携带 6 个质粒（pSE11-1 ~ pSE11-6，4 kb ~ 100 kb）的共生菌株，加上 4.8 Mb 染色体
共 7 个 replicon，适合验证 pgr 在多染色体 / 小 contig 上的行为：

- chain/net 的多染色体排序、`net split`、maf 合并；
- 反向 PSL（以多 replicon 基因组作 target）的字节级一致性验证；
- 小质粒（最小 4 kb）在链化中的边界处理。

下载后确认 replicon 数（应为 7）：

```bash
gzip -dc tests/genome/se11.fa.gz | grep -c '^>'
```

SE11 与 Sakai 的比对 LAV 已预存为 `tests/genome/se11-sakai.lastz.lav`（SE11 作 target，
约 3.8 MB），供 `scripts/verify-ucsc-pipeline.sh` 的反向多染色体验证使用，免去重跑
lastz。生成方式：lastz 的 `multiple` action 与 LAV 格式不兼容，因此按 SE11 每个
replicon 单独跑 lastz（target 单序列），再按 replicon 顺序合并 LAV：

```bash
# 对 SE11 每条序列（染色体 + 6 质粒）分别：
lastz <(gzip -dcf se11.fa.gz | awk '/^>/{p=($0 ~ /NC_011415/)} p') \
    <(gzip -dcf tests/genome/sakai.fa.gz) > NC_011415.lav
# 然后按 replicon 顺序 cat 合并：
cat NC_011415.lav NC_011419.lav NC_011413.lav NC_011416.lav \
    NC_011407.lav NC_011408.lav NC_011411.lav > se11-sakai.lastz.lav
```

实测 7 个 replicon 中 6 个与 Sakai 有链（仅 4 kb 的 pSE11-6 无比对），链数
13999+132+55+32+2+1 块。

## Alignment and visualization

```bash
# Pairwise alignment with FastGA
FastGA -v -pafx tests/genome/sakai.fa.gz tests/genome/mg1655.fa.gz > tmp.paf
FastGA -v -psl tests/genome/sakai.fa.gz tests/genome/mg1655.fa.gz > tmp.psl

# UCSC pipeline: chain → net → axt → maf
pgr pl ucsc -t="" tests/genome/mg1655.fa.gz tests/genome/sakai.fa.gz tmp.psl > tmp.chain.maf
pgr pl ucsc --syn -t="" tests/genome/mg1655.fa.gz tests/genome/sakai.fa.gz tmp.psl > tmp.syn.maf

# LASTZ alignment
# 1. Generate LAV with lastz (~2 min).  A pre-computed copy is saved in the
#    repo as tests/genome/mg1655-sakai.lastz.lav, so this step can be skipped.
lastz <(gzip -dcf tests/genome/mg1655.fa.gz) <(gzip -dcf tests/genome/sakai.fa.gz) \
    > tmp.lastz.lav
# 2. Convert LAV to PSL (use the saved LAV, or tmp.lastz.lav if regenerated).
lavToPsl tests/genome/mg1655-sakai.lastz.lav tmp.lastz.psl
pgr pl ucsc --syn -t="" tests/genome/mg1655.fa.gz tests/genome/sakai.fa.gz tmp.lastz.psl > tmp.lastz.maf

# Dotplot visualization
wgatools dotplot -f paf tmp.paf > tmp.html
wgatools dotplot tmp.chain.maf > tmp.chain.html
wgatools dotplot tmp.syn.maf > tmp.syn.html
wgatools dotplot tmp.lastz.maf > tmp.lastz.html
```

| ![paf](../images/paf.png) | ![chain](../images/chain.png) |
|:--------------------------:|:------------------------------:|
|            paf             |             chain              |

| ![syn](../images/syn.png) | ![lastz](../images/lastz.png) |
|:--------------------------:|:------------------------------:|
|            syn             |             lastz              |
