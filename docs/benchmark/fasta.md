# `pgr fa` vs `faops`

Benchmarks comparing the Rust implementation (`pgr fa`) against the C implementation (`faops`).

## Preparation

Ensure `pgr` is built in release mode for accurate benchmarking:

Install dependencies:

```shell
cbp install faops
cbp install hyperfine
cargo install neofetch
```

## System info

Benchmarks were run on the following system:

**Wangq@BD795 (Ryzen 9 7945HX)**

```text
OS: Ubuntu 24.04.3 LTS x86_64
Kernel: 6.6.87.2-microsoft-standard-WSL2
CPU: AMD Ryzen 9 7945HX with Radeon Graphics (8) 2.50 GHz
Memory: 1.3 GiB / 45.9 GiB
```

<details>
<summary>Full System Details</summary>

```text
            .-/+oossssoo+\-.               wangq@BD795
        ´:+ssssssssssssssssss+:`           -------
      -+ssssssssssssssssssyyssss+-         OS: Ubuntu 24.04.3 LTS x86_64
    .ossssssssssssssssssdMMMNysssso.       Kernel: 6.6.87.2-microsoft-standard-WSL2
   /ssssssssssshdmmNNmmyNMMMMhssssss\      Uptime: 5 days, 3 hours, 3 mins
  +ssssssssshmydMMMMMMMNddddyssssssss+     Packages: 479 (dpkg)
 /sssssssshNMMMyhhyyyyhmNMMMNhssssssss\    Shell: bash 5.2.21
.ssssssssdMMMNhsssssssssshNMMMdssssssss.   Display(rdp-0): 3840x2160 @ 60Hz (as 3840x2160)
+sssshhhyNMMNyssssssssssssyNMMMysssssss+   Terminal: Windows Terminal
ossyNMMMNyMMhsssssssssssssshmmmhssssssso   Disk(/): 78.3 GiB / 1006.9 GiB (8%)
ossyNMMMNyMMhsssssssssssssshmmmhssssssso   CPU: AMD Ryzen 9 7945HX with Radeon Graphics (8) 2.50 GHz
+sssshhhyNMMNyssssssssssssyNMMMysssssss+   Memory: 1.3 GiB / 45.9 GiB
.ssssssssdMMMNhsssssssssshNMMMdssssssss.   Local IP: 172.17.86.235
 \sssssssshNMMMyhhyyyyhdNMMMNhssssssss/    Network: eth0 (172.17.86.235)
  +sssssssssdmydMMMMMMMMddddyssssssss+     Locale: C.UTF-8
   \ssssssssssshdmNNNNmyNMMMMhssssss/
    .ossssssssssssssssssdMMMNysssso.
      -+sssssssssssssssssyyyssss+-
        `:+ssssssssssssssssss+:`
            .-\+oossssoo+/-.
```
</details>

## `pgr fa size`

* ufasta

```shell
hyperfine --warmup 10 --export-markdown size.md.tmp \
    -n "pgr fa size .fa" \
    'cat tests/fasta/ufasta.fa | pgr fa size stdin > /dev/null' \
    -n "faops size .fa" \
    'cat tests/fasta/ufasta.fa | faops size stdin > /dev/null' \
    -n "pgr fa size .fa.gz" \
    'pgr fa size tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "faops size .fa.gz" \
    'faops size tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "pgr 2bit size .2bit" \
    'pgr 2bit size tests/fasta/ufasta.2bit > /dev/null'

cat size.md.tmp

```

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `pgr fa size .fa` | 2.9 ± 0.3 | 2.4 | 4.7 | 1.19 ± 0.21 |
| `faops size .fa` | 2.6 ± 0.3 | 2.0 | 6.4 | 1.07 ± 0.18 |
| `pgr fa size .fa.gz` | 3.0 ± 0.3 | 2.5 | 4.8 | 1.23 ± 0.20 |
| `faops size .fa.gz` | 2.4 ± 0.3 | 1.8 | 4.4 | 1.00 |
| `pgr 2bit size .2bit` | 10.7 ± 0.7 | 9.6 | 14.3 | 4.41 ± 0.63 |

* mg1655

```shell
hyperfine --warmup 10 --export-markdown size.md.tmp \
    -n "pgr fa size .fa.gz" \
    'pgr fa size tests/genome/mg1655.fa.gz > /dev/null' \
    -n "faops size .fa.gz" \
    'faops size tests/genome/mg1655.fa.gz > /dev/null' \
    -n "pgr 2bit size .2bit" \
    'pgr 2bit size tests/genome/mg1655.2bit > /dev/null'

cat size.md.tmp

```

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `pgr fa size .fa.gz` | 19.0 ± 0.6 | 17.8 | 21.8 | 6.87 ± 0.64 |
| `faops size .fa.gz` | 36.7 ± 0.7 | 35.6 | 40.4 | 13.23 ± 1.19 |
| `pgr 2bit size .2bit` | 2.8 ± 0.2 | 2.4 | 4.2 | 1.00 |

* mg1655 protein

```shell
hyperfine --warmup 10 --export-markdown size.md.tmp \
    -n "pgr fa size .pro.fa.gz" \
    'pgr fa size tests/genome/mg1655.pro.fa.gz > /dev/null' \
    -n "faops size .pro.fa.gz" \
    'faops size tests/genome/mg1655.pro.fa.gz > /dev/null'

cat size.md.tmp

```

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `pgr fa size .pro.fa.gz` | 11.7 ± 0.5 | 10.8 | 14.3 | 1.00 |
| `faops size .pro.fa.gz` | 19.1 ± 0.6 | 18.2 | 24.6 | 1.63 ± 0.09 |

## `pgr fa some`

```shell
hyperfine --warmup 10 --export-markdown some.md.tmp \
    -n "pgr fa some" \
    'pgr fa some tests/fasta/ufasta.fa.gz tests/fasta/list.txt > /dev/null' \
    -n "faops some" \
    'faops some tests/fasta/ufasta.fa.gz tests/fasta/list.txt stdout > /dev/null' \
    -n "pgr fa some -i" \
    'pgr fa some -i tests/fasta/ufasta.fa.gz tests/fasta/list.txt > /dev/null' \
    -n "faops some -i" \
    'faops some -i tests/fasta/ufasta.fa.gz tests/fasta/list.txt stdout > /dev/null'

cat some.md.tmp

```

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `pgr fa some` | 4.1 ± 0.3 | 3.3 | 5.6 | 1.08 ± 0.14 |
| `faops some` | 3.9 ± 0.4 | 3.3 | 5.4 | 1.04 ± 0.16 |
| `pgr fa some -i` | 4.6 ± 0.5 | 3.6 | 8.1 | 1.21 ± 0.19 |
| `faops some -i` | 3.8 ± 0.4 | 3.2 | 6.2 | 1.00 |

## `pgr fa n50`

```shell
hyperfine --warmup 10 --export-markdown n50.md.tmp \
    -n "pgr fa n50 .gz" \
    'pgr fa n50 tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "faops n50 .gz" \
    'faops n50 tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "pgr fa n50 -E -S -A" \
    'pgr fa n50 -E -S -A tests/fasta/ufasta.fa > /dev/null' \
    -n "faops n50 -E -S -A" \
    'faops n50 -E -S -A tests/fasta/ufasta.fa > /dev/null'

cat n50.md.tmp

```

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `pgr fa n50 .gz` | 3.0 ± 0.3 | 2.5 | 5.7 | 1.29 ± 0.19 |
| `faops n50 .gz` | 2.4 ± 0.3 | 2.0 | 3.9 | 1.00 |
| `pgr fa n50 -E -S -A` | 3.4 ± 0.2 | 3.0 | 5.0 | 1.44 ± 0.19 |
| `faops n50 -E -S -A` | 2.4 ± 0.2 | 2.1 | 3.9 | 1.03 ± 0.15 |

## Conclusion

The `pgr` implementation demonstrates competitive performance compared to the highly optimized C implementation `faops`.

*   **Small datasets**: `pgr` shows a slight overhead (1.1x - 1.4x). This is likely due to the startup time of the larger executable (~3MB vs ~400KB for `faops`) and the initialization of Rust's runtime and command-line parser.
*   **Larger compressed datasets**: `pgr` significantly outperforms `faops` (e.g., ~2x faster for `mg1655.fa.gz` and ~1.6x faster for `mg1655.pro.fa.gz`), highlighting the efficiency of Rust's Gzip handling and I/O.
*   **2bit format**: `pgr` provides extremely fast metadata retrieval.
