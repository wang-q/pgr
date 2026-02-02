# `hnsm` and `faops`

Benchmarks between C and Rust implementations

```shell
cbp install faops

cbp install hyperfine

cargo install neofetch

```

## System info

* Ryzen 7 5800

```text
wangq@R7
--------
OS: Ubuntu 20.04.5 LTS on Windows 10 x86_64
Kernel: 5.15.153.1-microsoft-standard-WSL2
Uptime: 16 days, 6 hours, 30 mins
Packages: 1408 (dpkg), 191 (brew), 5 (snap)
Shell: bash 5.0.17
Theme: Adwaita [GTK3]
Icons: Adwaita [GTK3]
Terminal: Windows Terminal
CPU: AMD Ryzen 7 5800 (16) @ 3.393GHz
GPU: fb45:00:00.0 Microsoft Corporation Device 008e
Memory: 562MiB / 32030MiB

```

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

## `pgr fa size`

```shell
hyperfine --warmup 10 --export-markdown size.md.tmp \
    -n "pgr fa size .fa" \
    'cat tests/fasta/ufasta.fa | pgr fa size stdin > /dev/null' \
    -n "faops size .fa" \
    'cat tests/fasta/ufasta.fa | faops size stdin > /dev/null' \
    -n "pgr fa size .fa.gz" \
    'pgr fa size tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "faops size .fa.gz" \
    'faops size tests/fasta/ufasta.fa.gz > /dev/null'

cat size.md.tmp

```

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `pgr fa size .fa` | 3.0 ± 0.4 | 2.3 | 5.5 | 1.27 ± 0.21 |
| `faops size .fa` | 2.7 ± 0.2 | 2.2 | 4.2 | 1.13 ± 0.14 |
| `pgr fa size .fa.gz` | 3.0 ± 0.3 | 2.6 | 5.2 | 1.29 ± 0.17 |
| `faops size .fa.gz` | 2.4 ± 0.2 | 1.9 | 3.8 | 1.00 |

## `hnsm some`

```shell
hyperfine --warmup 10 --export-markdown some.md.tmp \
    -n "hnsm some" \
    'hnsm  some tests/fasta/ufasta.fa.gz tests/fasta/list.txt > /dev/null' \
    -n "faops some" \
    'faops some tests/fasta/ufasta.fa.gz tests/fasta/list.txt stdout > /dev/null' \
    -n "hnsm some -i" \
    'hnsm  some -i tests/fasta/ufasta.fa.gz tests/fasta/list.txt > /dev/null' \
    -n "faops some -i" \
    'faops some -i tests/fasta/ufasta.fa.gz tests/fasta/list.txt stdout > /dev/null'

cat some.md.tmp

```

| Command         | Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:----------------|----------:|---------:|---------:|------------:|
| `hnsm some`     | 3.9 ± 0.3 |      3.3 |      6.0 | 1.01 ± 0.12 |
| `faops some`    | 4.1 ± 0.3 |      3.6 |      5.9 | 1.06 ± 0.12 |
| `hnsm some -i`  | 3.9 ± 0.3 |      3.4 |      5.8 |        1.00 |
| `faops some -i` | 4.1 ± 0.4 |      3.6 |      6.3 | 1.07 ± 0.13 |

## `hnsm n50`

```shell
hyperfine --warmup 10 --export-markdown n50.md.tmp \
    -n "hnsm n50 .gz" \
    'hnsm  n50 tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "faops n50 .gz" \
    'faops n50 tests/fasta/ufasta.fa.gz > /dev/null' \
    -n "hnsm n50 -E -S -A" \
    'hnsm  n50 -E -S -A tests/fasta/ufasta.fa > /dev/null' \
    -n "faops n50 -E -S -A" \
    'faops n50 -E -S -A tests/fasta/ufasta.fa > /dev/null'

cat n50.md.tmp

```

| Command              | Mean [ms] | Min [ms] | Max [ms] |    Relative |
|:---------------------|----------:|---------:|---------:|------------:|
| `hnsm n50 .gz`       | 2.7 ± 0.2 |      2.3 |      3.7 | 1.00 ± 0.09 |
| `faops n50 .gz`      | 2.7 ± 0.2 |      2.3 |      3.6 |        1.00 |
| `hnsm n50 -E -S -A`  | 3.2 ± 0.2 |      2.9 |      4.2 | 1.22 ± 0.10 |
| `faops n50 -E -S -A` | 2.8 ± 0.2 |      2.4 |      3.8 | 1.04 ± 0.09 |
