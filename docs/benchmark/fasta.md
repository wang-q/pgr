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

The `pgr` implementation is currently slightly slower (1.1x - 1.4x) than the highly optimized C implementation `faops`. This is expected as `pgr` prioritizes safety and features over raw C-level optimization in its initial release. The overhead of command-line argument parsing and safety checks in Rust contributes to the difference, especially for small execution times in the millisecond range.
