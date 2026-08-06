# SIMD / HV / Jaccard 基准（自 hnsm 迁移）

> 2026-08-06 自 `hnsm/benches/{simd,jaccard,hd}.rs` 迁移至 `benches/`，
> 对齐 pgr 现役实现（`libs::linalg` / `libs::hv`，wide 1.6.0）。
> 机器：AMD Ryzen 9 7945HX（32 核，x86_64），stable 1.97.0，release 构建。

## 1. SIMD：L2 范数（`benches/simd_benchmark.rs`）

向量长度 10005（8-lane SIMD 的余数分支随之验证）。四个实现：
`map` / `fold` / pgr `linalg::norm_l2`（wide SIMD）/ `nalgebra::DVector::norm`。

| 实现 | 耗时 |
|---|---:|
| norm_map | 5.617 µs |
| norm_fold | 5.610 µs |
| **norm_simd_pgr** | **723.9 ns** |
| norm_nalgebra | 1.354 µs |

SIMD 比 map/fold 快 **~7.8×**，比 nalgebra 快 ~1.9×。
历史（Ryzen 7 8745HS）：map 6.31 µs / simd 810 ns / nalgebra 1.56 µs——趋势一致。

## 2. HV 编码（`benches/hv_benchmark.rs`，hv_d = 4096）

对比 pgr 现役 `hash_hv_bit`（位操作 SIMD）、`hash_hv_i8`（i8 累加 SIMD）
与三个标量 + RNG 对照（RapidRng / StdRng / SmallRng，i16 累加）。
seed 集合固定种子（StdRng seed 42，RapidHashSet）。

### small（1000 seeds）

| 实现 | 耗时 |
|---|---:|
| **hash_hv_i8** | **447.8 µs** |
| hash_hv_bit | 677.1 µs |
| 标量 + RapidRng | 906.9 µs |
| 标量 + SmallRng | 936.0 µs |
| 标量 + StdRng | 1.039 ms |

### medium（10000 seeds）

| 实现 | 耗时 |
|---|---:|
| **hash_hv_i8** | **4.452 ms** |
| hash_hv_bit | 6.818 ms |
| 标量 + RapidRng | 9.150 ms |
| 标量 + SmallRng | 9.376 ms |
| 标量 + StdRng | 10.445 ms |

i8 实现比 bit 实现快 ~1.53×、比最快标量对照快 ~2×。
历史（2026-01-30）：i8 421 µs / 4.20 ms，lib(bit) 670 µs / 6.73 ms——一致。

### 迁移中发现并修复的性能退化

`hash_hv_i8` 初版迁移用 `bytes.map(|b| b as i8 as i32)`（8 元素标量转换）
构造 SIMD 向量，导致 small 447 → 1789 µs（~4.3× 退化，medium 同样）。
修复：wide 无 u8→i8→i32 数值 lane 转换，改用
`u8x16 → u16x8::from_u8x16_low → i32x8::from_u16x8` 零扩展链 +
`(x << 24) >> 24` 算术移位还原有符号语义（等价 `b as i8 as i32`）。
`hash_hv_bit` 的 u32→i32 转换改 `bytemuck::cast`（0/1 bit pattern 重解释，零开销）。
正确性由 `test_hash_hv_i8_serial_vs_simd` 等对照测试保证。

## 3. 集合 Jaccard（`benches/jaccard_benchmark.rs`）

每个集合 5005 个 u64（0..u16::MAX），105 个集合，随机取两集合比较。

| 实现 | Jaccard 耗时 |
|---|---:|
| **rapidinlinehash** | **20.55 µs** |
| rapidhash | 20.66 µs |
| BTreeSet | 38.24 µs |
| HashSet（SipHash） | 54.46 µs |
| tinyset | 54.41 µs |
| nohash | 59.98 µs |

rapidhash 系最快，比 std HashSet 快 ~2.6×；intersection/union/access
（btree 38.4/47.6 µs、hashset 55.4/40.8 µs、btree_access 4.7 ns）。
历史（msvc）：rapidhash 22.7 µs 最快、hashset 75.6 µs——结论一致。

## 4. 结论

- SIMD 加速真实有效：norm ~7.8×、HV i8 编码 ~2× over 标量。
- pgr 现役选择（`norm_l2` wide SIMD、`hash_hv_i8`）均为基准最优/次优，
  迁移到 wide 后性能与 std::simd 时代持平。
- 集合场景若日后需要，rapidhash 的 HashSet 比 std 默认快 ~2.6×。

## 5. 遗留疑虑（待处理）

HV 编码矢量化后的提速不如预期明显（相对标量仅 ~1.5–2×，而 norm 有 ~7.8×）。
作者长期有此疑虑（2026-08-06 记录）：8-lane 向量理论上应有更高收益，
瓶颈可能在 RNG 生成（RapidRng 每 chunk 一次）、字节→向量转换或内存带宽，
后续单独深挖（见 `../todo.md` §4）。
