# Range 字符串解析基准：正则 vs 手写扫描器

> 目的：`libs::ds::Range` 的 `from_str` 用正则解析区间字符串（如
> `S288c.I(-):27070-29557`），对正则性能存疑，验证手写解析器能否等价替代
> 并提速。

## 方案（已切换）

* **生产路径**：`Range::from_str` → `decode`，已切为手写逐字节扫描器；
  复刻正则的全部语义：非锚定最左匹配、贪婪 name/chr/strand、`start` 缺
  `-end` 时 `end = start`、显式 `end = 0` 视为缺失（`c:911_0` → 911-911）、
  无匹配时回退为第一个空白 token；唯一差异：数字溢出 i32 时返回 0 而不是
  panic（修复了正则路径的 `parse::<i32>().unwrap()` panic）。
* **正则保留为文档**：原正则原文与语义说明在 `src/libs/ds/range.rs` 的
  模块文档里；测试模块保留正则解码器作为对拍 oracle，基准文件里也留了
  一份作为对比基线。

等价性由 `regex_and_manual_decoders_agree` 保证：固定语料（含
`foo I:1-100`、`a.b.c:1-2`、`1:-100` 等边界）+ 2 万条随机 fuzz 逐字段对拍
一致。

## 执行

```bash
cargo bench --offline --bench range_parse_benchmark
```

语料为 17 条真实格式混合（普通 / 带链向 / 带物种前缀 / 单坐标 /
斜杠下划线 contig / 回退用例）。

## 结果（median）

| 方案 | 17 条语料整体 |
| --- | ---: |
| 正则（基线） | 7.651 µs（复跑 5.066 µs） |
| 手写 `from_str`（生产） | 0.859 µs（复跑 1.015 µs） |

手写版快约 **5–9 倍**（单条约 50–60 ns vs 300–450 ns；两次 criterion 运行
中位数有波动），已作为生产实现。
