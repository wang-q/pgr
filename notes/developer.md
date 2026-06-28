# 开发者指南

本文档旨在为 `pgr` 的开发者提供技术背景、架构设计思路以及未来演进路线。

## changelog

```bash
# 用最新的 tag 作为起点（当前为 v0.2.0）
git tag | sort -V | tail -1   # 查看最新 tag
git log v0.2.0..HEAD > gitlog.txt
git diff v0.2.0 HEAD -- "*.rs" "*.md" > gitdiff.txt
```

## code coverage

```bash
rustup component add llvm-tools
cargo install cargo-llvm-cov

# 生成覆盖率报告
cargo llvm-cov
```

使用 `cargo llvm-cov` 生成覆盖率报告，找出需要提升测试覆盖率的代码路径，供我分析。

XXX 的测试覆盖度不高，使用 `cargo llvm-cov` 生成覆盖率报告，找出需要提升的地方.

为这些地方，添加单元测试与整合测试

为刚才的修改，添加单元测试与整合测试

## WSL

```bash
mkdir -p /tmp/cargo
export CARGO_TARGET_DIR=/tmp/cargo
cargo build
```
