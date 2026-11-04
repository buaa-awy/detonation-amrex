# 二维 H2/O2 爆震示范算例

两个算例使用相同的 9 组分、21 反应化学机理，以及约 25 μm 的最细网格。
它们用于演示程序用法，不代表已经完成网格无关性验证或燃烧模型验证。

## 快速检查算例

`inputs` 的计算域为 0.12 m x 0.01 m，停止时间为 5 μs。它适合在安装完成后
先检查 AMReX、CUDA、输入输出和监视器是否正常。

`inputs-smoke` 是固定 100 步、无 ParaView/Tecplot 输出的确定性回归输入，
用于判断版本更新是否保持基本数值结果。推荐先运行：

```bash
./scripts/run_smoke_test.sh
```

```bash
cd /path/to/detonation-amrex-github
DETONATION_RUN_DIR=runs/quick \
  ./scripts/run_case.sh examples/h2o2_2d_demo/inputs
```

## 长时间算例

`inputs-long` 的计算域为 0.30 m x 0.02 m，停止时间为 75 μs。笔记本电脑 GPU
通常需要数小时；建议先用快速算例完成环境和输出检查。

```bash
DETONATION_RUN_DIR=runs/long \
  ./scripts/run_case.sh examples/h2o2_2d_demo/inputs-long
```

## 初始条件

- 点火区：`x <= det.x_igniter`
- 点火区温度、压力：`1800 K`、`151987.5 Pa`
- 未燃混合气温度、压力：`300 K`、`8000 Pa`
- 初始组分：`Y_H2=0.1119067440`、`Y_O2=0.8880932560`
- 横向温度扰动：12%，4 个正弦波周期
- 左边界：反射壁面
- 右边界：一阶外推
- 上下边界：周期边界

修改算例时，建议先复制输入文件，不要直接覆盖示范算例。完整参数说明见
[`docs/INPUT_REFERENCE.md`](../../docs/INPUT_REFERENCE.md)。
