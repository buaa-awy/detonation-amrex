# 物理验证流程

本项目的二维爆震案例需要与原始 Python 实现和理论基准逐层对照。推荐按以下
顺序验证，每次代码修改都从 `t=0` 开始。

## 1. 化学基准

1. 运行 `./build/chemtest`，确认隐式化学迭代收敛；
2. 在同一机理和初值下，用 Cantera 计算 0D 点火延迟；
3. 比较 H2、O2、OH、H2O 和温度随时间的变化；
4. 检查反应速率单位、第三体效率和反应用热力学数据。

## 2. CJ 和平面爆震基准

1. 计算目标初值下的 Chapman–Jouguet 压力和速度；
2. 用平面入射激波初始化未燃区；
3. 检查压力锋与反应锋是否保持耦合；
4. 比较波前速度、峰值压力和平衡产物，而不是只看计算是否结束。

仓库提供了可直接运行的验证输入：

```bash
DETONATION_RUN_DIR=runs/cj-validation \
  ./scripts/run_case.sh examples/h2o2_2d_demo/inputs-cj-validation
```

当前结果的 5–30 μs 波速为 `2830.03 m/s`，相对 Cantera CJ 参考值误差为
`-0.30%`。完整参考状态和剖面见 [`CJ_VALIDATION.md`](CJ_VALIDATION.md)。

## 3. 原 Python 对照

AMReX 长算例使用与原 Python 配置一致的最细网格：

```text
dx = 25 μm
dy = 25 μm
base n_cell = 1500 x 100
refinement = three 2:1 levels
```

从 `t=0` 同时计算 Python 和 AMReX 的 5 μs、10 μs、20 μs 和 30 μs 场，
比较中心线 `p`、`T`、`Y_H2`、`Y_O2`，并检查波前位置和峰值。

## 4. 二维网格协调性

为了研究胞格结构，再计算接近正方形的最细网格，例如：

```text
25 μm x 25 μm
50 μm x 50 μm
100 μm x 100 μm
```

各向异性网格可以用来严格复现原 Python 算例，但最终结论应通过网格收敛性
检查。若 `dx` 和 `dy` 相差较大，必须检查压力、温度、反应进度和胞格尺寸
是否出现方向偏差。

## 5. checkpoint 版本

checkpoint 头部写入求解器版本。修改实现后，请从 `t=0` 生成新的 checkpoint，
并确认新 checkpoint 的版本号与当前程序一致。

只有确认旧 checkpoint 的物理状态仍兼容时，才显式设置：

```text
det.allow_legacy_restart = 1
```
