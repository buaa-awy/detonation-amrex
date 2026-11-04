# 从零建立一个二维爆震算例

当前求解器支持一个由两段均匀状态组成的二维 H2/O2 初值问题。点火区位于计算域
左侧，未燃混合气位于右侧；横向可以叠加正弦温度扰动。若需要 H2/空气、其他
燃料、分层混合气或随空间连续变化的初值，必须同时修改计算源码，只改输入文件
无法完成。

## 1. 复制一个可运行的输入文件

不要直接覆盖仓库中的示范输入。建议为新算例建立独立目录：

```bash
mkdir -p examples/my_case
cp examples/h2o2_2d_demo/inputs examples/my_case/inputs
```

然后只编辑 `examples/my_case/inputs`。输出写入 `runs/`，不会和源码混在一起。

## 2. 设置计算域和网格

计算域使用 SI 单位：

```text
geometry.prob_lo = 0.0 0.0
geometry.prob_hi = Lx Ly
amr.n_cell = Nx Ny
```

基础网格间距必须满足：

```text
dx = Lx / Nx
dy = Ly / Ny
```

如果 `amr.max_level = L`，且每层细化比均为 2，则最细网格间距为：

```text
dx_finest = dx / 2^L
dy_finest = dy / 2^L
```

仓库的快速示例使用 `Lx=0.12 m`、`Nx=600`、`Ly=0.01 m`、`Ny=50`，基础间距
为 `200 μm`；三层 2:1 加密后最细间距约为 `25 μm`。建立新算例时，应按照
爆震诱导区、反应区以及目标胞格尺度选择最细网格，而不是只依据基础网格。

`amr.blocking_factor` 必须是 2 的幂。`amr.n_cell`、`amr.max_grid_size` 和
AMR 网格布局必须满足 AMReX 的整除/对齐要求。初次建立算例时可以保留：

```text
amr.blocking_factor = 2
amr.max_grid_size = 64
amr.n_error_buf = 3
amr.do_reflux = 1
```

## 3. 设置点火区和未燃气体

求解器根据单元中心坐标判断 `x <= det.x_igniter` 还是 `x > det.x_igniter`：

```text
x <= det.x_igniter:
    T = det.T_igniter
    p = det.p_igniter

x > det.x_igniter:
    T = det.T_fuel
    p = det.p_fuel
```

示例输入为：

```text
det.x_igniter = 0.005
det.T_igniter = 1800.0
det.p_igniter = 151987.5
det.T_fuel = 300.0
det.p_fuel = 8000.0
```

`det.x_igniter` 必须严格位于 `geometry.prob_lo[0]` 和
`geometry.prob_hi[0]` 之间。速度初始值固定为 0。密度和总能量由求解器根据
理想气体状态方程、NASA7 热力学数据和温度计算，不需要在输入文件中给出。

点火区的温度为：

```text
T(x,y) = T_igniter
         * [1 + perturb_amp * sin(2*pi*perturb_waves*y/Ly)]
```

温度扰动只施加在点火区。`det.perturb_amp = 0` 表示无横向扰动。

## 4. 设置 H2/O2 初始组分

当前输入接口只有两个组分参数：

```text
det.Y_H2 = 0.1119067440
det.Y_O2 = 0.8880932560
```

其余 7 个组分在初始时刻均为 0。这两个质量分数必须满足：

```text
det.Y_H2 + det.Y_O2 = 1
```

当前值对应化学计量的 `2H2 + O2` 混合物。修改组分后，求解器不会自动归一化，
因此应自行保证质量分数和为 1。点火区与未燃区使用相同的初始组分，当前版本
没有单独的“点火区组分”输入参数。

## 5. 设置停止时间和输出

最小时间控制为：

```text
stop_time = 1.0e-5
max_step = 100000
```

计算会在物理时间达到 `stop_time` 或步数达到 `max_step` 时停止，以先到者为准。

常用的按物理时间输出设置为：

```text
amr.plot_int = 0
amr.chk_int = 0
amr.plot_time_int = 1.0e-6
amr.chk_time_int = 5.0e-6
```

`plot_time_int` 控制场数据输出，`chk_time_int` 控制用于续算的 checkpoint。
长时间计算应保留 checkpoint。若只想按固定步数输出，可设置 `plot_int` 和
`chk_int`，并将对应的 `*_time_int` 设为 0。

checkpoint 带有求解器版本。修改物理模型、化学积分或数值通量后，应从
`t=0` 重算；不要直接把旧 checkpoint 接到新代码上。旧格式 checkpoint
默认会拒绝，只有确认兼容并显式设置 `det.allow_legacy_restart=1` 才能读取。

## 6. 设置 AMR 加密

当前求解器会在以下位置加密：

- 压力、温度以及 H2/H2O 梯度较大的单元；
- 化学温度活动度超过 `det.chem_activity` 相对阈值的单元；
- 压力波前后 `det.front_pre` 和 `det.front_post` 距离内的波前窗口；
- `det.refine_zone_until` 之前、`det.refine_zone_lo/hi` 之间的点火区窗口。

一个可运行的起点为：

```text
det.grad_p = 0.15
det.grad_T = 0.10
det.grad_Y = 0.02
det.chem_activity = 0.05
det.front_pre = 1.5e-3
det.front_post = 2.5e-3
det.front_regrid_distance = 5.0e-4
det.refine_zone_lo = 0.004
det.refine_zone_hi = 0.006
det.refine_zone_until = 3.0e-6
```

数值越小通常会增加被标记的单元和显存/计算量。修改后应检查 AMR 层级图，
确认细网格覆盖波面、反应区以及目标物理结构，而不只是关注总加密单元数。

## 7. 检查边界条件和化学设置

当前二维版本固定采用：

- 左边界：反射壁面；
- 右边界：一阶外推；
- 上下边界：周期边界。

因此 `geometry.is_periodic` 应保持为：

```text
geometry.is_periodic = 0 1
```

如需改变边界类型，必须修改 `src/Detonation.cpp` 中 `Detonation` 构造函数里的
`bcs` 设置。正式物理计算应保持：

```text
det.chem_stride = 1
det.chem_Tmin = 650.0
```

`det.chem_stride` 大于 1 会降低化学积分频率，只适合快速调试。

## 8. 运行新算例

在仓库根目录运行：

```bash
DETONATION_RUN_DIR=runs/my_case \
  ./scripts/run_case.sh examples/my_case/inputs
```

另开一个终端查看实时状态：

```bash
python3 scripts/watch_status.py runs/my_case/status.json
```

临时覆盖输入参数时，把 AMReX 参数追加在输入文件后面：

```bash
DETONATION_RUN_DIR=runs/my_case_debug \
  ./scripts/run_case.sh examples/my_case/inputs \
  max_step=100 stop_time=1.0e-6 amr.plot_time_int=5.0e-7
```

## 9. 检查新算例是否正常

启动后至少检查以下内容：

1. 日志中 GPU 名称和 AMReX 版本是否正确；
2. `status.json` 中的 `step`、`t_us`、`dt_s` 是否持续更新；
3. `mass` 和 `energy` 是否出现异常跳变；
4. `finest_level` 是否达到预期，且加密区域是否跟随波面；
5. 输出目录中是否出现初始时刻和后续时刻的场数据；
6. 用 ParaView 打开波形，检查压力、温度、组分是否相互耦合；
7. 改变网格分辨率后检查关键结果是否发生不合理变化。

参数定义全集见 [INPUT_REFERENCE.md](INPUT_REFERENCE.md)，仓库中每个文件的
职责见 [FILE_OVERVIEW.md](FILE_OVERVIEW.md)。
