# AMReX 2-D Detonation Solver

一个基于 AMReX、CUDA 和 MPI 的二维可压缩反应流/爆震求解器。它包含 9 组分、
21 反应 H2/O2 机理、NASA7 热力学、GPU HLLC、AMR、粗/细界面 reflux、
checkpoint/restart，以及 ParaView/Tecplot 数据输出。

> 这是研究代码，不包含燃烧模型验证、数值方法验证或任何形式的安全保证。

## 功能

- 守恒二维 Euler 方程：`rho, rho*u, rho*v, rho*E, rho*Y_k`
- GPU MUSCL-minmod 重构和二维 HLLC
- SSP-RK2 时间推进
- 9 组分、21 反应化学，BDF1 隐式化学与自适应子步
- AMReX `AmrCore` 动态 AMR
- `FluxRegister` 粗/细界面通量修正
- 化学活性、压力/温度/组分梯度与波前窗口联合加密
- checkpoint/restart
- 实时 `status.json` / `status.txt`
- 按物理时间命名的 AMReX 场
- ParaView `.pvd/.vtm/.vtu` 和 Tecplot `.dat`

## 目录

```text
src/                    求解器、物理、化学和格式导出
scripts/build.sh        通用构建脚本
scripts/run_case.sh     通用算例启动脚本
scripts/watch_status.py 无绘图实时状态
tools/                  场提取、绘图和 ParaView 启动工具
examples/h2o2_2d_demo/  快速和长时二维爆震算例
docs/                   输入参数参考
```

## 系统要求

- Linux 或 WSL2；建议 Ubuntu 20.04 或更新版本
- CMake 3.24 或更新版本
- MPI（OpenMPI 或 MPICH）
- CUDA Toolkit
- 二维 CUDA AMReX，建议 AMReX 26.10
- Python 3，仅绘图工具需要 NumPy 和 Matplotlib

当前版本已在 AMReX 26.10、CUDA 12.9 和 Ubuntu 20.04 WSL2 上验证。
求解器只支持二维 CUDA 构建。CPU 构建分支未包含在发布版中。

## 安装 AMReX

如果还没有二维 CUDA AMReX，可从 AMReX 源码安装：

```bash
export AMREX_SOURCE=$HOME/amrex
export CUDA_HOME=/usr/local/cuda
export AMREX_PREFIX=$HOME/amrex-install-2d-cuda

cmake -S "$AMREX_SOURCE" -B "$HOME/amrex-build-2d-cuda" \
  -DCMAKE_BUILD_TYPE=Release \
  -DAMReX_SPACEDIM=2 \
  -DAMReX_GPU_BACKEND=CUDA \
  -DCMAKE_CUDA_ARCHITECTURES=89 \
  -DCMAKE_CUDA_COMPILER="$CUDA_HOME/bin/nvcc" \
  -DAMReX_MPI=ON \
  -DAMReX_OMP=OFF \
  -DAMReX_EB=OFF \
  -DAMReX_PARTICLES=OFF \
  -DCMAKE_INSTALL_PREFIX="$AMREX_PREFIX"

cmake --build "$HOME/amrex-build-2d-cuda" --parallel 4
cmake --install "$HOME/amrex-build-2d-cuda"
```

`CMAKE_CUDA_ARCHITECTURES` 要与本机 GPU 匹配。RTX 4070/4060 使用 `89`，
RTX 30 系列通常使用 `86`，RTX 20 系列通常使用 `75`，A100 使用 `80`。

## 下载源码

可以直接在 GitHub 页面选择 **Code → Download ZIP** 并解压，也可以在终端使用：

```bash
git clone https://github.com/buaa-awy/detonation-amrex.git detonation-amrex
cd detonation-amrex
```

## 构建

进入源码目录后设置 AMReX、CUDA 和 GPU 架构，然后运行构建脚本：

```bash
export AMREX_PREFIX=$HOME/amrex-install-2d-cuda
export CUDA_HOME=/usr/local/cuda
export CUDA_ARCH=89

./scripts/build.sh
./build/chemtest
```

如果 AMReX 或 CUDA 位于其他位置，可通过
`AMREX_PREFIX`、`CUDA_HOME`、`CUDA_ARCH` 指定。`CUDA_HOME` 也可以不设置，
脚本会自动检查 `/usr/local/cuda` 和 `/usr/local/cuda-12.9`。默认并行编译数为
`2`；
内存较小时可用 `BUILD_JOBS=1`。

绘图依赖：

```bash
python3 -m pip install -r requirements.txt
```

## 快速开始

运行二维爆震示范算例：

```bash
cd detonation-amrex
DETONATION_RUN_DIR=runs/quick \
  ./scripts/run_case.sh examples/h2o2_2d_demo/inputs
```

另开终端查看实时进度：

```bash
python3 scripts/watch_status.py runs/quick/status.json
```

状态包括模拟时间、完成百分比、步号、`dt`、AMR 层、活动单元、压力锋、
反应锋、二者间距、墙钟时间和 ETA。

## 初始条件

示范算例使用两段初值：

```text
x <= det.x_igniter:
    T = det.T_igniter
    p = det.p_igniter

x > det.x_igniter:
    T = det.T_fuel
    p = det.p_fuel
```

初始成分为：

```text
Y_H2 = 0.1119067440
Y_O2 = 0.8880932560
其它组分 = 0
```

点火区可叠加 `det.perturb_amp` 和 `det.perturb_waves` 控制的横向正弦温度
扰动。边界为左端反射壁、右端外推、上下周期。

修改初始条件时，至少检查：

1. `geometry.prob_lo/hi` 和 `amr.n_cell` 是否与目标尺寸一致。
2. `det.x_igniter` 是否位于计算域内。
3. 燃料与氧化剂质量分数之和是否为 1。
4. 化学温度阈值 `det.chem_Tmin` 是否覆盖实际反应温区。
5. `det.chem_stride` 是否保持为 `1`。

完整的输入参数说明见 [docs/INPUT_REFERENCE.md](docs/INPUT_REFERENCE.md)。

## 长时算例

长时二维算例位于 `examples/h2o2_2d_demo/inputs-long`，计算域为
`0.30 m x 0.02 m`，三级 2:1 AMR，最细约 25 μm，目标时间为 75 μs。

```bash
DETONATION_RUN_DIR=runs/long \
  ./scripts/run_case.sh examples/h2o2_2d_demo/inputs-long
```

单个 8 GB GPU 建议先使用 `NP=1`，并监控显存：

```bash
nvidia-smi
```

长时间计算应使用 checkpoint。续算时把 `amr.restart` 指向绝对路径：

```bash
DETONATION_RUN_DIR=runs/restart \
  ./scripts/run_case.sh examples/h2o2_2d_demo/inputs-long \
  amr.restart=/absolute/path/to/chkNNNNN
```

## 输出

每个算例目录包含：

```text
runs/<case>/
  status.json
  status.txt
  run.log
  chkNNNNN/
  fields/
    catalog.json
    output_index.csv
    t<seconds>/
    paraview/
      fields.series
      time_series.pvd
      t<seconds>.vtm
      t<seconds>/*.vtu
    tecplot/
      t<seconds>.dat
```

`catalog.json` 描述字段单位与运行参数；`output_index.csv` 汇总时间、步号、
AMR 单元数、波前位置和守恒量。

### ParaView

打开 `fields/paraview/time_series.pvd` 查看全部时间步，或打开
`fields/paraview/t<seconds>.vtm` 查看单时刻。首次打开时在
“Open Data With”中选择标准 XML/VTK reader。

### Tecplot

打开 `fields/tecplot/t<seconds>.dat`。每个 AMR 网格块对应一个 Tecplot
zone；zone 名称包含层级、块号和 MPI rank。

### 内置绘图

AMReX `fextract` 是二维切片提取所需的 AMReX 工具。若它不在默认位置：

```bash
export FEXTRACT=$AMREX_SOURCE/Tools/Plotfile/fextract.gnu.MPI.ex
```

生成温压密、组分和中心线图：

```bash
python3 tools/plot_amrex_snapshot.py \
  runs/quick/fields/t0.000001000000 --jobs 8
```

然后生成流场汇总和 P/T/O2 剖面：

```bash
python3 tools/plot_final_flow.py \
  runs/quick/fields/t0.000001000000
```

AMR 布局图：

```bash
python3 tools/plot_amr_levels.py \
  runs/quick/fields/t0.000001000000 --xlim 0.03
```

## 数值与性能说明

- 正式物理计算必须使用 `det.chem_stride=1`。
- HLLC 使用左右单元声速和 Roe 加权波速。
- 化学活性与波前窗口共同控制局部二维 AMR。
- 8 GB GPU 上 AMReX arena 可能预留数 GB 显存；进程退出后会释放。
- 长算例可能需要数小时，建议先在快速算例上确认环境。

## 许可

本项目以 MIT License 发布，见 [LICENSE](LICENSE)。
