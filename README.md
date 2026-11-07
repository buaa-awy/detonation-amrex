# AMReX 2-D Detonation Solver

一个基于 AMReX、CUDA 和 MPI 的二维可压缩反应流/爆震求解器。它包含 9 组分、
21 反应 H2/O2 机理、NASA7 热力学、GPU HLLC、AMR、粗/细界面 reflux、
checkpoint/restart，以及 ParaView/Tecplot 数据输出。

> 当前稳定版本：**v0.1.4**（2026-11-07）。版本变化见
> [CHANGELOG.md](CHANGELOG.md)。

> 这是研究代码。当前版本完成了反应 Euler、HLLC、化学、AMR、GPU 和 I/O 的
> 核心传播验证；不包含网格无关性、实验标定或任何形式的安全保证。

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
docs/                   参数参考、新建算例教程和完整文件说明
VERSION                 当前发布版本号
CHANGELOG.md            版本变更记录
```

## 系统要求

- Linux 或 WSL2；建议 Ubuntu 20.04 或更新版本
- CMake 3.24 或更新版本
- MPI（OpenMPI 或 MPICH）
- CUDA Toolkit
- 二维 CUDA AMReX，固定使用 AMReX 26.10
- Python 3，仅绘图工具需要 NumPy 和 Matplotlib

当前版本已在 AMReX 26.10、CUDA 12.9 和 Ubuntu 20.04 WSL2 上验证。
求解器只支持二维 CUDA 构建。CPU 构建分支未包含在发布版中。

完整的系统依赖、固定版本、GPU 架构表和干净环境步骤见
[docs/ENVIRONMENT.md](docs/ENVIRONMENT.md)。首次复现建议依次执行：

```bash
export CUDA_ARCH=89
./scripts/install_amrex.sh
export AMREX_PREFIX=$HOME/amrex-install-2d-cuda
export FEXTRACT=$AMREX_PREFIX/bin/amrex_fextract
./scripts/check_environment.sh
./scripts/build.sh
./build/chemtest
./scripts/run_smoke_test.sh
```

`CUDA_ARCH` 必须与 GPU 匹配：RTX 20 系列通常为 `75`，RTX 30 系列通常为
`86`，RTX 40 系列为 `89`，A100 为 `80`。

## 安装 AMReX

推荐直接使用脚本安装固定的 AMReX 26.10。脚本会从固定标签克隆源码，开启
AMReX plotfile 工具，并把 `amrex_fextract` 安装到 `AMREX_PREFIX/bin`：

```bash
export CUDA_ARCH=89
./scripts/install_amrex.sh
```

脚本默认使用 `$HOME/amrex`、`$HOME/amrex-build-2d-cuda` 和
`$HOME/amrex-install-2d-cuda`，可通过 `AMREX_SOURCE`、`AMREX_BUILD`、
`AMREX_PREFIX` 修改。

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

与当前验证环境一致的 Python 3.13 版本可使用
`requirements-lock-py313.txt`。系统依赖和版本详情见
[docs/ENVIRONMENT.md](docs/ENVIRONMENT.md)。

## 可复现回归测试

在长时间计算前先运行 100 步确定性测试：

```bash
./scripts/run_smoke_test.sh
```

测试检查最大 AMR 层、最终状态有限性、总质量和总能量。参考值和容差见
[docs/EXPECTED_RESULTS.md](docs/EXPECTED_RESULTS.md)。

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

## 新建一个算例

不要直接修改公共示范输入。建议复制一份输入文件：

```bash
mkdir -p examples/my_case
cp examples/h2o2_2d_demo/inputs examples/my_case/inputs
```

然后编辑 `examples/my_case/inputs` 中的四组参数：

1. `geometry.prob_lo/hi`、`amr.n_cell` 和 `amr.max_level`：计算域尺寸、基础
   网格和加密层数。基础网格间距为 `dx=Lx/Nx`、`dy=Ly/Ny`；每增加一层 2:1
   加密，间距减半。
2. `det.x_igniter`、`det.T_igniter`、`det.p_igniter`：左侧点火区边界、温度
   和压力。
3. `det.T_fuel`、`det.p_fuel`、`det.Y_H2`、`det.Y_O2`：右侧未燃混合气的
   温度、压力和质量分数。当前两段区域使用相同组分，且
   `det.Y_H2 + det.Y_O2` 必须等于 `1`。
4. `stop_time`、`amr.plot_time_int`、`amr.chk_time_int`：最终物理时间以及
   场数据、checkpoint 的时间间隔，全部使用秒。

初始速度固定为 0。密度和总能量由求解器根据压力、温度和 NASA7 热力学数据
计算。`det.perturb_amp` 和 `det.perturb_waves` 只在点火区叠加横向正弦温度
扰动。当前边界固定为左侧反射壁、右侧外推、上下周期。

启动新算例：

```bash
DETONATION_RUN_DIR=runs/my_case \
  ./scripts/run_case.sh examples/my_case/inputs
```

完整的尺寸选择、组分约束、AMR 设置、运行检查和常见错误说明见
[docs/NEW_CASE_GUIDE.md](docs/NEW_CASE_GUIDE.md)。

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
`0.30 m x 0.02 m`，三级 2:1 AMR，最细网格为 `25 μm x 25 μm`，与原 Python
参考算例一致，目标时间为 75 μs。

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

checkpoint 头部写入求解器版本。代码更新后，旧 checkpoint 默认会被拒绝，
避免把旧数值状态混入新求解器。只有在确认物理兼容时，才显式设置
`det.allow_legacy_restart=1` 读取旧格式 checkpoint。

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

`install_amrex.sh` 已构建 AMReX `fextract`。若使用已有 AMReX：

```bash
export FEXTRACT=$AMREX_PREFIX/bin/amrex_fextract
```

输出时间通常略晚于输出间隔，因此不要假设目录名一定是精确的
`t0.000001000000`。选择最新场目录：

```bash
PLOT=$(find runs/quick/fields -maxdepth 1 -type d -name 't*' | sort | tail -1)
echo "$PLOT"
python3 tools/plot_amrex_snapshot.py \
  "$PLOT" --jobs 8
```

然后生成流场汇总和 P/T/O2 剖面：

```bash
python3 tools/plot_final_flow.py \
  "$PLOT"
```

AMR 布局图：

```bash
python3 tools/plot_amr_levels.py \
  "$PLOT" --xlim 0.03
```

大 AMR 算例也可以用 ParaView 的 `pvpython` 直接读取 VTM，不依赖
`fextract`：

```bash
VTM=$(find runs/long/fields/paraview -maxdepth 1 -name '*.vtm' | sort -V | tail -1)
PLOT=$(find runs/long/fields -maxdepth 1 -type d -name 't*' | sort -V | tail -1)
pvpython tools/extract_vtm_fields.py \
  "$VTM" analysis/fields.npz --nx 3001 --ny 401
python3 tools/plot_final_flow.py \
  "$PLOT" --data analysis/fields.npz --outdir analysis/figures
```

## 文件说明

仓库的主要目录如下：

```text
README.md                项目入口和快速使用说明
CMakeLists.txt           CUDA/AMReX 构建配置
src/                     求解器、物理模型、化学机理和输出源码
scripts/                 构建、运行和实时监视脚本
examples/                可直接运行的二维爆震输入
docs/                    参数、新建算例和文件说明
mechanism/               化学机理来源、引用和确定性转换说明
tools/                   AMReX 提取、绘图的辅助工具
requirements.txt          Python 绘图依赖
requirements-lock-py313.txt  验证过的 Python 3.13 绘图依赖
LICENSE                  MIT 许可证
```

其中 `src/Detonation.cpp` 是 AMR 时间推进核心，`src/Physics.H` 包含热力学和
化学积分，`src/MechanismData.H` 保存 9 组分/21 反应机理，
`src/FieldExport.cpp` 负责 ParaView/Tecplot 输出，`src/ChemTest.cpp` 是化学
自检程序。

发布包中每个文件的职责、一次运行产生的 `status/chk/fields/paraview/tecplot`
文件说明见 [docs/FILE_OVERVIEW.md](docs/FILE_OVERVIEW.md)。

## 平面 CJ 爆震验证

`examples/h2o2_2d_demo/inputs-cj-validation` 使用 Cantera 平衡 CJ 状态初始化
平面化学计量 H2/O2 爆轰，不依赖经验 DDT 初值。实测 5–30 μs 波速为
`2830.03 m/s`，相对 Cantera 参考值 `2838.5 m/s` 的误差为 `-0.30%`。
详细参考状态、运行命令和 30 μs 剖面见
[docs/CJ_VALIDATION.md](docs/CJ_VALIDATION.md)。

## 数值与性能说明

- 正式物理计算必须使用 `det.chem_stride=1`。
- HLLC 使用左右单元声速和 Roe 加权波速。
- 化学活性与波前窗口共同控制局部二维 AMR。
- 8 GB GPU 上 AMReX arena 可能预留数 GB 显存；进程退出后会释放。
- 长算例可能需要数小时，建议先在快速算例上确认环境。

## 许可

本项目以 MIT License 发布，见 [LICENSE](LICENSE)。
