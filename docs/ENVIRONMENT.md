# 可复现环境

本页给出从干净 Linux/WSL2 系统建立求解器环境的步骤。求解器本身只要求
Linux、MPI、CUDA 和二维 CUDA AMReX；Python 只用于监视和绘图。

## 已验证环境

当前 `v0.1.2` 开发版本使用以下环境完成干净克隆、Release 编译、化学自检和
100 步回归测试：

| 组件 | 已测试版本 |
|---|---|
| 系统 | Windows 11 + WSL2，Ubuntu 20.04 |
| GPU | NVIDIA GeForce RTX 4070 Laptop GPU，8188 MiB |
| NVIDIA 驱动 | 621.01 |
| CUDA Toolkit | 12.9.86 |
| GPU 架构 | `sm_89`，项目参数 `CUDA_ARCH=89` |
| 主机编译器 | GCC 11.4.0 |
| CMake | 4.4.3；项目最低要求 3.24 |
| MPI | OpenMPI |
| AMReX | 标签 `26.10`，提交 `015c14b65715b277fbc54706ebe7455b196d2254` |
| 绘图验证 | Python 3.13.0、NumPy 2.2.6、Matplotlib 3.11.2 |

`chemtest` 和短算例不要求 Python。实时状态脚本仅使用 Python 标准库。

## 1. 安装基础工具

Ubuntu 20.04 的系统 CMake 版本低于项目要求，因此安装新版 CMake：

```bash
sudo apt update
sudo apt install -y build-essential git openmpi-bin libopenmpi-dev \
  python3 python3-pip python3-venv nvidia-cuda-toolkit
python3 -m pip install --user "cmake==4.4.3"
```

不同系统不一定提供合适的 `nvidia-cuda-toolkit`。实验机器实际使用 NVIDIA
CUDA Toolkit 12.9，`nvcc` 位于 `/usr/local/cuda-12.9/bin/nvcc`。如果没有
预装 CUDA，应按 NVIDIA 针对 WSL2 或原生 Linux 的 CUDA Toolkit 安装说明安装，
并确保：

```bash
export CUDA_HOME=/usr/local/cuda
test -x "$CUDA_HOME/bin/nvcc"
nvidia-smi
```

WSL2 需要 Windows 侧 NVIDIA 驱动支持 CUDA，通常不在 Linux 内另装显示驱动。

## 2. 安装固定版本 AMReX

项目脚本固定使用 AMReX 26.10，并开启 AMReX plotfile 工具，因此会同时生成
`amrex_fextract`：

```bash
export CUDA_ARCH=89       # 按实际 GPU 修改
./scripts/install_amrex.sh
```

默认路径为：

```text
$HOME/amrex                     AMReX 26.10 源码
$HOME/amrex-build-2d-cuda      AMReX 构建目录
$HOME/amrex-install-2d-cuda    AMReX 安装目录
```

脚本不会覆盖一个版本不符或非 Git 的已有源码目录。常用 GPU 架构：

| GPU | `CUDA_ARCH` |
|---|---|
| RTX 20 系列 | `75` |
| RTX 30 系列 | `86` |
| RTX 40 系列 | `89` |
| A100 | `80` |
| H100 | `90` |

## 3. 构建求解器

```bash
export AMREX_PREFIX=$HOME/amrex-install-2d-cuda
export CUDA_HOME=/usr/local/cuda
export CUDA_ARCH=89
export FEXTRACT=$AMREX_PREFIX/bin/amrex_fextract

./scripts/check_environment.sh
./scripts/build.sh
./build/chemtest
```

`check_environment.sh` 不应报告 `git`、`g++`、`cmake`、`mpirun`、`nvcc` 或
AMReX 缺失。

## 4. 运行确定性回归

```bash
./scripts/run_smoke_test.sh
```

成功时会输出类似：

```text
Smoke test PASSED: step=100 t=0.225629 us finest=3
mass=1260.95732991 energy=193540164.212
```

这个测试固定运行 100 步，检查 AMR 达到第 3 层、最终状态有限、总质量与总能量
在容差内。完整参考值见 [EXPECTED_RESULTS.md](EXPECTED_RESULTS.md)。

## 5. Python 绘图环境

可选的精确 Python 绘图环境：

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements-lock-py313.txt
```

`requirements-lock-py313.txt` 与当前验证环境一致。Ubuntu 20.04 默认 Python
通常不是 3.13；可以只安装 `requirements.txt` 中的最低兼容版本。绘图不影响
求解器构建和主要输出。

## 6. ParaView 与 `fextract`

- 打开 `.vtm`、`.pvd` 和 `.vtu`：仅需 ParaView。
- 传统逐行切片图：使用 AMReX `fextract`，由 `install_amrex.sh` 一并构建。
- 大 AMR 算例二维合成：使用 `tools/extract_vtm_fields.py`，直接在 ParaView
  `pvpython` 中读取已有 VTM，不再需要 `fextract`。

运行大 AMR 合成工具：

```bash
VTM=$(find runs/long/fields/paraview -maxdepth 1 -name '*.vtm' | sort -V | tail -1)
pvpython tools/extract_vtm_fields.py \
  "$VTM" analysis/fields.npz --nx 3001 --ny 401
```

然后：

```bash
PLOT=$(find runs/long/fields -maxdepth 1 -type d -name 't*' | sort -V | tail -1)
python3 tools/plot_final_flow.py \
  "$PLOT" --data analysis/fields.npz --outdir analysis/figures
```
