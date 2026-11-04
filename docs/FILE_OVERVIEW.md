# 文件与目录说明

本文说明 GitHub 发布包中每个源码、脚本和配置文件的用途。`build/`、`runs/`
和 `fields/` 是本地生成目录，不会提交到 Git。

## 根目录

| 文件 | 用途 |
|---|---|
| `.gitattributes` | 固定 Shell、Python、C++ 和 Markdown 文件使用 LF 换行，避免从 Windows 克隆后在 WSL/Linux 中出现脚本执行错误。 |
| `.gitignore` | 排除 `build/`、`runs/`、checkpoint、场数据、图片、日志和 Python 缓存。 |
| `CMakeLists.txt` | 配置 C++20/CUDA 构建，查找 AMReX，并生成 `detonation2d` 和 `chemtest` 两个可执行程序。 |
| `README.md` | 项目首页，包含安装、构建、下载、运行、实时监视、输出查看和简单建算例流程。 |
| `LICENSE` | MIT 开源许可证，说明复制、修改和再发布的条件。 |
| `requirements.txt` | Python 绘图工具需要的 `numpy` 和 `matplotlib`。 |
| `requirements-lock-py313.txt` | 已验证 Python 3.13 绘图环境的精确 NumPy/Matplotlib 版本。 |

## 文档

| 文件 | 用途 |
|---|---|
| `docs/INPUT_REFERENCE.md` | 全部输入参数及单位的中文参考表。 |
| `docs/NEW_CASE_GUIDE.md` | 从复制输入文件到运行、检查结果的新建算例完整教程。 |
| `docs/FILE_OVERVIEW.md` | 当前文件，逐个解释仓库文件及运行输出目录。 |
| `docs/ENVIRONMENT.md` | 固定版本、系统依赖、AMReX、CUDA、MPI、Python 和 `fextract` 的安装步骤。 |
| `docs/EXPECTED_RESULTS.md` | 化学自检、100 步回归测试和当前物理验证状态的参考值。 |
| `docs/VALIDATION_70US.md` | 70 μs 长算例的场图、波前速度、守恒量和爆震耦合诊断。 |
| `docs/images/*.png` | 70 μs 算例的温度、压力、AMR、密度梯度、中心线剖面和波前历史图。 |

## 化学机理

| 文件 | 用途 |
|---|---|
| `mechanism/README.md` | 9 组分/21 反应机理的 Cantera 来源、论文引用、转换方法、反应列表和 SHA256。 |

## 示例

| 文件 | 用途 |
|---|---|
| `examples/h2o2_2d_demo/README.md` | H2/O2 二维示范算例的物理设置和运行说明。 |
| `examples/h2o2_2d_demo/inputs-smoke` | 固定 100 步的确定性安装/回归测试输入。 |
| `examples/h2o2_2d_demo/inputs` | 5 μs 快速检查算例，基础网格 600 x 50，三层 AMR。 |
| `examples/h2o2_2d_demo/inputs-long` | 75 μs 长算例，基础网格 1504 x 104，输出 ParaView 和 Tecplot 数据。 |

输入文件不是可执行程序。用户复制其中一个文件并修改参数即可建立新算例。

## 构建和运行脚本

| 文件 | 用途 |
|---|---|
| `scripts/build.sh` | 根据 `AMREX_PREFIX`、`CUDA_HOME`、`CUDA_ARCH` 和 `BUILD_JOBS` 配置 Release 构建。 |
| `scripts/install_amrex.sh` | 从固定标签安装二维 CUDA AMReX 26.10，并构建/安装 `amrex_fextract`。 |
| `scripts/check_environment.sh` | 检查 Git、编译器、CMake、MPI、CUDA、AMReX、GPU 和绘图提取工具。 |
| `scripts/run_case.sh` | 在指定运行目录启动 MPI/GPU 求解器，把终端输出追加到 `run.log`，并允许附加 AMReX 参数覆盖。 |
| `scripts/watch_status.py` | 读取 `status.json`，在终端显示完成百分比、物理时间、步数、时间步长、AMR 层、网格数、波前位置、墙钟时间和 ETA。 |
| `scripts/run_smoke_test.sh` | 创建独立运行目录并执行 100 步确定性回归。 |
| `scripts/check_smoke.py` | 校验回归的步数、AMR 层、完成标记、总质量和总能量。 |

## C++/CUDA 源码

| 文件 | 用途 |
|---|---|
| `src/main.cpp` | 程序入口。初始化 AMReX 和 CUDA，打印 GPU 名称，创建求解器对象并调用初始化与演化流程。 |
| `src/Detonation.H` | `Detonation` 类声明和全部运行参数默认值，包括 AMR 状态、化学状态、波前诊断、输出时间和 restart 路径。 |
| `src/Detonation.cpp` | 核心实现：参数读取、初值、时间推进、GPU MUSCL/HLLC、SSP-RK2、化学子步、AMR 标记、reflux、状态输出、plotfile/checkpoint 读写和 restart。 |
| `src/Physics.H` | 物理模型实现：守恒量/原始量转换、混合物状态、NASA7 热力学、声速、反应速率、化学隐式积分和温度反解。 |
| `src/MechanismData.H` | 9 组分、21 反应 H2/O2 机理数据，包括组分名、相对分子质量、NASA7 系数和反应参数。 |
| `src/FieldExport.H` | ParaView/Tecplot 外部格式导出的函数声明。 |
| `src/FieldExport.cpp` | 把各 AMR 层的 MultiFab 写成 `.vtu/.vtm/.pvd` 和 Tecplot ASCII `.dat`，处理 MPI rank 临时文件和单时刻/时间序列索引。 |
| `src/ChemTest.cpp` | 独立的化学刚性和热力学自检程序。运行 `./build/chemtest` 可在正式算例前检查 GPU 化学求解是否正常。 |

## Python 工具

| 文件 | 用途 |
|---|---|
| `tools/plot_amrex_snapshot.py` | 使用 AMReX `fextract` 从 plotfile 提取二维切片，绘制密度、压力、温度、反应进度和组分图。 |
| `tools/plot_amr_levels.py` | 读取 AMReX `Cell_H` 网格布局，绘制各层加密区域和网格数量。 |
| `tools/plot_final_flow.py` | 从 `plot_amrex_snapshot.py` 或 `extract_vtm_fields.py` 生成的 NPZ 及 AMReX plotfile 生成流场汇总、波前缩放图及中心线剖面图。 |
| `tools/extract_vtm_fields.py` | 使用 ParaView `pvpython` 把多级 VTM 合成为规则二维 NPZ，适合大 AMR 算例且不依赖 `fextract`。 |
| `tools/plot_diagnostics.py` | 从 `output_index.csv` 绘制前导波/反应波位置、速度、间距、守恒漂移和 AMR 单元历史。 |
| `tools/open_paraview_snapshot.py` | 在 ParaView Python 环境中打开指定 `.vtm`，按温度着色并显示。 |

## 一次运行生成的文件

`DETONATION_RUN_DIR` 指定运行目录。以 `runs/my_case` 为例：

| 路径 | 内容 |
|---|---|
| `runs/my_case/run.log` | 终端日志，包含 GPU、AMReX、每条 PROGRESS/CONSERVED 信息和错误信息。 |
| `runs/my_case/status.json` | 实时状态的主文件，供监视或脚本读取。 |
| `runs/my_case/status.txt` | 同一状态的纯文本版本。 |
| `runs/my_case/chkNNNNN/` | 按步号命名的 AMReX checkpoint，可用 `amr.restart=` 续算。 |
| `runs/my_case/fields/t<seconds>/` | AMReX 原生多级 plotfile，每个时刻一个目录。 |
| `runs/my_case/fields/paraview/t<seconds>.vtm` | 单时刻 ParaView 多块数据集入口。 |
| `runs/my_case/fields/paraview/time_series.pvd` | ParaView 时间序列入口，打开后可一次加载全部已输出时刻。 |
| `runs/my_case/fields/tecplot/t<seconds>.dat` | 单时刻 Tecplot ASCII 数据，每个 AMR 网格块是一个 zone。 |
| `runs/my_case/fields/catalog.json` | 输出变量单位、停止时间、AMR 和输出配置。 |
| `runs/my_case/fields/output_index.csv` | 每次输出的时间、步号、总网格数、各层网格数、波前位置、总质量和总能量。 |
| `runs/my_case/fields/figures/` | 使用 `tools/` 脚本后生成的 PNG 图片和中间数据。 |

## 下载包中没有包含的内容

- AMReX 源码或安装目录：由用户通过 `AMREX_PREFIX` 指向自己的二维 CUDA 安装。
- CUDA Toolkit 和 GPU 驱动：由用户在系统中安装。
- `build/`：应在目标机器上重新构建。
- checkpoint 和大体积场数据：使用 `runs/` 保存，不进入 Git。
- CPU 测试分支：当前发布版只支持二维 CUDA 构建。
