# 输入参数参考

除特别说明外，所有有量纲参数均使用 SI 单位。参数名区分大小写。

## 时间和几何

| 参数 | 含义 |
|---|---|
| `stop_time` | 最终物理时间，单位 s |
| `max_step` | 最大水动力步数 |
| `geometry.prob_lo/hi` | 二维计算域下界/上界，单位 m |
| `geometry.is_periodic` | x、y 方向是否采用周期边界；`0` 为否，`1` 为是 |
| `geometry.coord_sys` | 坐标系类型；二维笛卡尔坐标使用 `0` |

## AMR 和输出

| 参数 | 含义 |
|---|---|
| `amr.n_cell` | 基础网格在 x、y 方向的单元数 |
| `amr.max_level` | 最大加密层号；`0` 表示不加密 |
| `amr.ref_ratio` | 各加密层相对上一层的细化比 |
| `amr.regrid_int` | 两次重网格之间允许的最大步数 |
| `amr.blocking_factor` | AMReX 网格块的最小对齐因子，必须是 2 的幂 |
| `amr.max_grid_size` | 单个 AMReX 网格块的最大单元数 |
| `amr.n_error_buf` | 标记加密区域向外扩充的单元数 |
| `amr.do_reflux` | 是否执行粗/细网格通量修正；正式计算应设为 `1` |
| `amr.plot_file` | 场数据文件名前缀，可包含目录 |
| `amr.chk_file` | checkpoint 文件名前缀 |
| `amr.restart` | 从指定 checkpoint 续算；不设置时从初始条件开始 |
| `amr.plot_vtk` | 是否写出 ParaView VTM/VTU/PVD 文件 |
| `amr.plot_tecplot` | 是否写出 Tecplot ASCII 文件 |
| `amr.plot_export_precision` | 文本格式输出精度 |
| `amr.plot_int` | 按固定步数输出场数据的间隔；`0` 表示关闭 |
| `amr.chk_int` | 按固定步数输出 checkpoint 的间隔；`0` 表示关闭 |
| `amr.plot_time_int` | 场数据输出间隔，单位 s |
| `amr.chk_time_int` | checkpoint 输出间隔，单位 s |

## 爆震模型和初始条件

| 参数 | 含义 |
|---|---|
| `det.cfl` | 水动力 CFL 数 |
| `det.change_max` | 相邻时间步允许的最大增长倍率 |
| `det.x_igniter` | 初始高温点火区右边界，单位 m |
| `det.T_igniter` | 点火区温度，单位 K |
| `det.p_igniter` | 点火区压力，单位 Pa |
| `det.T_fuel` | 未燃混合气温度，单位 K |
| `det.p_fuel` | 未燃混合气压力，单位 Pa |
| `det.perturb_amp` | 横向正弦温度扰动的相对幅值 |
| `det.perturb_waves` | 横向温度扰动的波数 |
| `det.Y_H2`, `det.Y_O2` | 兼容旧输入：点火区和未燃区共同的 H2/O2 质量分数 |
| `det.Y_igniter` | 点火区 9 组分质量分数数组；设置后覆盖 `Y_H2`/`Y_O2` |
| `det.Y_fuel` | 未燃区 9 组分质量分数数组；设置后覆盖 `Y_H2`/`Y_O2` |
| `det.u_igniter`, `det.u_fuel` | 两段初始流向速度，单位 m/s |
| `amr.left_wall` | 左边是否为固壁；0 表示透射边界 |
| `amr.right_wall` | 右边是否为固壁；0 表示透射边界 |
| `amr.slab_refine` | 周期展向是否强制按列使用统一 AMR 层级 |

未指定 `Y_igniter`/`Y_fuel` 时，两段仍使用 `Y_H2 + Y_O2 = 1` 的旧行为。
指定组分数组后可以初始化 CJ 平衡产物等一般状态。左右边界默认分别为固壁和
透射；上下边界为周期。

## 自适应加密判据

| 参数 | 含义 |
|---|---|
| `det.grad_p` | 压力相对梯度加密阈值 |
| `det.grad_T` | 温度相对梯度加密阈值 |
| `det.grad_Y` | H2/H2O 相对梯度加密阈值 |
| `det.chem_activity` | 化学温度活动度相对于峰值的加密阈值 |
| `det.front_pre` | 压力波前前方细网格覆盖距离，单位 m |
| `det.front_post` | 压力波前后方细网格覆盖距离，单位 m |
| `det.front_regrid_distance` | 波前移动触发重网格的附加距离，单位 m |
| `det.refine_zone_lo/hi` | 点火初期临时加密区间的下界/上界 |
| `det.refine_zone_until` | 临时点火加密区的结束时间，单位 s |

## 化学和状态输出

| 参数 | 含义 |
|---|---|
| `det.chem_Tmin` | 低于该温度时跳过化学积分，单位 K |
| `det.chem_hmax` | 化学积分允许的最大子步长，单位 s |
| `det.chem_accum_tol` | 累计多少秒的化学时间后至少尝试一次化学推进 |
| `det.chem_drift_T` | 累计等待期间允许的温度变化，超过后触发化学推进，单位 K |
| `det.chem_tol_y` | 化学子步的组分相对变化容差 |
| `det.chem_tol_T` | 化学子步的温度相对变化容差 |
| `det.chem_stride` | 每多少个水动力步调用一次分块化学；正式计算必须使用 `1` |
| `det.progress_every` | 控制台进度输出间隔，单位：水动力步 |
| `det.status_every` | `status.json/status.txt` 更新间隔，单位：水动力步 |
| `det.status_file` | JSON 状态文件名 |
| `det.status_text_file` | 纯文本状态文件名 |
| `det.allow_legacy_restart` | 默认 `0`；拒绝不带版本元数据的旧 checkpoint。仅在确认兼容时设为 `1`。 |

完整的建算例步骤和示例见
[`NEW_CASE_GUIDE.md`](NEW_CASE_GUIDE.md)。
