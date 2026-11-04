# H2/O2 化学机理来源

求解器使用 9 组分、21 反应记录的 H2/O2 机理。数据来源为 Cantera 3.2.0
官方 `h2o2.yaml`，原始研究为：

> M. Ó Conaire, H. J. Curran, J. M. Simmie, W. J. Pitz, and C. K. Westbrook,
> “A comprehensive modeling study of hydrogen oxidation,”
> *International Journal of Chemical Kinetics*, 36:603–622, 2004.
> DOI: 10.1002/kin.20036

原始 Cantera 文件包含 Ar 和 29 条反应记录。项目做了以下确定性转换：

- 删除 Ar，保留 `H2 H O O2 OH H2O HO2 H2O2 N2` 9 个组分；
- 将同一化学方程的两条重复 Arrhenius 记录合并为一项的两条速率通道；
- 将 `H+O2+O2/H+O2+H2O/H+O2+N2` 和 `2H+H2/2H+H2O` 的增强通道
  转为等效第三体效率；
- 生成 C++ 头文件时把 cm-mol-cal 单位转换为 SI；
- 保留 NASA7 低温和高温热力学多项式。

`src/MechanismData.H` 是当前求解器实际编译的数据，使用前建议核对：

```text
SHA256 42EFA78C5E8B147CB7133193C1E4090E868146670FF4635DD4C87360E3D6737D
```

生成中间 YAML 的校验值：

```text
mech_9sp21r.yaml
SHA256 C9F8CA8498AC61412385E3D432CFE4FF56A8BE46C69E8069A910B1282F40D416
```

21 条反应记录为：

1. `H + O2 <=> O + OH`
2. `O + H2 <=> H + OH`
3. `OH + H2 <=> H + H2O`
4. `2 OH <=> O + H2O`
5. `H + O2 + M <=> HO2 + M`
6. `H + O2 + H2O <=> HO2 + H2O`
7. `H + HO2 <=> 2 OH`
8. `H + HO2 <=> O2 + H2`
9. `H + HO2 <=> O + H2O`
10. `HO2 + O <=> O2 + OH`
11. `OH + HO2 <=> O2 + H2O`
12. `2 HO2 <=> O2 + H2O2`
13. `2 OH (+M) <=> H2O2 (+M)`
14. `OH + H2O2 <=> HO2 + H2O`
15. `H + H2O2 <=> HO2 + H2`
16. `H + H2O2 <=> OH + H2O`
17. `H2O2 + O <=> HO2 + OH`
18. `2 H + M <=> H2 + M`
19. `H + OH + M <=> H2O + M`
20. `O + H + M <=> OH + M`
21. `2 O + M <=> O2 + M`

速率常数形式为：

```text
k = A T^b exp(-Ea/RT)
```

原始机理使用 cm-mol-cal；求解器中已转换为 SI。
