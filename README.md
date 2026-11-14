# Cl / Si 刻蚀仿真

逐个追踪入射粒子，通过表面捕获、反射和材料更新推进刻蚀。代码分为粒子源、网格、曲线拟合、轨迹追踪和表面反应组件。

当前保留二维局部曲线原型。论文采用三维元胞和局部二次曲面；三维部分尚未实现。

## 结构

```text
src/
  source.hpp          入射位置、角度和能量分布
  grid.hpp/.cpp       材料网格和表面位点库存
  surface.hpp/.cpp    局部最小二乘二次曲线
  transport.hpp/.cpp  网格遍历、曲线求交和反射追踪
  reaction.hpp/.cpp   产额测试模型、单 Cl 逐步氯化模型
  simulation.hpp/.cpp 逐粒子推进和统计
  vector.hpp          向量运算
  main.cpp            参数读取和组件装配
utils/
  build.py            构建和运行测试
  results.py          结果读取与库存核对
  plot.py             原始网格和局部曲线绘图
tests/                单元测试和几何追踪回归测试
docs/etching/         参考文献
outputs/              仿真输出
```

## 角度函数和表面反应

令 `d` 为粒子的单位运动方向，`n` 为固体指向气体的局部曲线法线：

$$
\mu=\cos\theta=-\mathbf d\cdot\mathbf n.
$$

于骁论文图 4-5、4-6 的产额测试函数为：

$$
F_A(\theta)=\cos\theta=\mu,
\qquad F_B(\theta)=(1+4\sin^2\theta)\cos\theta=(5-4\mu^2)\mu.
$$

这两个函数是产额，`F_B` 可以超过 1。

### `yield`：产额测试

按照于骁第 4.1.1 节的先后顺序，先以 `cos(theta)` 判断粒子是否渗入，未渗入则反射；渗入后直接以 `F(theta)` 计算刻蚀产额并扣除硅库存。代码对非整数产额采用无偏随机取整，单个网格的去除量受剩余库存限制。两步分别描述渗入和渗入后的产额，不使用 `F/cos(theta)` 补偿。

`--angular cosine`、`--angular enhanced` 分别选择两条产额函数。

### `chlorination`：单 Cl 逐步氯化

保留用户提出的简化反应链：

$$
\mathrm{Si}\rightarrow\mathrm{SiCl}\rightarrow\mathrm{SiCl_2}
\rightarrow\mathrm{SiCl_3}\rightarrow\mathrm{SiCl_4(g)}.
$$

每个被捕获的入射粒子只提供一个 Cl。从碰撞网格的硅位点中随机选择一个，将其氯化程度推进一级；形成 SiCl₄ 后脱附，硅库存减一。

这个模型不使用产额函数倍增 Cl，也不自动设置为四倍入射剂量。它是上述简化反应链的实现，不是论文第 3 章完整混合层反应模型。`--angular` 只应用于 `yield` 模型。

## 网格与存储

| 项目 | 设置 | 依据 |
|---|---|---|
| 网格边长 | 2.5 nm | 于骁，第 4.1 节 |
| 每格初始粒子数 | 30 | 于骁，第 4.1 节 |
| 元胞移除条件 | 剩余数量小于 5 | 于骁，第 4.1.1 节 |
| 局部拟合邻域 | 目标格周围 ±3 格 | 原文局部邻域设置的二维实现 |
| 横向边界 | 周期边界 | 原文结构示意与粒子追踪设置 |
| 反射后能量 | 保持不变 | 于骁，第 4.3.1 节 |

`x` 为横向坐标，`z` 向下。初始硅表面平坦，上方是固定掩膜及其开口。整体网格数、掩膜尺寸、开口尺寸和入射剂量由调用者明确给出。

网格使用连续数组，索引为 `z * nx + x`：

- `material`：`uint8`；0 为气体、1 为硅、2 为掩膜。
- `sites`：每格四个 `uint16`，记录 Si、SiCl、SiCl₂、SiCl₃ 的数量。
- `boundary` 与 `surface_cells`：记录暴露于气体的固体单元。
- `primary_hits`、`secondary_hits`：`uint32`，累计首次撞击和反射后的撞击。

移除元胞时，剩余的 Si、Cl 分别计入 `pruned_si`、`pruned_cl`。这些是删格残量，不计作化学脱附产物。库存关系为：

$$
N_{\mathrm{Si,initial}}=N_{\mathrm{Si,remaining}}+N_{\mathrm{Si,reacted}}+N_{\mathrm{Si,pruned}}.
$$

逐步氯化模型另外满足：

$$
N_{\mathrm{Cl,captured}}=N_{\mathrm{Cl,surface}}+4N_{\mathrm{SiCl_4}}+N_{\mathrm{Cl,pruned}}.
$$

## 入射角与能量

粒子在整个源平面均匀发射。离子发射角采用高斯分布，中性粒子角分布可采用二维余弦分布；发射角相对于整体向下方向定义，碰撞角相对于局部法线定义。

- `--angle gaussian`：截断到 ±90° 的正态分布；`sigma = FWHM / 2.35482`。
- `--angle cosine`：二维密度 `cos(alpha)/2`。
- `--angle normal`：法向单束，用于几何测试。
- `--energy constant`：指定单能。
- `--energy gaussian`：正能量范围内的高斯抽样。

FWHM、能量及能量宽度由调用者给出。当前产额测试和单 Cl 反应链没有附加能量阈值、饱和函数或经验能量权重，反射保留入射能量。

## 局部拟合与碰撞

先沿轨迹定位粗网格上的表面元胞，再收集其附近的表面元胞中心。几何取点包含不同材料的气固边界，在局部旋转坐标系中作二次最小二乘拟合：

$$
v=f(u)=au^2+bu+c.
$$

**碰撞位置和反射法线都由该拟合曲线决定。** 在所选局部拟合邻域内求轨迹与曲线的前向交点，由曲线导数计算法线，按镜面反射公式更新方向。

反应元胞按刘婉第 3.5 节的规则定位：若交点所在格是表面元胞，直接使用该格；若交点位于气体格或内部非表面格，沿交点法线的两个方向查找该局部邻域内最近的表面元胞。没有符合条件的元胞时记录映射失败，不直接使用拟合片所属格。

元胞被移除后，重算附近受影响的拟合；只有表面氯化程度改变时复用原拟合。`rank_deficient_fits` 统计局部二次拟合的秩不足次数。

## 构建和测试

需要 C++17 编译器、Python、NumPy 和 Matplotlib。

```bash
python -m pip install -r requirements.txt
python -m utils.build
```

追踪和元胞归属的独立回归测试：

```bash
python -m utils.build --geometry
```

运行时明确提供几何和源参数，以下是命令格式，其中大写占位符替换为输入值：

```text
build/etch_sim.exe --model yield --angular cosine --nx NX --nz NZ --surface-z Z0 --opening W --mask-thickness H --fwhm FWHM_DEG --energy-ev E_EV --particles N --output OUTPUT_DIR
```

`--model chlorination` 使用单 Cl 反应链；`--reflect 0` 用于检查二次追踪。默认不设反射次数上限；`--max-bounces` 只用于数值诊断，并单独统计截断次数。

绘制已有结果：

```text
python -m utils.plot RUN_DIRECTORY --output FIGURE.png
```

## 数值诊断

`summary.json`、`arguments.txt` 和原始数组保存每次运行的输入与统计。当前诊断重点是：

- `interior_exits`：没有找到曲线碰撞却从硅片底部离开的轨迹，属于漏检。
- `normal_mapped_hits`：交点不是表面元胞、经过法线投影定位的次数。
- `mapping_failures`：已找到前向交点，但无法沿法线定位局部表面元胞的候选次数。
- `no_local_intersections`：候选局部拟合与前向轨迹没有有效交点的次数，包含粒子反射后离开原表面的情况，不等同于漏检轨迹数。
- `rank_deficient_fits`：局部拟合矩阵不满秩的次数。
- `pruned_si`、`pruned_cl`：按论文删格条件产生的残量。

`tests/test_geometry.cpp` 单独检查追踪和元胞归属。基础反应与守恒测试通过后，仍需通过这些几何测试。

## 扩展

新增表面模型实现 `Reaction::apply`；新增入射角或能量分布提供 `double(Random&)` 函数，装入 `Source`；产额模型接收独立的 `AngularYield` 函数。源粒子和表面事件使用独立随机数流。

## 参考文献

- [于骁，等离子体刻蚀的动力学模型及三维模拟](docs/etching/等离子体刻蚀的动力学模型及三维模拟_于骁.pdf)，第 4.1、4.3 节。
- [刘婉，等离子刻蚀微观过程的动态蒙特卡洛仿真](docs/etching/等离子刻蚀微观过程的动态蒙特卡洛仿真_刘婉.pdf)，第 3.3～3.5 节。
