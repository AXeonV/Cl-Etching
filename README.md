# Cl / Si 刻蚀仿真

在二维材料网格上逐个追踪入射粒子，通过表面捕获、反射、反应和元胞更新推进刻蚀。粒子源、几何追踪、表面反应和元胞更新相互独立。

## 结构

```text
src/
  source.hpp          入射位置、角度和能量分布
  grid.hpp/.cpp       材料网格和表面位点库存
  surface.hpp/.cpp    局部最小二乘二次曲线及法向投影
  transport.hpp/.cpp  网格遍历、曲线求交和反射追踪
  reaction.hpp/.cpp   产额模型、单 Cl 逐步氯化模型
  cell_update.hpp/.cpp 元胞删除位置选择与组分转移
  simulation.hpp/.cpp 逐粒子推进和统计
  vector.hpp          向量运算
  main.cpp            参数读取和组件装配
utils/
  build.py            构建和运行测试
  results.py          结果读取与库存核对
  plot.py             材料网格和局部曲线绘图
tests/                单元测试和几何追踪回归测试
outputs/              仿真输出
```

## 碰撞角、捕获与表面反应

令 $\mathbf d$ 为粒子的单位运动方向，$\mathbf n$ 为碰撞点处从固体指向气体的曲线单位法线，$\theta$ 为相对内法线的入射角：

$$
\mu=\cos\theta=-\mathbf d\cdot\mathbf n,\qquad 0\leq\theta\leq\frac{\pi}{2}.
$$

硅表面采用捕获概率 $p(\theta)=\mu$。捕获后执行所选表面反应；未捕获的粒子作镜面反射，继续追踪后续碰撞。固定掩膜表面只发生反射。

### `yield`：产额模型

`--angular cosine` 和 `--angular enhanced` 分别选择于骁论文图 4-5、4-6 中的角度函数：

$$
F_A(\theta)=\cos\theta=\mu,
$$

$$
F_B(\theta)=(1+4\sin^2\theta)\cos\theta=(5-4\mu^2)\mu.
$$

$F(\theta)$ 表示每个入射离子的平均硅原子去除数。捕获概率为 $p(\theta)$ 时，捕获事件的条件期望去除量为：

$$
Y_{\mathrm{capture}}(\theta)=\frac{F(\theta)}{p(\theta)}.
$$

对条件期望进行随机取整：

$$
K=\lfloor Y_{\mathrm{capture}}\rfloor+
\mathrm{Bernoulli}(Y_{\mathrm{capture}}-\lfloor Y_{\mathrm{capture}}\rfloor).
$$

实际去除数为 $\min(K,N_{\mathrm{Si,cell}})$。库存充足时，每次表面撞击的期望反应去除量为 $pY_{\mathrm{capture}}=F$。在 $p=\cos\theta$ 的设置下，两条角度函数对应的捕获后期望去除量分别为 $1$ 和 $1+4\sin^2\theta$。删格残量另行统计。

### `chlorination`：Cl 化反应模型

每个被捕获的粒子提供一个 Cl，在碰撞元胞的全部硅位点中等概率选择一个，将其氯化程度推进一级：

$$
\mathrm{Si}\rightarrow\mathrm{SiCl}\rightarrow\mathrm{SiCl_2}
\rightarrow\mathrm{SiCl_3}\rightarrow\mathrm{SiCl_4(g)}.
$$

形成 SiCl₄ 后脱附，硅库存减一。四种固相状态分别保存数量；该模型采用用户指定的四步反应链，`--angular` 参数用于 `yield` 模型。

## 网格与存储

| 项目 | 设置 |
|---|---|
| 网格边长 | 2.5 nm |
| 每格初始硅原子数 | 30 |
| 元胞更新条件 | 剩余硅位点数小于 5 |
| 局部拟合邻域 | 目标格周围 ±3 格 |
| 元胞删除候选邻域 | 目标格及周围 ±1 格的硅表面元胞 |
| 法线推进与平滑权重比 | $m/n=0.3$ |
| 横向边界 | 周期边界 |

网格尺度、初始库存和更新阈值采用于骁第 4.1 节的设置；元胞删除位置按刘婉第 3.6.2 节的方法选择。程序使用二维曲线和八邻域表示其局部几何。

$x$ 为横向坐标，$z$ 向下。初始硅表面平坦，上方是固定掩膜及其开口。网格数、掩膜尺寸、开口尺寸和入射剂量通过运行参数指定。

网格使用连续数组，索引为 `z * nx + x`：

- `material`：`uint8`；0 为气体、1 为硅、2 为掩膜。
- `sites`：每格四个 `uint16`，记录 Si、SiCl、SiCl₂、SiCl₃ 的数量。
- `boundary` 与 `surface_cells`：记录与气体共边的固体元胞。
- `primary_hits`、`secondary_hits`：`uint32`，累计首次撞击和反射后的撞击。

## 元胞更新与悬挂元胞处理

当反应元胞的剩余库存小于 5 时，在该元胞及其八邻域内的硅表面元胞中选择删除位置。选择同时考虑沿法线推进和减少表面凸起，以降低直接删格产生的悬挂、漂浮结构。

对每个候选元胞计算偏离局部外法线的角度 $\alpha_i$，以及中心到更新前拟合曲线的有符号法向距离 $d_i$，选择评分最小的元胞：

$$
S_i=m\alpha_i+nd_i,\qquad m/n=0.3.
$$

二维实现中，$\alpha_i$ 用弧度，$d_i$ 以元胞边长为单位；$d_i$ 在曲线固体侧为正、气体侧为负，使气体侧凸起优先参与移除。原待删除元胞的 $\alpha_i=0$，其余候选的角度由两格中心连线与外法线计算。法向距离通过求候选中心在二次曲线上的最近点获得。

若选中邻格，将其四种固相库存完整转移到原待删除元胞，再删除选中格。掩膜保持固定。实际删除位置附近的边界和拟合缓存随之更新。

原待删除元胞的阈值残量分别计入 `pruned_si`、`pruned_cl`，组分转移保持 Si 和结合 Cl 数量不变。库存关系为：

$$
N_{\mathrm{Si,initial}}=N_{\mathrm{Si,remaining}}+N_{\mathrm{Si,reacted}}+N_{\mathrm{Si,pruned}}.
$$

逐步氯化模型另外满足：

$$
N_{\mathrm{Cl,captured}}=N_{\mathrm{Cl,surface}}+4N_{\mathrm{SiCl_4}}+N_{\mathrm{Cl,pruned}}.
$$

## 入射角与能量

粒子在整个源平面均匀发射。发射角 $\alpha$ 相对于整体向下方向定义，碰撞角 $\theta$ 相对于局部内法线定义。

- `--angle gaussian`：截断到 ±90° 的正态分布；$\sigma=\mathrm{FWHM}/2.35482$。
- `--angle cosine`：二维余弦分布，密度为 $\cos\alpha/2$。
- `--angle normal`：法向单束。
- `--energy constant`：指定单能。
- `--energy gaussian`：正能量范围内的高斯抽样。

离子采用高斯角分布，中性粒子可采用余弦角分布。FWHM、能量及能量宽度通过运行参数指定。当前两种表面反应模型的反应量由角度或氯化状态决定；粒子携带抽样能量，镜面反射后能量保持不变。

## 局部拟合与碰撞

局部坐标系的 $v$ 轴取邻域固体元胞几何质心指向目标元胞中心的方向，按刘婉第 3.4.2 节旋转坐标。收集附近气固边界元胞的中心，作二次最小二乘拟合：

$$
v=f(u)=au^2+bu+c.
$$

**碰撞位置和反射法线都由该拟合曲线决定。** 轨迹采用覆盖网格角点的 DDA 遍历；每段检索附近拟合片，在局部邻域及采样点的切向覆盖范围内选择最近的前向交点。由曲线导数计算法线，粒子从该交点按镜面方向继续运动：

$$
\mathbf d'=\mathbf d-2(\mathbf d\cdot\mathbf n)\mathbf n.
$$

反应元胞按刘婉第 3.5 节定位：若交点所在格是表面元胞，直接使用该格；否则沿交点法线的两个方向查找局部邻域内最近的表面元胞。

拟合曲线与候选表面元胞列表分别缓存。元胞删除后更新受影响的拟合与候选列表；仅库存或氯化状态改变时复用几何缓存。

## 构建和测试

需要 C++17 编译器、Python、NumPy 和 Matplotlib。

```bash
python -m pip install -r requirements.txt
python -m utils.build
python -m utils.build --geometry
```

运行命令中的大写占位符替换为输入值：

```text
build/etch_sim.exe --model yield --angular cosine --nx NX --nz NZ --surface-z Z0 --opening W --mask-thickness H --fwhm FWHM_DEG --energy-ev E_EV --particles N --output OUTPUT_DIR
```

`--model chlorination` 使用单 Cl 反应链。`--reflect 0` 关闭反射追踪；默认不设反射次数上限，`--max-bounces` 设置诊断用的截断次数。

绘制已有结果：

```text
python -m utils.plot RUN_DIRECTORY --output FIGURE.png
```

## 输出与诊断

`summary.json`、`arguments.txt` 和原始数组保存输入、库存和轨迹统计。

- `deleted_cells`：删除的元胞数。
- `relocated_cells`：选择邻格删除并转移库存的次数。
- `pruned_si`、`pruned_cl`：阈值删格残量。
- `interior_exits`：未找到曲线碰撞却从硅片底部离开的漏检轨迹数。
- `normal_mapped_hits`：交点经过法线投影定位反应元胞的次数。
- `mapping_failures`：找到交点后无法定位反应元胞的候选次数。
- `no_local_intersections`：局部拟合没有有效前向交点的候选次数，包含反射后离开原表面的情况。
- `rank_deficient_fits`：局部拟合矩阵不满秩的次数。
- `undefined_normal_fits`：邻域质心与目标中心重合、无法定义旋转轴的拟合次数；对应元胞保留库存，暂不参与曲线碰撞。

## 扩展

新增表面模型实现 `Reaction::apply`，其参数包含碰撞角余弦和捕获概率；入射角和能量分布通过 `double(Random&)` 函数装入 `Source`；产额模型接收独立的 `AngularYield` 函数。元胞更新由 `CellUpdater` 负责。源粒子和表面事件使用独立随机数流。
