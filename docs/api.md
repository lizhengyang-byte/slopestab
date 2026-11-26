# API 参考

`slopestab` 的公开接口都在 `import slopestab as sst` 之后可用。本页按模块列出签名、
参数含义与常见用法；每个类/函数的 docstring 在源码里更细。

约定：

- 长度 m、应力 kPa、重度 kN/m³、角度**度**；
- 所有构造器都做范围校验，非法输入立刻抛 `ValueError`，不会算出一个漂亮的错数；
- `numpy` 是唯一必需依赖，`matplotlib` / `PyYAML` / `pandas` 都在用到时才导入。

---

## 1. 几何

### `Slope(profile)`

坡面折线。`profile` 是 `[(x, y), ...]`，x 严格递增。

```python
slope = sst.Slope([(0, 0), (10, 0), (30, 10), (50, 10)])
slope.height      # 最大高差
slope.x_min, slope.x_max
slope.y_at(15.0)  # 折线插值，两端做钳制
```

`Slope.simple(height=10.0, ratio=2.0, crest_length=10.0, toe_length=5.0)` 生成
"坡脚平台 + 单坡 + 坡顶平台"的标准断面。`ratio` 是水平:垂直，`ratio=2.0` 即 2H:1V。

### `SoilLayer(name, c, phi, gamma, gamma_sat=None, boundary=None)`

一层土。`c` 是**有效**黏聚力（kPa），`phi` 是有效内摩擦角（度），`gamma` 是水上重度，
`gamma_sat` 是水下重度（缺省时全用 `gamma`，即总应力/干燥分析）。

`boundary` 是这一层的**底面**折线；层序自上而下排列，最后一层不需要 `boundary`
（向下无限延伸）。

```python
soil = sst.SoilModel([
    sst.SoilLayer("填土", c=5, phi=25, gamma=18, gamma_sat=19.5, boundary=[(0, 3), (50, 3)]),
    sst.SoilLayer("黏土", c=15, phi=18, gamma=19, gamma_sat=20, boundary=[(0, 7), (50, 7)]),
    sst.SoilLayer("基岩", c=50, phi=30, gamma=21),
])
```

`SoilModel` 提供 `layer_at(x, y)`、`gamma_at(x, y, below_water)`、`boundaries()`，
条分时按位置取值，不做整条平均。

### 滑面

```python
sst.SlipCircle(cx, cy, r)                    # 圆弧
sst.PolylineSurface([(x0, y0), (x1, y1), ...])  # 折线（任意形状）
```

两者都实现同一套接口：

| 方法/属性 | 说明 |
|---|---|
| `intersections(slope)` | 与坡面的交点（按 x 排序），滑动方向据此判定 |
| `base_y(x)` | 滑面高程，支持数组 |
| `base_slope(x)` | 滑面切线斜率 |
| `length(x1, x2)` | 两点间的弧长/折线长 |
| `sample(n, span=(x0, x1))` | 采样点，供出图与积分 |
| `moment_center()` | 力矩参考点（圆心 / 折线均值点） |
| `x_span()` | 滑面覆盖的 x 范围 |

`SlipCircle` 上还能直接读 `cx`、`cy`、`r`；`PolylineSurface` 读 `points`。

### `build_slices(slope, soil, surface, water=None, surcharges=None, seismic=None, n_slices=30)`

生成条分表，返回 `SliceSet`。一般不用直接调，`analyze()` 内部会调；需要自己改切点时用它。

`SliceSet` 是纯数据对象（NumPy 数组），字段：

| 字段 | 含义 |
|---|---|
| `x_left, x_right, b, x_mid` | 条的边界与宽度 |
| `y_base, y_ground, h` | 条底高程、地表高程、条高 |
| `l, c, phi, tan_phi` | 底边长、该条黏聚力与内摩擦角（按层取值） |
| `w, w_eff` | 条重与含超载的有效条重（kN/m） |
| `kh_force` | 水平惯性力（kN/m，已含滑动方向符号；竖向分量并入 `w_eff`） |
| `u` | 条底孔压（kPa） |
| `sin_a, cos_a` | 条底倾角的正余弦（见 `docs/methods.md` §2 的符号约定） |
| `ux, uy, nx, ny` | 单位切向/法向量 |
| `x_ref, y_ref` | 力矩参考点（圆弧取圆心，折线取均值点） |
| `y_c, area` | 条重心纵坐标与条面积（m²） |
| `layer_index, water_table` | 每条落在哪一层、该位置的浸润线高程（无水为 NaN） |
| `slide_dir` | 滑动方向符号，`-1` 表示物质向左下运动 |
| `notes` | 切点加密等过程中的说明 |
| `surface_kind` | `"circle"` 或 `"polyline"` |

辅助方法：`.as_table()`（dict 列表）、`.to_dataframe()`（pandas）、`.total_weight`、`.span`。

---

## 2. 荷载

| 类 | 用途 |
|---|---|
| `NoWater()` | 无水（默认） |
| `RuCoefficient(ru)` | 孔压比 \(u=\rho_u\gamma z\)，快速估算 |
| `PhreaticLine(points, gamma_w=9.81)` | 浸润线，静水压 \(u=\gamma_w(h_{water}-y_{base})\) |
| `Seismic(kh, kv=0.0)` | 拟静力地震，`kh` 沿滑动方向，`kv` 向下为正 |
| `UniformSurcharge(q, x_start=-inf, x_end=+inf)` | 均布超载（kPa），可给多个；两端缺省即全坡面 |
| `SurchargeSet([...])` | 超载集合，多个时用 |

```python
water = sst.PhreaticLine([(10, 0), (30, 9), (50, 9)])
seismic = sst.Seismic(kh=0.15)
surcharge = sst.UniformSurcharge(q=30, x_start=25, x_end=50)
```

孔压约定：压为正，作用在条底法向。**浸润线不得高于地面**：一旦某条的 \(u b\) 超过条重 \(W\)，
条分阶段就抛 `ValueError`，报出具体 x、\(u b\) 与 \(W\) 的数值——这是有意设计的护栏，
因为此时有效法向力已为负，算出的 Fs 没有物理意义。

---

## 3. 分析

### `analyze(method="bishop", slope=None, soil=None, surface=None, water=None, surcharge=None, surcharges=None, seismic=None, n_slices=30, **kwargs)`

一次分析，返回 `LEMResult`。

```python
res = sst.analyze("spencer", slope, soil, circle, water=water, n_slices=30)
```

- `method`：`"fellenius"`（别名 `"ordinary"`、`"swedish"`）、`"bishop"`、`"janbu"`、
  `"spencer"`、`"morgenstern_price"`（别名 `"mp"`、`"gle"`）。
- `kwargs` 透传给方法：`janbu` 的 `f0`、`morgenstern_price` 的 `f`（形状名或数组）、
  `lam_range`、`n_scan`，以及迭代控制 `max_iter`、`tol`。

### `compare(methods=METHOD_NAMES, slope, soil, surface, **kwargs)`

对**同一套条分**跑多种方法，返回 `{method: LEMResult}`。这是横向互校的正确姿势：

```python
results = sst.compare(["bishop", "spencer", "morgenstern_price"], slope, soil, circle)
for name, r in results.items():
    print(f"{name:22s} {r.fs:.4f}  {r.residuals['M']:+.2f}%")
```

### `run_on_slices(method, slices, **kwargs)`

跳过几何构建，直接对已有的 `SliceSet` 求解。搜索或自定义切点时会用到。

### `LEMResult`

| 属性 | 说明 |
|---|---|
| `method, fs` | 方法名与安全系数 |
| `slices` | 参与计算的 `SliceSet` |
| `normal_forces, shear_forces` | 每条的有效法向力与剪力（kN/m） |
| `residuals` | `{'Fx','Fy','M','Fx_abs','Fy_abs','M_abs'}`，前三个是相对总重的百分比 |
| `converged, iterations` | 是否收敛、迭代次数 |
| `info` | 方法特有信息：`lambda`、`theta_deg`、`f0`、`force_moment_gap`、`note` |
| `slope, soil, surface, water` | 可选挂载，出图/出报告时用 |

方法：`.summary()`、`.to_csv(path)`、`.to_json(path)`。

`analyze()` **不会**自动把 `slope/soil/surface` 挂到结果上（避免隐式依赖），
要出图或出报告时手动赋值，或直接用 `search_circular` 的返回对象。

---

## 4. 五种方法（直接调用）

```python
sst.fellenius(slices)
sst.bishop(slices, max_iter=200, tol=1e-4)
sst.janbu(slices, f0=1.05, max_iter=200, tol=1e-4)
sst.spencer(slices, fs0=None, lam0=0.3, tol=1e-6, max_iter=80)
sst.morgenstern_price(slices, f="half_sine", lam_range=(-0.2, 1.2), n_scan=31, tol=1e-5)
```

它们都接收 `SliceSet`，返回 `LEMResult`。`spencer` / `morgenstern_price` 在非圆弧滑面上
抛 `GeometryError`（v0.1 限制，见 `docs/validation.md` §7）。

M-P 的 `f` 支持 `"constant"`、`"half_sine"`、`"trapezoid"`、`"sine"` 或长度等于条数的数组。
`f="constant"` 时结果与 `spencer` 完全一致，可作为自检。

---

## 5. 滑面搜索

### `SearchGrid(cx_range, cy_range, r_range, n_cx=24, n_cy=14, n_r=12)`

`SearchGrid.default(slope)` 由坡高 \(H\) 推导：`cx_range = (x_min − 0.1H, x_min + 2.6H)`，
`cy_range = (y_max + 0.05H, y_max + 2.8H)`，`r_range = (0.5H, 3.2H)`。默认网格 24×14×12 = 4032 个圆。

### `search_circular(slope, soil, method="bishop", grid=None, n_slices=20, refine=True, refine_candidates=3, refine_n_slices=30, max_refine_iter=60, water=None, surcharges=None, seismic=None, progress=None)`

返回 `SearchResult`：

| 属性 | 说明 |
|---|---|
| `critical` | 临界滑面的 `LEMResult` |
| `fs` | 临界安全系数（快捷属性） |
| `candidates` | `[(cx, cy, r, fs), ...]`，粗网格排序结果 |
| `grid, method, elapsed, n_evaluated` | 搜索配置与开销 |
| `notes` | 预筛/精修过程中的说明 |

`progress(done, total)` 可挂进度回调。搜索是确定性的：同样输入同样输出。

```python
best = sst.search_circular(slope, soil, method="spencer", n_slices=25)
print(best.fs, best.critical.surface)
best.to_dataframe().to_csv("candidates.csv")
```

注意 `n_evaluated` 是**有效试算数**（几何预筛之后），不是网格总数。

---

## 6. 闭式解

### `infinite_slope_fs(beta_deg, z, c=0.0, phi=30.0, gamma=18.0, gamma_sat=None, gamma_w=9.81, seepage=False)`

无限边坡：`seepage=False` 用干燥式，`True` 用顺坡稳定渗流式。`z` 是滑面深度（m）。

### `critical_depth(beta_deg, c, phi, gamma, gamma_sat=None, gamma_w=9.81, seepage=False)`

反解 \(F=1\) 的深度；强度足够（永远不破坏）时返回 `None`。

### `planar_wedge(slope, soil, x_toe, x_top, crack_depth=0.0, water_depth=0.0, seismic=None, gamma_w=9.81)`

平面楔体，返回 `WedgeResult`（`fs`、`weight`、`plane_length`、`water_force_plane`、
`water_force_crack`、`summary()`）。`crack_depth` 是坡顶张裂缝深度，`water_depth` 是缝内水深。

---

## 7. 出图与计算书

### `plot(result, path=None, *, slope=None, soil=None, water=None, figsize=(9, 5.5), dpi=150, show_slices=True, title=None, ax=None)`

`result` 可以是 `LEMResult` 或 `SearchResult`。不给 `path` 时返回 matplotlib 的 `Figure`
（便于嵌到自己的版面）；给 `path` 则按扩展名存 PNG/PDF/SVG 并关闭 figure。

剖面里画：土体、层界、浸润线与水压力区、条分线、滑面、滑动体高亮。matplotlib 在函数内
导入，未安装时抛 `ImportError` 而不是污染核心依赖。

### `report_html(result, path, *, title="边坡稳定分析计算书", include_plot=True, max_rows=60)`

生成**自包含**的一页 HTML 计算书：结论卡片（Fs 与分级）、方法参数、条分表（可截断）、
平衡残差、剖面图（内嵌 base64）。不依赖网络，可以邮件发出去或直接归档。

---

## 8. 案例文件

```python
case = sst.load_case("case.yaml")        # YAML 或 JSON → dict
problem = sst.build_problem(case)        # dict → {slope, soil, surface, water, ...} 对象
sst.dump_case(case, "copy.yaml")         # 回写
```

`build_problem` 的返回可直接解包给 `analyze`：

```python
p = sst.build_problem(sst.load_case("case.yaml"))
res = sst.analyze("bishop", p["slope"], p["soil"], p["surface"], water=p.get("water"))
```

字段格式见 README 的"命令行"一节。

---

## 9. 命令行

```
slopestab run    <case> --method bishop --plot out.png --report r.html [--csv s.csv] [--json r.json]
slopestab search <case> --method spencer --csv all_circles.csv [--plot best.png] [--report b.html]
slopestab methods
```

`run` 走单滑面分析，`search` 走滑面搜索；两者都能同时出图与出报告。返回值即进程退出码，
便于脚本串联。

---

## 10. 异常

| 异常 | 何时抛 | 怎么办 |
|---|---|---|
| `SlopestabError` | 基类 | — |
| `GeometryError` | 几何非法：滑面与坡面不相交、非圆弧配严格法 | 检查滑面与坡线范围 |
| `ConvergenceError` | 迭代到非物理 Fs（≤ 0）或 M-P 找不到 λ | 调大 `max_iter`、加条数，或换方法/重选滑面 |
| `ArithmeticError` | 净驱动力 ≤0，滑面不是可动机构 | 换滑面，或确认滑动方向 |
| `ValueError` | 构造参数越界（c<0、phi≥90、profile 未递增等）；孔压超过条重（浸润线高于地面，报 x 与 `u*b`、`W`） | 按提示修正输入，或降低水位 |

`ConvergenceError` 带 `last_fs` 属性（最后一次迭代值），便于判断是"差一点"还是"完全跑飞"。
