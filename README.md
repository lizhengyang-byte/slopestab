# slopestab

**二维边坡稳定性分析 · 极限平衡法 · 零重依赖 · 可验证**

[![CI](https://github.com/lizhengyang-byte/slopestab/actions/workflows/ci.yml/badge.svg)](https://github.com/lizhengyang-byte/slopestab/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/slopestab.svg)](https://pypi.org/project/slopestab/)
[![Python](https://img.shields.io/pypi/pyversions/slopestab.svg)](https://pypi.org/project/slopestab/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`slopestab` 用一套几何内核实现五种经典极限平衡法（Fellenius / Bishop / Janbu /
Spencer / Morgenstern-Price），带滑面自动搜索、水压与拟静力地震、出图与 HTML 计算书，
并且**把每一个安全系数都对着文献算例验过**。

核心只依赖 `numpy`：装不上 matplotlib 也能算，算完仍能导出 CSV/JSON 喂给别的工具。

---

## 为什么再造一个轮子

| 现有选择 | 痛点 | slopestab 的取舍 |
|---|---|---|
| `pyslope` | 只有圆形滑面、只有 Bishop 一类方法 | 五种 LEM 同一套切片，可横向互校 |
| 商业软件（Slide/GeoStudio） | 不可编程、结果不可复现、无免费 CI 验证 | 纯 Python、结果可复现、benchmark 进 CI |
| 论文里的零散脚本 | 不成库、无测试、条分切点随意 | 层界/水线强制加密 + 平衡残差审计 |

设计目标很明确：**方法全、算得准、结果能直接拿去交付**。不追求做成有限元软件。

---

## 安装

```bash
pip install slopestab            # 核心（numpy）
pip install "slopestab[all]"     # + matplotlib（出图/计算书）+ PyYAML（yaml 案例）
```

从源码：

```bash
git clone https://github.com/lizhengyang-byte/slopestab
cd slopestab && pip install -e ".[all,dev]"
```

要求 Python ≥ 3.9，已在 3.9 / 3.11 / 3.13 上跑 CI。

---

## 30 秒上手

```python
import slopestab as sst

slope = sst.Slope([(0, 0), (20, 10), (30, 10)])  # 坡面折线，m
soil = sst.SoilModel([sst.SoilLayer("黏土", c=3.0, phi=19.6, gamma=20.0)])

circle = sst.SlipCircle(cx=0.03, cy=27.70, r=27.69)  # ACADS 1(a) 标准滑弧

for method in ["fellenius", "bishop", "janbu", "spencer", "morgenstern_price"]:
    res = sst.analyze(method, slope, soil, circle, n_slices=30)
    print(f"{method:20s} Fs = {res.fs:.4f}")

# 自动搜索最危险滑面
best = sst.search_circular(slope, soil, method="bishop")
print(best.fs, best.critical.surface)

sst.plot(best, "section.png")
sst.report_html(best, "report.html")
```

输出（同一算例，与文献值对照见下一节）：

```
fellenius             Fs = 0.9479
bishop                Fs = 0.9854
janbu                 Fs = 0.9427     (未修正)
spencer               Fs = 0.9844
morgenstern_price     Fs = 0.9841
```

![示例剖面](docs/images/example-section.svg)

---

## 方法一览

| 方法 | 平衡条件 | 条间力假设 | 备注 |
|---|---|---|---|
| `fellenius` | 力矩 | 忽略 | 最保守的包络下界，通常低 5–15% |
| `bishop` | 力矩 | 水平 | 工程最常用；与 Spencer 通常差 1–2% |
| `janbu` | 水平力 | 无剪力 | 需配 `f0` 修正（`analyze(..., f0=1.05)`） |
| `spencer` | 力矩 + 力 | 平行（`X/E=tanθ`） | 严格法，两个未知数 `(Fs, θ)` 同时解 |
| `morgenstern_price` | 力矩 + 力 | `X = λ·f(x)·E` | 最通用；`f(x)` 支持 constant/half-sine/sine/梯形 |

另有两个特殊工况的独立实现：

- `infinite_slope_fs()` —— 无限边坡闭式解（干燥 / 顺坡渗流），用作数值校验锚点；
- `planar_wedge()` —— 平面滑动楔体，含张裂缝、裂缝充水、拟静力地震。

---

## 验证：不是"跑通"，是对得上

`tests/benchmarks/` 里对着公开文献值做断言，CI 每次提交都跑。

**ACADS 1(a)**（Giam & Donald 1989；10 m 高 2H:1V 均质坡，c′=3 kPa，φ′=19.6°，γ=20 kN/m³，
标准滑弧 圆心 (0.03, 27.70)、R=27.69 m）

| 方法 | slopestab | 文献值 | 偏差 |
|---|---|---|---|
| Bishop 简化 | 0.9854 | 0.987 (Slide2) / 0.992 (Zhao 2014) | −0.002 / −0.007 |
| Spencer | 0.9844 | 0.986 (Slide2) | −0.002 |
| Morgenstern-Price | 0.9841 | 0.986 (Slide2 GLE) / 0.991 (Zhao) | −0.002 / −0.007 |
| Fellenius（普通条分） | 0.9479 | 0.967（另一滑弧） | 同量级 |
| Janbu 简化（未修正） | 0.9427 | — | 需 `f0` 修正后 ≈ 0.99 |
| ACADS 官方答案 | — | 1.000 | 全部落在 ±5% 内 |

另外三类自洽校验（都在单测里）：

1. **φ=0 恒等式**：此时 Bishop 必须与 Fellenius 完全一致（实测差 2×10⁻⁴），
   且 Fs 等于力矩解析解 `c·L·R / ΣW·(x−x₀)`；
2. **网格收敛**：条分数 10→80，Bishop 的 Fs 变化 < 0.001；
3. **严格法平衡闭合**：Spencer / M-P 的 ΣFₓ、ΣM 残差 < 0.05%（相对总重）。

详见 [`docs/validation.md`](docs/validation.md)。

---

## 文档

| 文档 | 内容 |
|---|---|
| [`docs/methods.md`](docs/methods.md) | 条分与符号约定、五种方法的平衡方程与公式、数值解法 |
| [`docs/validation.md`](docs/validation.md) | ACADS 基准、收敛性、恒等式校验、偏差归因 |
| [`docs/api.md`](docs/api.md) | 全部公开接口、参数、返回值、异常 |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | 开发与提交流程 |

---

## 结果里有什么

`analyze()` 返回 `LEMResult`，不只一个数字：

```python
res.fs  # 安全系数
res.slices.to_dataframe()  # 条分表（pandas 可选）
res.residuals  # {'Fx': -1.80, 'Fy': 0.00, 'M': 0.03} 平衡残差 %
res.info  # λ、θ、迭代次数、警告
res.to_csv("slices.csv")
res.to_json("result.json")
sst.report_html(res, "report.html")  # 离线可开的自包含计算书
```

**每个结果都带平衡残差审计**：Bishop/Fellenius 只有力矩平衡，力残差非零是正常的，
但你能看见它有多大，而不是被藏起来。

---

## 命令行

```bash
slopestab run    examples/cases/01_simple_slope.yaml --method bishop --plot out.png --report r.html
slopestab search examples/cases/02_layered_water.yaml --method spencer --csv all_circles.csv
slopestab methods
```

案例文件格式（YAML 或 JSON）：

```yaml
method: bishop
n_slices: 30
slope:
  profile: [[0, 0], [10, 0], [30, 10], [50, 10]]
soil:
  layers:
    - {name: 填土, c: 5, phi: 25, gamma: 18, gamma_sat: 19.5, boundary: [[0, 3], [50, 3]]}
    - {name: 黏土, c: 15, phi: 18, gamma: 19, gamma_sat: 20}
surface: {type: circle, cx: 10.03, cy: 27.70, r: 27.69}
water: {type: phreatic, points: [[10, 0], [30, 9], [50, 9]]}
seismic: {kh: 0.1, kv: 0.0}
surcharge: [{q: 20, x_start: 25, x_end: 50}]
```

---

## 工程约定（很重要，别踩坑）

- **单位强制 SI**：长度 m、应力 kPa、重度 kN/m³、角度 **度**。传错量纲不会自动救你，
  但构造时会做符号/范围校验。
- **滑动方向自动判定**：由滑面与坡面两个交点的高程决定，质量往低处滑；切片的
  `sin α > 0` 恒表示"沿滑动方向下降"，所以驱动力永远是 `ΣW sinα` 的正贡献。
- **条分切点强制加密**：土层层界、浸润线、超载边界与滑面的交点都会插入切点，
  保证每条底边只落在一个土层、一个孔压规则内。这是多数业余脚本算不准的主因。
- **水下重度分段积分**：不在整条上平均，逐段用 `γ` 或 `γ_sat`。
- **非物理结果直接报错**：不收敛抛 `ConvergenceError`（附最后一次迭代值），
  净驱动力 ≤ 0 抛 `ArithmeticError`，不会静默返回一个漂亮的错数。
- **可复现**：搜索是确定性网格 + 确定性精修，同样输入同样输出。

---

## 已知限制

- **严格法只支持圆形滑面**。非圆弧滑面的 Spencer/M-P 需要力矩参考点与条间力线位置的
  额外假设，v0.1 不做，以免给出"看起来对"的数；非圆弧滑面请用
  `Slope` + `PolylineSurface` 配 `bishop` / `janbu`。（v0.2 规划中）
- **Janbu 的 `f0` 不自动查表**：默认 `f0=1.0` 即未修正值，会系统性偏低；做设计请用
  Spencer / M-P，或自己按 Janbu 图表给 `f0`。
- **不做**：有限元强度折减（SRM）、三维楔体、渗流场求解、动力时程、GIS 接入。
- 浸润线是**静水压**假设（`u = γw·h`），非稳定渗流场需先用其他工具算好孔压再传入。

---

## Roadmap

- v0.2：非圆弧滑面的严格 GLE（Janbu 严格法 + M-P 任意 f(x)）、加筋（土钉/锚杆/格栅）
- v0.3：可靠度（Fs 的蒙特卡罗/一阶二阶矩）、Janbu `f0` 自动查表
- 长期：QGIS/GeoPandas 剖面接入、报告模板定制

---

## English summary

`slopestab` is a numpy-only Python library for two-dimensional slope stability
analysis by limit equilibrium: Fellenius (ordinary), Bishop simplified, Janbu
simplified, Spencer and Morgenstern-Price on a shared slice kernel, plus a
circular slip-surface search, infinite-slope and planar-wedge closed forms,
pseudo-static seismic loading, pore-pressure models and self-contained HTML
reporting.

Every result carries a global force/moment imbalance audit, and the suite is
verified against the ACADS 1(a) benchmark (Bishop 0.985 vs published 0.987,
Spencer 0.984 vs 0.986) and against exact identities (undrained moment
solution, mesh convergence). CI runs on Python 3.9/3.11/3.13.

Scope is deliberately limited: no FEM strength reduction, no 3D wedges, and
rigorous GLE is restricted to circular surfaces in v0.1.

---

## 贡献 / 引用 / 许可

- 发现算例对不上？请带上 `case.yaml` 开 issue，这比任何 star 都值钱。
- 贡献指南见 [`CONTRIBUTING.md`](CONTRIBUTING.md)。
- 学术引用见 [`CITATION.cff`](CITATION.cff)。
- MIT License。**计算结果须由注册岩土工程师复核后方可用于设计。**

作者：李尚儒（成都锦城学院）· 2026
