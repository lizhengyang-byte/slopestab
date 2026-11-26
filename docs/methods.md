# 计算方法

本文说明 `slopestab` 里每个安全系数是怎么算出来的：几何与条分、荷载、平衡方程、
五种方法的假设与公式，以及数值求解与残差审计。目标是让人能拿着这份文档复现每一个
数字，也便于指出哪一步有问题。

---

## 1. 总体流程

一次分析固定走这几步，五种方法只在第 5 步分岔：

1. 读入坡面折线 `Slope`（x 严格递增）与土层 `SoilModel`；
2. 给定滑面（圆弧 `SlipCircle` 或折线 `PolylineSurface`）；
3. 求滑面与坡面的交点，用这两个点确定滑动方向，并布置条分切点；
4. 逐条算几何量、重度、孔压、超载与地震力（`SliceSet`）；
5. 按方法组装平衡方程，解出安全系数 `Fs`；
6. 用同一套条分复核整体力/力矩残差，写进 `LEMResult.residuals`。

第 3–4 步是五种方法共用的。**同一算例换方法时条分表完全一致**，所以方法之间的差异
只来自平衡假设，不来自切点位置——`compare()` 就是靠这一点做横向互校的。

---

## 2. 符号与方向约定

这一节是最容易出错的地方，实现里所有公式都按它写。

- 单位强制 SI：长度 m，应力 kPa，重度 kN/m³，角度输入用**度**，内部转弧度。
- `cos α > 0`：条底法向总是指向"土体内部"一侧。
- `sin α > 0` 表示沿滑动方向**下降**。滑动方向由滑面两端与坡面交点的高程决定：
  右端更高则物质向左下滑，`slide_dir = -1`，反之为 `+1`。
- 单位切向量 \(\mathbf{u} = (\ sd\cos\alpha,\ -\sin\alpha\ )\)，单位法向量
  \(\mathbf{n} = (\ sd\sin\alpha,\ \cos\alpha\ )\)，`sd` 即 `slide_dir`。
- 于是驱动力恒为 \(\sum W\sin\alpha\)（`sin α > 0`），抗滑力恒沿 \(-\mathbf{u}\)。
- 孔压 \(u\) 以**压为正**，沿 \(+\mathbf{n}\) 作用在条底，合力 \(U = u\,l\)。

条间力约定：\(E\) 为水平条间力（指向滑动方向为正），\(X\) 为竖向条间力（向下为正），
上标 \(i,i+1\) 表示条的左右侧。

---

## 3. 条分与几何

### 3.1 强制切点

滑面按等宽切分只是起点，真正的切点集合是下面这些 x 坐标的并集：

- 等分点；
- **坡面折点**（坡顶、坡脚、平台端）；
- **土层层界与滑面的交点**——把层界折线当成一个 `Slope`，复用
  `surface.intersections()` 求交，避免自己写重复的求交逻辑；
- **浸润线与滑面的交点**；
- **超载边界**（`x_start`、`x_end`）。

理由很实际：一条底边若跨越两个土层，`c`、`φ` 该取谁就说不清了；跨过水线则孔压不连续。
加密之后每条底边只落在一个土层、一个孔压区间内，这是多数业余脚本算不准的主因。

### 3.2 条的几何量

对第 \(i\) 条，宽 \(b_i = x_{i+1}-x_i\)，底边与滑面重合并近似为直线：

\[
l_i = \sqrt{b_i^2 + (y_{i+1}-y_i)^2},\qquad
\sin\alpha_i = \frac{-(y_{i+1}-y_i)}{l_i}\cdot(\pm),\qquad
\cos\alpha_i>0
\]

条高 \(h_i\) 由坡面（或其上的超载面）与条底的高差给出；条重按**层界与水线分段积分**：

\[
W_i = b_i\sum_k \gamma_k\,h_{ik},\qquad W_i^{eff}=W_i + q_i b_i
\]

其中水下部分用 \(\gamma_{sat}\)、水上用 \(\gamma\)，不在整条上取平均。重心纵坐标由
\(\gamma A y\) 加权得到，供力矩方程使用。

### 3.3 条重与孔压的相互作用

计算有效法向力时到处出现 \(u b\) 这一项（\(b\) 是条宽，不是底边长），因为孔压竖直分量
等于 \(u\,l\cos\alpha = u\,b\)。实现里所有方法都按这个量写，五种方法因此可以互相校验。

---

## 4. 荷载

| 荷载 | 实现 | 进入方程的方式 |
|---|---|---|
| 自重 | 分层积分 | 竖向 \(W\)，作用在条重心 |
| 孔压 | `NoWater` / `RuCoefficient` / `PhreaticLine` | 法向合力 \(U=u\,l\)，压力为正 |
| 超载 | `UniformSurcharge`，可多条 | 并入 \(q_i b\)，切点在其边界加密 |
| 地震 | `Seismic(kh, kv)` | 水平惯性力 \(k_h W\) 沿滑动方向，竖向 \(k_v W\) 向下为正（不利） |

地震力取绝对值符号约定最省事：水平惯性力永远指向滑动方向（拟静力法的保守取法），
所以 `Fellenius` 的 \(\sum W\sin\alpha + k_h\cos\alpha\) 与实现一致。

---

## 5. 条间力与通用力平衡

除 Fellenius 外，四种方法都给条间力一个假设。把第 \(i\) 条单独取出来，作用在它上面的力是
自重（含超载）、条底有效法向力 \(N'\)、条底剪力 \(S\)、孔压合力 \(U\)、水平惯性力，
以及左右两侧的 \(E\)、\(X\)。

剪力按莫尔–库仑强度 mobilized 到安全系数 \(F\)：

\[
S_i = \frac{c_i' l_i + N_i'\tan\varphi_i'}{F}
\]

竖向与水平平衡写作（`sd` 为滑动方向符号）：

\[
\Delta X_i = W^{eff}_i - (N_i'+U_i)\cos\alpha_i - S_i\sin\alpha_i
\]

\[
\Delta E_i = sd\left[\,-k_hW_i - (N_i'+U_i)\sin\alpha_i + S_i\cos\alpha_i\right]
\]

对严格法再引入条间力形状假设

\[
\Delta X_i = \lambda\, f(x_i)\,\Delta E_i
\]

把 \(S_i\) 展开、以 \(N_i'\) 为未知整理，就得到显式的有效法向力（这是 Spencer / M-P 共用的一步，
见 `methods/gle.py`）：

\[
N_i'=\frac{
W^{eff}_i-u_ib_i+\lambda f_i\,sd\,k_hW_i+\lambda f_i\,sd\,u_il_i\sin\alpha_i
-\dfrac{c_i'l_i}{F}\left(\sin\alpha_i+\lambda f_i\,sd\cos\alpha_i\right)
}{
\cos\alpha_i-\lambda f_i\,sd\sin\alpha_i
+\dfrac{\tan\varphi_i'}{F}\left(\sin\alpha_i+\lambda f_i\,sd\cos\alpha_i\right)
}
\]

`λ = 0` 时它严格退化为 Bishop 的闭式（单测里直接比对过，误差 0）。这条式子是全库
最容易出错的地方——最初漏掉 \(u b\) 与 \(u l\sin\alpha\) 两个孔压项，导致有水的算例
Spencer 比 Bishop 高 22%，现已修正并加进 benchmark。

---

## 6. 五种方法

### 6.1 Fellenius（普通条分 / 瑞典条分）

假设条间力整体忽略，每条只做径向平衡：

\[
N_i' = W^{eff}_i\cos\alpha_i - k_hW_i\sin\alpha_i - u_il_i
\]

\[
F=\frac{\sum\left(c_i'l_i+N_i'\tan\varphi_i'\right)}{\sum\left(W_i\sin\alpha_i+k_hW_i\cos\alpha_i\right)}
\]

只满足力矩平衡，系统性偏低 5–15%（孔压大时更多）。定位是**保守筛查下界**，不做设计值。

### 6.2 Bishop 简化

假设条间力水平（\(X=0\)），对圆心取矩。经典写法把 \(N'\) 藏进 \(m_\alpha\)：

\[
m_\alpha=\cos\alpha_i+\frac{\sin\alpha_i\tan\varphi_i'}{F}
\]

\[
F=\frac{\sum\dfrac{c_i'b_i+\left(W^{eff}_i-u_ib_i-k_hW_i\sin\alpha_i\right)\tan\varphi_i'}{m_\alpha}}
{\sum\left(W_i\sin\alpha_i+k_hW_i\cos\alpha_i\right)}
\]

这里分子用 \(c'b\) 而不是 \(c'l\)：\(m_\alpha\) 的推导已经把 \(l\) 吸收进去了，写成 \(c'l/m_\alpha\)
会多算一个 \(1/\cos\alpha\)，早期版本就是这么错成 0.8765 的。收敛后按

\[
N_i'=\frac{W^{eff}_i-u_ib_i-k_hW_i\sin\alpha_i-\dfrac{c_i'l_i}{F}\sin\alpha_i}{m_\alpha}
\]

回算有效法向力，用于残差与出图。迭代是定点迭代加 0.5 欠松弛，初值取 Fellenius。

### 6.3 Janbu 简化（含 \(f_0\) 修正）

假设条间力水平（无剪力），做**水平力**平衡，得

\[
n_\alpha=\cos^2\alpha_i\left(1+\tan\alpha_i\tan\varphi_i'/F\right)
\]

\[
F=f_0\,\frac{\sum\dfrac{c_i'b_i+\left(W^{eff}_i-u_ib_i-k_hW_i\sin\alpha_i\right)\tan\varphi_i'}{n_\alpha}}
{\sum\left(W_i\tan\alpha_i+k_hW_i\right)}
\]

未修正的简化 Janbu 一般比严格法低 5–20%。\(f_0\) 取决于滑面形状与强度发挥程度，
本库**不自动查表**：默认 `f0=1.0`（即未修正），要修正就显式传 `f0`。ACADS 1(a)
上 `f0=1.05` 得到 0.9987，与官方答案 1.000 吻合。

### 6.4 Spencer

假设条间合力方向处处相同：\(X_i/E_i=\tan\theta=\lambda\)，等价于 \(f(x)\equiv1\)。
未知量是 \((F,\theta)\)，需要同时满足**力平衡**与**力矩平衡**：

\[
\sum\Delta E_i(F,\lambda)=0,\qquad \sum M_i(F,\lambda)=0
\]

实现上把两式当成二维非线性方程组，用数值 Jacobian 的 Newton 迭代（带线搜索阻尼）求解，
初值取 Bishop 的 \(F\) 与 \(\theta=\arctan(\lambda)\approx 0\)。ACADS 1(a) 三次迭代收敛，
\(\sum F_x\)、\(\sum M\) 残差均到 0。

### 6.5 Morgenstern-Price

把条间力形状放成一般式 \(X=\lambda f(x)E\)，\(f(x)\) 支持：

| 名称 | 形状 | 说明 |
|---|---|---|
| `constant` | \(f\equiv1\) | 等价于 Spencer，可用来互校 |
| `half_sine` | \(\sin(\pi t)\)，\(t=(x-x_{left})/L\) | 默认；两端为 0，最常用 |
| `trapezoid` | 中段 1、两端线性过渡 | 折减两端条间力 |
| `sine` | \(-\sin(2\pi t)\) | 全波，部分区间反号；仅作对照 |

\(f(x)\) 也可直接传长度等于条数的数组。

对每个 \(\lambda\) 分别解两个一元方程：力矩平衡给 \(F_m(\lambda)\)，力平衡给 \(F_f(\lambda)\)，
再找使 \(F_m-F_f=0\) 的 \(\lambda\)。实现流程：

1. 在 `lam_range`（默认 \(-0.2\sim1.2\)）上扫 31 个点，剔除无解点；
2. 找变号区间，用**割线/二分混合**迭代收敛（试值落到无解区时自动回退，不硬碰）；
3. 若整段都不变号但 \(|F_m-F_f|\) 的最小值 \(\le 0.02\)，接受该 \(\lambda\) 并在
   `info["note"]` 写明"没有严格夹逼"；否则抛 `ConvergenceError`。

第 3 条是实证结果：三层土 + 浸润线的算例里，`half_sine` 的残差最小约 0.002，
再细扫也夹不出零点——这时给一个诚实的近似解比报错更合理，但必须让它可见
（`converged=False` + `note`）。

---

## 7. 特殊工况的闭式解

### 7.1 无限边坡

干燥：

\[
F=\frac{c'}{\gamma z\sin\beta\cos\beta}+\frac{\tan\varphi'}{\tan\beta}
\]

顺坡稳定渗流（水位与坡面平行、水头沿坡向线性下降）：

\[
F=\frac{c'}{\gamma_{sat}z\sin\beta\cos\beta}
+\frac{\gamma'}{\gamma_{sat}}\cdot\frac{\tan\varphi'}{\tan\beta}
\]

`critical_depth()` 反解 \(F=1\) 的深度；若强度足够则返回 `None`。

### 7.2 平面楔体

给定坡脚与坡顶横坐标，楔体几何由坡面与滑面两条直线围成，可选张裂缝深度与缝内水深。
力系为 \(W\)、滑面 \(N\)、\(S\)、缝内水压力 \(P_w=\tfrac12\gamma_w h_w^2\) 与惯性力 \(k_hW\)，
沿滑面方向与法向分解后求 \(F\)。与条分法同构，因此可用来验证"条数取一半"时量级不漂移。

---

## 8. 数值实现

全库只依赖 NumPy，不引 SciPy，自己实现的三个小工具都在 `methods/` 内：

- `bisect`：一元方程求根，要求端点变号；不变号直接报错而不是悄悄返回端点值；
  这是 M-P 依赖的行为，不许被"宽容化"。
- `_newton_2d`：Spencer 的二维 Newton，数值 Jacobian + 线搜索（步长对半折减直到残差下降）。
- `_nelder_mead`：搜索阶段的局部精修，单形反射/扩张/收缩/压缩，无梯度。

不收敛时有三种明确结果，不会静默给错数：

| 情形 | 抛出 |
|---|---|
| 迭代到非物理 \(F\le0\)（孔压超过条重） | `ConvergenceError`（提示检查水位是否高于地面） |
| 净驱动力 \(\le0\)（不是可动机构） | `ArithmeticError` |
| 力矩方程无定义（非圆弧配严格法） | `GeometryError` |

---

## 9. 平衡残差审计

`LEMResult.residuals` 里有两组数：以"总重量"为基准的百分比 `Fx`、`Fy`、`M`，
以及绝对值 `Fx_abs`、`Fy_abs`、`M_abs`（kN/m、kN·m/m）。

\[
F_x=\sum\left[(N_i'+U_i)n_{ix}-S_iu_{ix}\right]+\sum k_hW_i,\qquad
F_y=\sum\left[(N_i'+U_i)n_{iy}-S_iu_{iy}\right]-\sum W^{eff}_i
\]

力矩以 \((x_{ref},y_{ref})\) 为极点算：

\[
M=\sum\Big[(N_i'+U_i)\big[(r_{ix})(n_{iy})-(r_{iy})(n_{ix})\big]
-S_i\big[(r_{ix})(u_{iy})-(r_{iy})(u_{ix})\big]+k_hW_i\,(y_{ref}-y_{c,i})\Big]
\]

力矩参考点对圆弧取圆心，对折线滑面取折线均值点。Bishop / Fellenius 只满足力矩平衡，
\(F_x\) 非零是正常的——但你能看见它有多大，而不是被藏起来。严格法要求 \(F_x\)、\(M\)
都落到 0 附近，benchmark 里以相对总重 <0.05% 作为断言。

---

## 10. 参考

- Fellenius, W. (1936). *Calculation of the stability of earth dams.*
- Bishop, A. W. (1955). The use of the slip circle in the stability analysis of slopes. *Géotechnique* 5(1).
- Janbu, N. (1954, 1973). *Slope stability computations.*
- Spencer, E. (1967). A method of analysis of the stability of embankments assuming parallel interslice forces. *Géotechnique* 17(1).
- Morgenstern, N. R., & Price, V. E. (1965). The analysis of the stability of general slip surfaces. *Géotechnique* 15(1).
- Fredlund, D. G., & Krahn, J. (1977). Comparison of slope stability methods of analysis. *Canadian Geotechnical Journal* 14(3).
- Giam, S. K., & Donald, I. B. (1989). *ACADS* 1(a) benchmark problem.
