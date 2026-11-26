# 贡献指南

感谢你愿意花时间。这个库的价值不在功能多，而在**每一个数都能被复核**，所以下面几条
要求都是为了守住这一点。

## 1. 开发环境

```bash
git clone https://github.com/lizhengyang-byte/slopestab
cd slopestab
python -m venv .venv && source .venv/bin/activate
pip install -e ".[all,dev]"
```

要求 Python ≥ 3.9。核心依赖只有 `numpy`；`matplotlib`（出图/计算书）、`PyYAML`（yaml 案例）、
`pandas`（`to_dataframe`）都是可选的，代码里必须在函数内部延迟导入。

## 2. 本地校验

```bash
pytest -q                  # 单测 + 回归，约 2 s
pytest -m slow -q          # ACADS 基准、网格收敛
ruff check .               # 静态检查
ruff format --check .      # 格式
python -m build && twine check dist/*   # 打包自检
```

提交前这五条要全绿。CI 会跑同样的命令（3.9 / 3.11 / 3.13 三个版本）。

## 3. 改动数值结果时必须做的事

任何影响安全系数的改动（公式、条分、孔压、重度取值）都要同时更新：

1. `tests/benchmarks/test_acads.py` 里的期望值；
2. `docs/validation.md` 里对应的表格，以及 `README.md` / `README.en.md` 的验证表；
3. 如果偏差超过 ±1%，在 PR 描述里写清楚**为什么**是改动正确而不是引入了新错误
   （最好附一条能区分新旧实现的算例）。

只改期望值让它"变绿"的 PR 会被打回。

## 4. 加一个新方法

按现有方法的骨架走，五步：

1. **`methods/xxx.py`**：只接收 `SliceSet`，返回 `LEMResult`。不要碰几何、不要自己建切片。
2. **共用公式**：如果新方法也需要"由 X=λ f E 推有效法向力"，直接用
   `methods/gle.py` 的 `normal_forces` / `shear_forces`，别复制一份公式——库内一份公式
   只维护一处错误源。
3. **注册**：`analysis.py` 的 `_METHODS` 与 `METHOD_NAMES` 加上键（含常用别名）。
4. **残差**：调用 `equilibrium_residuals(slices, n_eff, shear)`，让别人能看见平衡差多少。
5. **测试**：至少三条——ACADS 上跑通、与最近的方法（Bishop/Spencer）差在合理带内、
   退化情形（φ=0 或 λ=0）能对上闭式。

## 5. 代码风格

- 注释与 docstring 用**英文**；用户文档用**中文**，但 README 必须双语：
  `README.md`（中文）与 `README.en.md`（英文）结构一一对应，改一处要同步另一处。
- 行宽 100，`ruff` 说了算；`l`（条底长）在 `geometry/slices.py` 里是行业习惯符号，
  该文件的 `E741` 已单独豁免。
- 不做"静默纠正"：几何非法、迭代发散、净驱动力 ≤0 一律抛异常，不许返回一个看起来
  正常的数。新增的容错必须同时给出可见的 `info["note"]`。
- 新依赖要极力避免。确实需要就放进可选 extras，并在函数内导入。

## 6. 提交与 PR

```
<type>: <一句话>

feat: 增加 Bishop 的孔压比输入
fix: 修正 GLE 中孔压项漏乘 cos α
docs: 补 ACADS 偏差归因
test: 增加三层土回归算例
```

PR 描述里请带：

- 一个最小可复现算例（`case.yaml` 或代码片段）；
- 改动前后的 Fs 对比；
- 如果与商业软件/文献不一致，说明差异归因。

## 7. 报告"算不对"

这比 star 更有价值。请附上：

1. 完整输入（坡面、土层、滑面、水位、条数）；
2. 你的手算或商业软件结果，以及出处（软件名+版本，或文献页码）；
3. `res.summary()` 的输出（含残差）。

如果只是"觉得不对"，请先说清楚期望值和依据——这样排查会快很多。

## 8. 不要做的事

- 不要把 `matplotlib`、`pandas` 提成必需依赖；
- 不要在核心模块里写 `print`（示例和 CLI 例外）；
- 不要给严格法加上"非圆弧也能算"的宽松分支——v0.1 的报错是有意的，见
  `docs/validation.md` §7；
- 不要提交 `examples/output/` 下的运行产物（已在 `.gitignore` 里）。
