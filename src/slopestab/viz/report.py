"""Self-contained HTML calculation report (inline CSS, base64 figure)."""

from __future__ import annotations

import base64
import datetime as _dt
from pathlib import Path
from typing import Any

_CSS = """
:root { --ink:#1d2433; --muted:#6b7280; --line:#d7dbe4; --accent:#0f5fa6; }
* { box-sizing: border-box; }
body { margin:0; padding:32px 28px 56px; font-family:-apple-system,"Segoe UI",
       "Helvetica Neue","PingFang SC","Microsoft YaHei",sans-serif;
       color:var(--ink); background:#fff; line-height:1.65; }
.wrap { max-width: 940px; margin: 0 auto; }
h1 { font-size: 24px; margin: 0 0 4px; }
h2 { font-size: 17px; margin: 30px 0 10px; padding-bottom:6px;
     border-bottom:1px solid var(--line); }
.meta { color:var(--muted); font-size:13px; margin-bottom:22px; }
.cards { display:flex; gap:14px; flex-wrap:wrap; margin:18px 0 6px; }
.card { flex:1 1 150px; border:1px solid var(--line); border-radius:10px;
        padding:12px 14px; background:#fafbfd; }
.card .k { font-size:12px; color:var(--muted); letter-spacing:.02em; }
.card .v { font-size:20px; font-weight:600; margin-top:4px; }
.card .v small { font-size:12px; font-weight:400; color:var(--muted); }
.fs-high { color:#1b7f4b; } .fs-low { color:#c62828; } .fs-mid { color:#b26a00; }
img { width:100%; border:1px solid var(--line); border-radius:10px; margin-top:12px; }
table { border-collapse:collapse; width:100%; font-size:12.5px; margin-top:10px; }
th, td { border-bottom:1px solid var(--line); padding:6px 8px; text-align:right; }
th { background:#f2f5f9; font-weight:600; }
th:first-child, td:first-child { text-align:left; }
.note { background:#f6f8fb; border-left:3px solid var(--accent); padding:10px 14px;
        font-size:13px; color:#3b4452; margin-top:18px; }
footer { margin-top:36px; color:var(--muted); font-size:12px;
         border-top:1px solid var(--line); padding-top:12px; }
"""


def report_html(
    result: Any,
    path: str | Path,
    *,
    title: str = "边坡稳定分析计算书",
    include_plot: bool = True,
    max_rows: int = 60,
) -> Path:
    """Write a standalone HTML calculation report for one analysis."""
    from .. import __version__

    lem, search = _unpack(result)
    slices = lem.slices
    fs = float(lem.fs)
    fs_class = "fs-high" if fs >= 1.5 else ("fs-mid" if fs >= 1.0 else "fs-low")

    img_tag = ""
    if include_plot:
        img_tag = _inline_plot(result, lem)

    cards = _cards(lem, slices, fs, fs_class)
    summary_rows = _summary_rows(lem)
    residual_rows = "".join(
        f"<tr><td>{key}</td><td>{value:+.4f}</td></tr>"
        for key, value in lem.residuals.items()
        if not key.endswith("_abs")
    )
    slice_rows = _slice_rows(slices, max_rows)

    geometry = ""
    slope = getattr(search, "slope", None) or getattr(lem, "slope", None)
    if slope is not None:
        pts = " ".join(f"({x:.2f}, {y:.2f})" for x, y in slope.profile)
        geometry = f"<p class='meta'>坡面折线（m）：{pts}</p>"

    search_note = ""
    if search is not None:
        search_note = (
            f"<div class='note'>滑面搜索：{search.n_evaluated} 个圆弧试算，"
            f"耗时 {search.elapsed:.2f} s，目标函数为 {search.method} 法的 Fs 最小值。</div>"
        )

    html = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title><style>{_CSS}</style></head>
<body><div class="wrap">
<h1>{title}</h1>
<div class="meta">方法 {lem.method} ｜ 生成于 {_dt.datetime.now():%Y-%m-%d %H:%M} ｜ slopestab {__version__}</div>
{geometry}
<div class="cards">{cards}</div>
{img_tag}
<h2>计算摘要</h2>
<table><tbody>{summary_rows}</tbody></table>
<h2>平衡校核（相对总重的百分比）</h2>
<table><thead><tr><th>量</th><th>残差 %</th></tr></thead><tbody>{residual_rows}</tbody></table>
<h2>条分结果</h2>
<table><thead><tr><th>条号</th><th>x 左 (m)</th><th>宽 b (m)</th><th>高 h (m)</th>
<th>重 W (kN/m)</th><th>孔压 u (kPa)</th><th>sin α</th><th>cos α</th></tr></thead>
<tbody>{slice_rows}</tbody></table>
{search_note}
<footer>slopestab 为开源极限平衡法实现，计算结果需由注册岩土工程师复核后方可用于设计。</footer>
</div></body></html>
"""
    path = Path(path)
    path.write_text(html, encoding="utf-8")
    return path


def _unpack(result: Any) -> tuple[Any, Any]:
    from ..search.circular import SearchResult

    if isinstance(result, SearchResult):
        return result.critical, result
    if hasattr(result, "critical"):
        return result.critical, result
    return result, None


def _cards(lem: Any, slices: Any, fs: float, fs_class: str) -> str:
    lam = lem.info.get("lambda")
    extra = (
        f"<div class='card'><div class='k'>λ</div><div class='v'>{float(lam):.3f}</div></div>"
        if isinstance(lam, (int, float))
        else ""
    )
    return (
        f"<div class='card'><div class='k'>安全系数 Fs</div>"
        f"<div class='v {fs_class}'>{fs:.3f}</div></div>"
        f"<div class='card'><div class='k'>条分数</div><div class='v'>{len(slices)}</div></div>"
        f"<div class='card'><div class='k'>滑体总重</div>"
        f"<div class='v'>{slices.total_weight:.0f} <small>kN/m</small></div></div>"
        f"<div class='card'><div class='k'>滑面跨度</div>"
        f"<div class='v'>{slices.span[1] - slices.span[0]:.1f} <small>m</small></div></div>"
        f"{extra}"
    )


def _summary_rows(lem: Any) -> str:
    rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in lem.summary().items()
        if k != "residual_percent"
    )
    residuals = lem.summary().get("residual_percent", {})
    res = "".join(f"<tr><td>residual {k}</td><td>{v}</td></tr>" for k, v in residuals.items())
    return rows + res


def _slice_rows(slices: Any, max_rows: int) -> str:
    rows = []
    for i in range(min(len(slices), max_rows)):
        rows.append(
            "<tr>"
            f"<td>{i + 1}</td><td>{slices.x_left[i]:.2f}</td><td>{slices.b[i]:.2f}</td>"
            f"<td>{slices.h[i]:.2f}</td><td>{slices.w[i]:.1f}</td><td>{slices.u[i]:.1f}</td>"
            f"<td>{slices.sin_a[i]:+.3f}</td><td>{slices.cos_a[i]:.3f}</td>"
            "</tr>"
        )
    if len(slices) > max_rows:
        rows.append(f"<tr><td colspan='8'>… 其余 {len(slices) - max_rows} 条从略</td></tr>")
    return "".join(rows)


def _inline_plot(result: Any, lem: Any) -> str:
    try:
        from .plot import plot
    except Exception:  # pragma: no cover - matplotlib missing
        return ""
    try:
        fig = plot(result, show_slices=True)
        import io as _io

        buf = _io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        import matplotlib.pyplot as plt

        plt.close(fig)
        buf.seek(0)
        b64 = base64.b64encode(buf.read()).decode("ascii")
        return f"<img src='data:image/png;base64,{b64}' alt='slope section'>"
    except Exception:  # pragma: no cover - plotting is best effort
        return ""
