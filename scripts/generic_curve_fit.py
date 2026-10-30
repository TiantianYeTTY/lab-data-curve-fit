# Designed by TiantianYe
# -*- coding: utf-8 -*-
"""
通用实验数据函数拟合图生成模板

适用：任何能通过变量代换化为 y = a + b·X 的实验公式
      （y=a+bx、y=a+b·x^k、y=a+b·ln x、y=a·e^(bx)、y=a·x^b 等）

统计量口径与大学物理实验教材 / 实验数据处理软件一致：
    r、S_y = sqrt(Σ(y-ŷ)²/(n-2))、s_a、s_b、u = t·s（P = 0.95，ν = n-2）

【复用适配指南】换实验时只改第 0 段的四处配置，其余不动：
1. DATA      —— x、y 数组。单位必须自洽：b 的量纲 = y·x^(−k)，换单位 b 要重算
2. MODELS    —— 每个模型四项：
        tag       输出文件名后缀（如 "k-1"）
        expr      图例与结果框里显示的公式文本
        x_trans   线性化变换 X = f(x)，直线拟合用
        y_trans   y 的变换及其逆变换 (f, f⁻¹)，不需要就写 (lambda v: v, lambda v: v)
        color     数据点描边色
3. LABELS    —— 横轴、纵轴、图标题。含单位；Å 之类缺字形符号用 mathtext
4. STUDENT   —— 署名（学号+姓名），显示在结果框首行

不要改：统计量口径（ν = n−2、t 表、u = t·s）、出图套路（600 点连续曲线、
空心圆 s=62）、结果框九项（与实验数据处理软件显示保持一致）。
"""

import os

import numpy as np
import matplotlib.pyplot as plt

# ==================== 0. 配置段（换实验只改这里） ====================

# --- 1) 数据 ---
x = np.array([5790.7, 5460.7, 4916.1, 4358.3, 4046.6], dtype=float)   # 波长 λ / Å
y = np.array([1.6516, 1.6544, 1.6636, 1.6728, 1.6882], dtype=float)  # 折射率 n

# --- 2) 待拟合模型 ---
# y = a + b·x^k，本例拟合三个 k
MODELS = [
    {"tag": "k-1", "expr": "n = a + b·λ^(-1)",
     "x_trans": lambda v: v ** -1, "y_trans": (lambda v: v, lambda v: v),
     "color": "#1769aa"},
    {"tag": "k-2", "expr": "n = a + b·λ^(-2)",
     "x_trans": lambda v: v ** -2, "y_trans": (lambda v: v, lambda v: v),
     "color": "#2ca02c"},
    {"tag": "k-4", "expr": "n = a + b·λ^(-4)",
     "x_trans": lambda v: v ** -4, "y_trans": (lambda v: v, lambda v: v),
     "color": "#d62728"},
]

# 其它常见模型（按需替换上面 MODELS）：
#   指数 y = a·e^(bx)：  "x_trans": lambda v: v,  "y_trans": (np.log, np.exp)
#   幂律 y = a·x^b  ：  "x_trans": np.log,         "y_trans": (np.log, np.exp)
#   对数 y = a + b·ln x："x_trans": np.log,         "y_trans": (lambda v: v, lambda v: v)
#   线性 y = a + bx  ：  "x_trans": lambda v: v,   "y_trans": (lambda v: v, lambda v: v)
# ⚠️ 对 y 取变换时，r 与 S_y 都是变换后空间的值，交付说明里要写明"Y = ln y"。

# --- 3) 标签 ---
LABELS = {
    "xlabel": r"波长 $\lambda$ / $\mathrm{\AA}$",   # 缺字形符号请用 mathtext
    "ylabel": "折射率 n",
    "title": "棱镜折射率色散拟合图",
}

# --- 4) 署名 ---
STUDENT = "学号+姓名"

# ==================== 1. 拟合与不确定度 ====================

n_pts = len(x)
nu = n_pts - 2                          # 自由度
# t 分布置信系数（P = 0.95）；ν 超出下表范围需另查教材 t 表
T_TABLE = {1: 12.7, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57,
           6: 2.45, 7: 2.36, 8: 2.31, 9: 2.26, 10: 2.23}
T_FACTOR = T_TABLE[nu]


def fit_linearized(x, y, x_trans, y_trans):
    """把 (x, y) 变换到 (X, Y) 后作一元线性最小二乘，返回全部统计量。"""
    y_fwd, y_inv = y_trans
    X = x_trans(x)
    Y = y_fwd(y)

    b, a = np.polyfit(X, Y, 1)          # 先斜率 b，后截距 a
    Y_fit = a + b * X
    r = np.corrcoef(X, Y)[0, 1]

    residual = Y - Y_fit
    s_y = np.sqrt(np.sum(residual ** 2) / (n_pts - 2))

    s_xx = np.sum((X - X.mean()) ** 2)
    s_a = s_y * np.sqrt(1.0 / n_pts + X.mean() ** 2 / s_xx)
    s_b = s_y / np.sqrt(s_xx)

    u_a = T_FACTOR * s_a
    u_b = T_FACTOR * s_b

    return dict(a=a, b=b, X=X, Y=Y, Y_fit=Y_fit, r=r,
                s_y=s_y, s_a=s_a, s_b=s_b, u_a=u_a, u_b=u_b,
                residual=residual, y_inv=y_inv)


def fmt_e(x):
    """科学计数、6 位小数、去尾随 0（与实验数据处理软件显示一致）。"""
    mant, exp = f"{x:.6E}".split("E")
    return f"{mant.rstrip('0').rstrip('.')}E{exp}"


results = [fit_linearized(x, y, m["x_trans"], m["y_trans"]) for m in MODELS]

# ==================== 2. 控制台输出 ====================

print(f"数据点数 n = {n_pts}，自由度 ν = {nu}，置信概率 P = 0.95，t = {T_FACTOR}")
for m, R in zip(MODELS, results):
    print(f"\n———— {m['tag']}：{m['expr']} ————")
    print(f"  相关系数 r            {R['r']:.7f}")
    print(f"  Y的标准偏差 S_y       {fmt_e(R['s_y'])}")
    print(f"  系数 a                {R['a']:.6f}")
    print(f"  a 的标准偏差          {fmt_e(R['s_a'])}")
    print(f"  a 的 A 类不确定度分量 {fmt_e(R['u_a'])}")
    print(f"  系数 b                {fmt_e(R['b'])}")
    print(f"  b 的标准偏差          {fmt_e(R['s_b'])}")
    print(f"  b 的不确定度          {fmt_e(R['u_b'])}")
    print("  各点残差              "
          + ", ".join(f"{v:+.5f}" for v in R["residual"]))

best = max(zip(MODELS, results), key=lambda mr: abs(mr[1]["r"]))
print(f"\n结论：|r| 最大者为 {best[0]['tag']}，该形式对本组数据拟合最优。")

# ==================== 3. 绘图风格 ====================

# 中文字体自动回退，避免标题/标签显示为方框
plt.rcParams["font.sans-serif"] = [
    "Hiragino Sans GB", "Microsoft YaHei", "SimHei", "PingFang SC",
    "Noto Sans CJK SC", "Source Han Sans SC", "Arial Unicode MS"
]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 12
plt.rcParams["axes.linewidth"] = 1.1

script_dir = os.path.dirname(os.path.abspath(__file__))
x_line = np.linspace(x.min() * 0.985, x.max() * 1.015, 600)


def style(ax):
    ax.grid(True, which="major", linestyle="--", linewidth=0.65, alpha=0.45)
    ax.tick_params(direction="in", length=5, width=1)


def result_text(m, R):
    """图内结果框文本，项名与实验数据处理软件保持一致。"""
    return (
        f"{STUDENT}\n"
        f"函数关系式：{m['expr']}\n"
        "拟合结果：\n"
        f"相关系数 r = {R['r']:.7f}\n"
        f"Y的标准偏差 = {fmt_e(R['s_y'])}\n"
        f"系数 a = {R['a']:.6f}\n"
        f"a 的标准偏差 = {fmt_e(R['s_a'])}\n"
        f"a 的 A 类不确定度分量 = {fmt_e(R['u_a'])}\n"
        f"系数 b = {fmt_e(R['b'])}\n"
        f"b 的标准偏差 = {fmt_e(R['s_b'])}\n"
        f"b 的不确定度 = {fmt_e(R['u_b'])}"
    )


# ==================== 4. 每个模型单独一张图 ====================

for m, R in zip(MODELS, results):
    fig, ax = plt.subplots(figsize=(9.4, 6.4), constrained_layout=True)
    fig.suptitle(f"{LABELS['title']}（{m['tag']}）", fontsize=17)

    ax.scatter(
        x, y, s=62, marker="o",
        facecolors="white", edgecolors=m["color"], linewidths=1.8,
        label="实验数据", zorder=3
    )
    # 原始变量空间的平滑拟合曲线：y = y_inv(a + b·x_trans(x))
    ax.plot(
        x_line, R["y_inv"](R["a"] + R["b"] * m["x_trans"](x_line)),
        color="#d62728", linewidth=2.0,
        label=f"拟合曲线 {m['expr']}", zorder=2
    )

    ax.set_xlabel(LABELS["xlabel"], fontsize=13, labelpad=8)
    ax.set_ylabel(LABELS["ylabel"], fontsize=13, labelpad=8)
    y_pad = (y.max() - y.min()) * 0.10
    ax.set_ylim(y.min() - y_pad, y.max() + y_pad * 1.6)
    style(ax)
    ax.legend(loc="lower left", frameon=True, framealpha=0.95, fontsize=11)
    ax.text(
        0.97, 0.97, result_text(m, R),
        transform=ax.transAxes, ha="right", va="top", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="0.75", alpha=0.92)
    )

    stem = os.path.join(script_dir, f"fit_{m['tag']}")
    fig.savefig(stem + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    print(f"已保存：{stem}.png / {stem}.pdf")
    plt.close(fig)

# ==================== 5. 多模型对比图 ====================

fig, ax = plt.subplots(figsize=(9.4, 6.4), constrained_layout=True)
ax.scatter(
    x, y, s=62, marker="o",
    facecolors="white", edgecolors="#333333", linewidths=1.8,
    label="实验数据", zorder=4
)
for m, R in zip(MODELS, results):
    ax.plot(
        x_line, R["y_inv"](R["a"] + R["b"] * m["x_trans"](x_line)),
        color=m["color"], linewidth=2.0,
        label=f"{m['tag']}：{m['expr']}，r = {R['r']:.7f}", zorder=3
    )
ax.set_title(f"{LABELS['title']}多模型对比", fontsize=17, pad=14)
ax.set_xlabel(LABELS["xlabel"], fontsize=13, labelpad=8)
ax.set_ylabel(LABELS["ylabel"], fontsize=13, labelpad=8)
y_pad = (y.max() - y.min()) * 0.10
ax.set_ylim(y.min() - y_pad, y.max() + y_pad * 1.6)
style(ax)
ax.legend(loc="lower left", frameon=True, framealpha=0.95, fontsize=11)

table = ["  tag      r            S_y           a           b"]
for m, R in zip(MODELS, results):
    table.append(
        f"  {m['tag']:<5}  {R['r']:.7f}   {R['s_y']:.3E}   "
        f"{R['a']:.6f}   {R['b']:.4E}"
    )
ax.text(
    0.97, 0.97, "\n".join(table),
    transform=ax.transAxes, ha="right", va="top", fontsize=10.5,
    family="monospace",
    bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
              edgecolor="0.75", alpha=0.92)
)

cmp_stem = os.path.join(script_dir, "fit_compare")
fig.savefig(cmp_stem + ".png", dpi=600, bbox_inches="tight")
fig.savefig(cmp_stem + ".pdf", bbox_inches="tight")
print(f"已保存：{cmp_stem}.png / {cmp_stem}.pdf")
plt.close(fig)

plt.show()
