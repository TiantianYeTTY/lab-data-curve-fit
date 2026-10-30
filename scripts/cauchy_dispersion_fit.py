# Designed by TiantianYe
# -*- coding: utf-8 -*-
"""
分光计测棱镜玻璃折射率 —— 柯西色散公式 n = a + b·λ^k 拟合（k = -1, -2, -4）

数据来源：实验记录 5 组测量数据
    λ：汞灯可见光谱线波长，单位 Å（1 Å = 0.1 nm）
    n：棱镜材料对应波长下的折射率

拟合方法（与实验数据处理软件完全一致）：
    令 X = λ^k，则 n = a + b·X 为线性关系，用一元线性最小二乘求解。
    相关系数 r、Y 的标准偏差 S_y = sqrt(Σ(n-ŷ)²/(n-2))、
    a、b 的标准偏差 s_a、s_b，以及 P=0.95 下的 A 类不确定度 u = t·s
    （本实验数据点数 n_pts = 5 → 自由度 ν = 3 → t = 3.18）。

运行：先 python3 -m pip install numpy matplotlib，再运行本脚本
    python3 cauchy_dispersion_fit.py

输出（均保存在脚本所在目录，600 dpi PNG + 矢量 PDF）：
    dispersion_k-1_n_lambda.png / .pdf      k = -1 的 n～λ 拟合曲线图
    dispersion_k-2_n_lambda.png / .pdf      k = -2 的 n～λ 拟合曲线图
    dispersion_k-4_n_lambda.png / .pdf      k = -4 的 n～λ 拟合曲线图
    dispersion_compare.png / .pdf           三条拟合曲线对比图
    拟合结果汇总.txt                         全部数值结果

【复用适配指南】换实验时必须改的四处：
1. 数据段：lam / n_obs 两个数组（单位要与系数 b 的量纲自洽）
2. K_LIST：待拟合的指数集合；COLORS：每个 k 的配色
3. 图例与标题文案、横纵坐标标签（含单位）
4. STUDENT：署名；stem：输出文件名前缀
不要改：统计量口径（ν = n−2、t 表、u = t·s）、出图套路（600 点连续曲线、
空心圆 s=62）、结果框的九项显示（与实验数据处理软件保持一致）。
"""

import os

import numpy as np
import matplotlib.pyplot as plt

# -------------------- 1. 实验数据 --------------------
# 汞灯可见光谱线波长 λ（单位 Å）与对应折射率 n
lam = np.array([5790.7, 5460.7, 4916.1, 4358.3, 4046.6], dtype=float)
n_obs = np.array([1.6516, 1.6544, 1.6636, 1.6728, 1.6882], dtype=float)

# 待拟合的三个指数
K_LIST = [-1, -2, -4]

# 每个 k 对应的配色（蓝 / 绿 / 红），区分三张图与对比图
COLORS = {-1: "#1769aa", -2: "#2ca02c", -4: "#d62728"}

STUDENT = "学号+姓名"        # 报告署名，改成你自己的学号与姓名

# -------------------- 2. 拟合与不确定度计算 --------------------
n_pts = len(lam)
nu = n_pts - 2                      # 自由度 = 3
# t 分布置信系数（P = 0.95）
T_FACTOR = {1: 12.7, 2: 4.30, 3: 3.18, 4: 2.78, 5: 2.57,
            6: 2.45, 7: 2.36, 8: 2.31, 9: 2.26, 10: 2.23}[nu]


def fit_power(lam, y, k):
    """对 y = a + b·λ^k 作线性化最小二乘拟合，返回全部统计量。"""
    X = lam ** k                              # 线性化自变量
    b, a = np.polyfit(X, y, 1)                # 斜率 b，截距 a
    y_fit = a + b * X
    r = np.corrcoef(X, y)[0, 1]               # 相关系数
    residual = y - y_fit
    s_y = np.sqrt(np.sum(residual ** 2) / (n_pts - 2))   # Y 的标准偏差

    s_xx = np.sum((X - X.mean()) ** 2)
    s_a = s_y * np.sqrt(1.0 / n_pts + X.mean() ** 2 / s_xx)  # a 的标准偏差
    s_b = s_y / np.sqrt(s_xx)                                # b 的标准偏差
    u_a = T_FACTOR * s_a                      # a 的 A 类不确定度分量
    u_b = T_FACTOR * s_b                      # b 的不确定度

    return dict(k=k, X=X, a=a, b=b, y_fit=y_fit, r=r, s_y=s_y,
                s_a=s_a, s_b=s_b, u_a=u_a, u_b=u_b, residual=residual)


results = {k: fit_power(lam, n_obs, k) for k in K_LIST}

# 科学计数格式：与处理软件显示一致（6 位小数、去掉尾随 0）
def fmt_e(x):
    mant, exp = f"{x:.6E}".split("E")
    return f"{mant.rstrip('0').rstrip('.')}E{exp}"


# -------------------- 3. 数值结果输出 --------------------
lines = []
lines.append("柯西色散公式 n = a + b·λ^k 拟合结果汇总")
lines.append(f"数据点数 n = {n_pts}，自由度 ν = {nu}，置信概率 P = 0.95，t = {T_FACTOR}")
lines.append("")
lines.append("原始数据：")
lines.append("  序号    λ / Å        n")
for i, (L, nn) in enumerate(zip(lam, n_obs), start=1):
    lines.append(f"  {i:>2d}   {L:>9.1f}   {nn:.4f}")
lines.append("")

for k in K_LIST:
    R = results[k]
    lines.append(f"———— k = {k}，函数关系式 n = a + b·λ^({k}) ————")
    lines.append(f"  线性关系式            n = a + b·X，X = λ^({k})")
    lines.append(f"  相关系数 r            {R['r']:.7f}")
    lines.append(f"  Y 的标准偏差 S_y      {fmt_e(R['s_y'])}")
    lines.append(f"  系数 a                {R['a']:.6f}")
    lines.append(f"  a 的标准偏差 s_a      {fmt_e(R['s_a'])}")
    lines.append(f"  a 的 A 类不确定度 u_a {fmt_e(R['u_a'])}")
    lines.append(f"  系数 b                {fmt_e(R['b'])}")
    lines.append(f"  b 的标准偏差 s_b      {fmt_e(R['s_b'])}")
    lines.append(f"  b 的不确定度 u_b      {fmt_e(R['u_b'])}")
    lines.append(f"  各点残差 n - n_fit    "
                 + ", ".join(f"{v:+.5f}" for v in R["residual"]))
    lines.append("")

best_k = max(K_LIST, key=lambda k: abs(results[k]["r"]))
lines.append(f"结论：|r| 最大者为 k = {best_k}，该形式对本组数据拟合最优。")

summary_txt = "\n".join(lines)
print(summary_txt)

# -------------------- 4. 绘图风格与字体 --------------------
# 按系统可用字体自动选择常见中文字体，避免中文标题显示为方框
plt.rcParams["font.sans-serif"] = [
    "Hiragino Sans GB", "Microsoft YaHei", "SimHei", "PingFang SC",
    "Noto Sans CJK SC", "Source Han Sans SC", "Arial Unicode MS"
]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 12
plt.rcParams["axes.linewidth"] = 1.1

script_dir = os.path.dirname(os.path.abspath(__file__))


def result_text(k, R):
    """生成图内标注文本（与实验数据处理软件显示项一致）。"""
    return (
        f"{STUDENT}\n"
        f"函数关系式：n = a + b·λ^({k})\n"
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


def style(ax):
    ax.grid(True, which="major", linestyle="--", linewidth=0.65, alpha=0.45)
    ax.tick_params(direction="in", length=5, width=1)


# -------------------- 5. 每个 k 单独出一张 n～λ 拟合曲线图 --------------------
lam_line = np.linspace(lam.min() * 0.985, lam.max() * 1.015, 600)

for k in K_LIST:
    R = results[k]
    color = COLORS[k]

    fig, ax = plt.subplots(figsize=(9.4, 6.4), constrained_layout=True)
    fig.suptitle(f"棱镜折射率色散拟合图（k = {k}）", fontsize=17)

    # 原始变量空间 n ～ λ
    ax.scatter(
        lam, n_obs, s=62, marker="o",
        facecolors="white", edgecolors=color, linewidths=1.8,
        label="实验数据", zorder=3
    )
    ax.plot(
        lam_line, R["a"] + R["b"] * lam_line ** k,
        color="#d62728", linewidth=2.0,
        label=f"拟合曲线 n = a + b·λ^({k})", zorder=2
    )
    ax.set_xlabel(r"波长 $\lambda$ / $\mathrm{\AA}$", fontsize=13, labelpad=8)
    ax.set_ylabel("折射率 n", fontsize=13, labelpad=8)
    y_pad = (n_obs.max() - n_obs.min()) * 0.10
    ax.set_ylim(n_obs.min() - y_pad, n_obs.max() + y_pad * 1.6)
    style(ax)
    ax.legend(loc="lower left", frameon=True, framealpha=0.95, fontsize=11)
    ax.text(
        0.97, 0.97, result_text(k, R),
        transform=ax.transAxes, ha="right", va="top", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="0.75", alpha=0.92)
    )

    stem = os.path.join(script_dir, f"dispersion_k{k}_n_lambda")
    fig.savefig(stem + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    print(f"已保存：{stem}.png")
    print(f"已保存：{stem}.pdf")
    plt.close(fig)

# -------------------- 6. 三条拟合曲线对比图 --------------------
fig, ax = plt.subplots(figsize=(9.4, 6.4), constrained_layout=True)
ax.scatter(
    lam, n_obs, s=62, marker="o",
    facecolors="white", edgecolors="#333333", linewidths=1.8,
    label="实验数据", zorder=4
)
for k in K_LIST:
    R = results[k]
    ax.plot(
        lam_line, R["a"] + R["b"] * lam_line ** k,
        color=COLORS[k], linewidth=2.0,
        label=f"k = {k}：n = a + b·λ^({k})，r = {R['r']:.7f}", zorder=3
    )
ax.set_title("三种 k 值的色散拟合曲线对比", fontsize=17, pad=14)
ax.set_xlabel(r"波长 $\lambda$ / $\mathrm{\AA}$", fontsize=13, labelpad=8)
ax.set_ylabel("折射率 n", fontsize=13, labelpad=8)
y_pad = (n_obs.max() - n_obs.min()) * 0.10
ax.set_ylim(n_obs.min() - y_pad, n_obs.max() + y_pad * 1.6)
style(ax)
ax.legend(loc="lower left", frameon=True, framealpha=0.95, fontsize=11)

table_lines = ["   k       r            S_y           a           b"]
for k in K_LIST:
    R = results[k]
    table_lines.append(
        f"  {k:>2d}   {R['r']:.7f}   {R['s_y']:.3E}   {R['a']:.6f}   {R['b']:.4E}"
    )
ax.text(
    0.97, 0.97, "\n".join(table_lines),
    transform=ax.transAxes, ha="right", va="top", fontsize=10.5,
    family="monospace",
    bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
              edgecolor="0.75", alpha=0.92)
)

cmp_stem = os.path.join(script_dir, "dispersion_compare")
fig.savefig(cmp_stem + ".png", dpi=600, bbox_inches="tight")
fig.savefig(cmp_stem + ".pdf", bbox_inches="tight")
print(f"已保存：{cmp_stem}.png")
print(f"已保存：{cmp_stem}.pdf")
plt.close(fig)

# -------------------- 7. 结果文本落盘 --------------------
txt_path = os.path.join(script_dir, "拟合结果汇总.txt")
with open(txt_path, "w", encoding="utf-8") as f:
    f.write(summary_txt + "\n")
print(f"已保存：{txt_path}")

plt.show()
