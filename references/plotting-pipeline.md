# 出图管线详解、质检清单与踩坑点

## 一、脚本骨架（照此分七段，段间用 `# ------ n. 标题 ------` 分隔）

```
1. 实验数据        x、y 数组 + 待拟合参数列表 + 配色 + 署名常量
2. 拟合与不确定度   自由度 ν、t 系数、fit 函数（返回全部统计量）
3. 数值结果输出     print + 汇总文本
4. 绘图风格与字体   rcParams（中文字体、unicode_minus、字号、线宽）
5. 每个公式单独出图  循环 MODELS → scatter + plot + 标签 + style + legend + text 框 → savefig
6. 多公式对比图     所有曲线叠加 + 等宽字体参数表
7. 结果文本落盘
```

每段只做一件事，换实验时只改第 1 段（数据/公式/标签）即可复跑。

## 二、出图固定套路

```python
fig, ax = plt.subplots(figsize=(9.4, 6.4), constrained_layout=True)
fig.suptitle(f"...（k = {k}）", fontsize=17)

ax.scatter(lam, n_obs, s=62, marker="o",
           facecolors="white", edgecolors=color, linewidths=1.8,
           label="实验数据", zorder=3)
ax.plot(lam_line, a + b * lam_line ** k,
        color="#d62728", linewidth=2.0, label="拟合曲线", zorder=2)

ax.set_xlabel(r"波长 $\lambda$ / $\mathrm{\AA}$", fontsize=13, labelpad=8)
ax.set_ylabel("折射率 n", fontsize=13, labelpad=8)
y_pad = (n_obs.max() - n_obs.min()) * 0.10
ax.set_ylim(n_obs.min() - y_pad, n_obs.max() + y_pad * 1.6)

ax.grid(True, which="major", linestyle="--", linewidth=0.65, alpha=0.45)
ax.tick_params(direction="in", length=5, width=1)
ax.legend(loc="lower left", frameon=True, framealpha=0.95, fontsize=11)
ax.text(0.97, 0.97, result_text, transform=ax.transAxes,
        ha="right", va="top", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="0.75", alpha=0.92))
fig.savefig(stem + ".png", dpi=600, bbox_inches="tight")
fig.savefig(stem + ".pdf", bbox_inches="tight")
```

## 三、"空白象限"摆放规则（遮挡问题的根治办法）

**先判断数据在坐标系里的走向，再选空白处放文字框和图例**，不要凭感觉放。

| 数据走向 | 空闲象限 | 典型场景 |
|---|---|---|
| 左上 → 右下（单调递减，如 n～λ） | **右上**、左下 | 结果框放右上（`ha="right"`），图例放左下 |
| 左下 → 右上（单调递增，如 n～λ^k，k<0 时 X=λ^k） | **左上**、右下 | 说明框放左上，图例放右下 |

两个图元（结果框、图例）**必须分别占用不同象限**，同侧必然相撞。
出图后逐张目视确认：文字框下缘应高于/低于曲线，不与任何点重叠。

## 四、质检清单（交付前逐项打勾）

1. 运行无 `Traceback`；`grep -i glyph` 无命中（无方框缺字）。
2. 所有统计量（r、S_y、a、b、s_a、s_b、u_a、u_b）与用户提供的软件数值逐位一致；
   若用户未提供，则口径（ν、t、P）已在脚本与交付说明中写明。
3. 每个公式一张图，无多余面板/子图。
4. 结果框、图例均落在空白象限，不遮挡曲线与数据点（**已看图确认，非推断**）。
5. 坐标轴有物理量名称与单位；单位符号渲染正常（如 Å）。
6. 拟合曲线为连续平滑曲线，两端略超出数据范围。
7. PNG（600 dpi）与 PDF 均已生成，路径在脚本所在目录。
8. `拟合结果汇总.txt` 含原始数据、全部统计量、残差与结论。
9. `scripts/` 内代码文件首行 `# Designed by TiantianYe`，头部有适配指南。
10. 记忆日志已更新（本次沉淀/改动）。

## 五、matplotlib 踩坑点

| 现象 | 根因 | 解法 |
|---|---|---|
| 中文显示为方框 | 默认字体无 CJK 字形 | 设 `font.sans-serif` 回退列表（顺序：Hiragino Sans GB → … → Arial Unicode MS） |
| `Glyph 197 missing` 告警 | 中文字体不含 `Å`(U+00C5) | 单位改用 mathtext `r"$\mathrm{\AA}$"` |
| `Font 'default' does not have a glyph for '\u2212'` | mathtext/字体缺减号 | `axes.unicode_minus = False`；本告警在本案例中不影响输出，但须目视确认 |
| 拟合线看起来是折线 | 只连接了数据点 | 用 600 点 `linspace` 画连续曲线 |
| 文字框压住曲线 | 凭感觉摆放 | 按 §三 空白象限规则，且两图元分占不同象限 |
| 无图形界面时报错/挂起 | 后端问题 | 本机运行可 `MPLBACKEND=Agg python3 xxx.py` |
| 从别的目录运行，图片写错位置 | 用了相对路径 | `os.path.dirname(os.path.abspath(__file__))` |

## 六、单位与量级注意

- 波长常用 **Å**（1 Å = 0.1 nm）；若用 nm，系数 b 会差 10^k 倍，**必须与用户数据单位一致**。
- `b` 的量纲是 `y·x^(−k)`；换单位时 b 必须重算，不能沿用。
- 报告里若引用 `b` 的数值，必须同时给出 x 的单位。
