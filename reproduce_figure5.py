"""
================================================================================
Figure 5A / 5B / 5D 复现 —— 单文件最终版
================================================================================
严格按论文作者官方代码复现：
  仓库   : https://github.com/ZhuangLab/SpatialBrainAgingCell22
  notebook: notebooks/merfish_spatial_celltype_org.ipynb (CELL 41-43)

作者原始代码（逐条照做，未做任何简化）：
  CELL41:
    sc.tl.score_genes(adata_annot, gene_list=[...'Ctss'],
                      score_name='activate_micro', use_raw=False)
    sc.tl.score_genes(adata_annot, gene_list=[...'Hif3a'],
                      score_name='activate_astro', use_raw=False)
  CELL42:
    adata_micro = adata_annot[adata_annot.obs.cell_type=="Micro"]
    adata_annot.obs.activate_micro -= np.mean(
        adata_micro[adata_micro.obs.age=='4wk'].obs.activate_micro)
  CELL43:
    adata_astro = adata_annot[adata_annot.obs.cell_type=="Astro"]
    adata_annot.obs.activate_astro -= np.mean(
        adata_astro[adata_astro.obs.age=='4wk'].obs.activate_astro)

用法:
    python reproduce_figure5.py
输入:
    C:\\Users\\yyj\\Downloads\\c93d78c2-ee17-4504-8d1c-17cf093ad7b5.h5ad
输出:
    ./figure5_out/5A_activation_violin.png
    ./figure5_out/5B_cluster_violin.png
    ./figure5_out/5D_region_boxplot.png
    ./figure5_out/summary_tables.txt
    ./figure5_out/astro_micro_scores.h5ad
================================================================================
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

# ============================== 路径配置 ==============================
IN = r"C:\Users\yyj\Downloads\c93d78c2-ee17-4504-8d1c-17cf093ad7b5.h5ad"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figure5_out")
os.makedirs(OUT, exist_ok=True)

# ============================== 基因集（作者原文）==============================
GENES_MICRO = ['B2m', 'Trem2', 'Ccl2', 'Apoe', 'Axl', 'Itgax',
               'Cd9', 'C1qa', 'C1qc', 'Lyz2', 'Ctss']
GENES_ASTRO = ['C4b', 'C3', 'Serpina3n', 'Cxcl10', 'Gfap', 'Vim',
               'Il18', 'Hif3a']

# ============================== 绘图参数 ==============================
# 配色来自论文原图 PDF 400dpi 采样（非目测）
C_JUV, C_YA, C_OLD = "#739ADA", "#D7C2D6", "#E08E8F"   # 蓝 / 淡紫 / 鲑红
C_ASTRO = "#367F5A"                                    # Astro 深绿
C_MICRO = ["#892B5D", "#B42E77", "#E23290"]             # Micro-1/2/3 品红渐变
AGE_LABELS = ["Juvenile", "Young Adult", "Old"]

REGION_ORDER = ["Pia", "L2/3", "L5", "L6",
                "Corpus Callosum", "Olf", "Striatum", "Ventricle"]
TISSUE_MAP = {
    "pia mater": "Pia", "cortical layer II/III": "L2/3",
    "cortical layer V": "L5", "cortical layer VI": "L6",
    "corpus callosum": "Corpus Callosum", "olfactory region": "Olf",
    "striatum": "Striatum", "brain ventricle": "Ventricle",
}
# 纵轴留白：下限 -1.3 让小提琴底端悬空、不碰轴线
Y_LIM = {"astro": (-1.3, 3.1), "micro": (-1.3, 2.1)}


# =============================================================================
# 第一部分：按作者代码计算激活分数
# =============================================================================
def compute_scores(in_path):
    print("[1/3] 读取数据 ...", flush=True)
    src_ad = ad.read_h5ad(in_path)

    # 用 raw（原始计数）矩阵，var_names 换成基因符号
    raw = src_ad.raw.to_adata() if src_ad.raw is not None else src_ad
    raw.var_names = pd.Index(raw.var["feature_name"].astype(str))
    obs_names = src_ad.obs_names.copy()

    adata_annot = ad.AnnData(X=raw.X.copy(), obs=src_ad.obs.copy(),
                             var=src_ad.var.copy())
    adata_annot.var_names = pd.Index(
        adata_annot.var["feature_name"].astype(str))
    adata_annot.var["feature_name"] = adata_annot.var_names.astype(str)
    adata_annot.var.index.name = None
    adata_annot.obs_names = obs_names
    del raw, src_ad
    print(f"      shape = {adata_annot.shape}", flush=True)

    # ---- 预处理 ----
    sc.pp.normalize_total(adata_annot, target_sum=1e4)
    sc.pp.log1p(adata_annot)
    print("[2/3] normalize_total(1e4) + log1p 完成", flush=True)

    # ---- CELL41：算分。use_raw=False，不设 random_state/ctrl_size（用scanpy 默认）----
    sc.tl.score_genes(adata_annot, gene_list=GENES_MICRO,
                      score_name='activate_micro', use_raw=False)
    sc.tl.score_genes(adata_annot, gene_list=GENES_ASTRO,
                      score_name='activate_astro', use_raw=False)
    print("      score_genes 完成（4wk 归零前）", flush=True)

    # ---- CELL42/43：归零 ----
    # 注意：本数据集 obs.cell_type 用全称（astrocyte / microglial cell），
    # 而作者代码写的是 "Astro" / "Micro"。需语义映射，否则取到空集→ 全NaN。
    CT_MAP = {"microglial cell": "Micro", "astrocyte": "Astro"}
    ct = adata_annot.obs.cell_type.astype(str)
    adata_annot.obs["cell_type_author"] = ct.map(CT_MAP).fillna(ct)

    adata_micro = adata_annot[adata_annot.obs.cell_type_author == "Micro"]
    bm = np.mean(adata_micro[adata_micro.obs.age == '4wk'].obs.activate_micro)
    adata_annot.obs.activate_micro = adata_annot.obs.activate_micro - bm

    adata_astro = adata_annot[adata_annot.obs.cell_type_author == "Astro"]
    ba = np.mean(adata_astro[adata_astro.obs.age == '4wk'].obs.activate_astro)
    adata_annot.obs.activate_astro = adata_annot.obs.activate_astro - ba
    print(f"[3/3] 4wk 归零完成: micro -= {bm:.6f} | astro -= {ba:.6f}",
          flush=True)

    adata_annot.obs["activate_micro"] = \
        adata_annot.obs.activate_micro.astype("float32")
    adata_annot.obs["activate_astro"] = \
        adata_annot.obs.activate_astro.astype("float32")
    return adata_annot


# =============================================================================
# 第二部分：绘图
# =============================================================================
def style(ax, ylab, ylim):
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylab, fontsize=10)
    ax.set_xlabel("")
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)
    ax.tick_params(axis="x", length=0, labelsize=9.5)
    ax.tick_params(axis="y", length=3, labelsize=8.5)


def violins(ax, groups, colors, tick_labels):
    """sns.violinplot(cut=0, bw_adjust=1) + 黑点收在中央窄竖带 (jitter=0.1)"""
    df = pd.DataFrame({
        "grp": np.concatenate([[tick_labels[i]] * len(g)
                               for i, g in enumerate(groups)]),
        "val": np.concatenate([g for g in groups if len(g)]),
    })
    order = [t for t, g in zip(tick_labels, groups) if len(g)]
    sns.violinplot(data=df, x="grp", y="val", order=order, palette=colors,
                   linewidth=0.5, cut=0, bw_adjust=1, density_norm="width",
                   ax=ax, legend=False)
    ax.set_xlabel("")
    # 手绘散点：该 seaborn 版本的 stripplot 不接受 linewidths/width 参数
    rng = np.random.default_rng(0)
    for i, g in enumerate(groups):
        if len(g) == 0:
            continue
        ax.scatter(i + rng.uniform(-0.1, 0.1, len(g)), g, s=1.5, c="black",
                   alpha=0.4, linewidths=0, zorder=3, rasterized=True)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.set_xlim(-0.6, len(order) - 0.4)


def plot_5A(astro, micro):
    fig, axes = plt.subplots(2, 1, figsize=(8, 5))
    for ax, sub, col, ylab, yl in [
        (axes[0], astro, "activate_astro", "Astro activation score",
         Y_LIM["astro"]),
        (axes[1], micro, "activate_micro", "Micro activation score",
         Y_LIM["micro"]),
    ]:
        g = [sub.loc[sub.Age == a, col].values.astype(float)
             for a in AGE_LABELS]
        violins(ax, g, [C_JUV, C_YA, C_OLD],
                ["Juvenile", "Young\nAdult", "Old"])
        style(ax, ylab, yl)
    fig.tight_layout(h_pad=1.5)
    p = os.path.join(OUT, "5A_activation_violin.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("      5A ->", p, flush=True)


def plot_5B(astro, micro):
    fig, axes = plt.subplots(2, 1, figsize=(8, 4.5))
    specs = [
        (axes[0], astro, "activate_astro", "Astro activation score",
         ["Astro-1", "Astro-2"], [C_ASTRO, C_ASTRO], Y_LIM["astro"]),
        (axes[1], micro, "activate_micro", "Micro activation score",
         ["Micro-1", "Micro-2", "Micro-3"], C_MICRO, Y_LIM["micro"]),
    ]
    for ax, sub, col, ylab, order, colors, yl in specs:
        g = [sub.loc[sub.clust_annot == c, col].values.astype(float)
             for c in order]
        violins(ax, g, colors, order)
        style(ax, ylab, yl)
    fig.tight_layout(h_pad=1.5)
    p = os.path.join(OUT, "5B_cluster_violin.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("      5B ->", p, flush=True)


def plot_5D(astro, micro):
    """箱线图 +背后密集黑点带。
    实测原图(500dpi): 点带宽 19~23px、彩色箱体宽仅 11~13px → 点带 ≈1.7× 箱宽，
    所以 jitter 必须 >0；jitter=0 会退化成 1px 细线，丢掉"密密麻麻"效果。
    zorder: 黑点 5 > 中位线 4 > 箱体 3 → 点压在箱体之上（与原图一致）。
    """
    sns.set_style("white")
    fig, axes = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
    for ax, sub, col, ylab, yl in [
        (axes[0], astro, "activate_astro", "Astro activation score",
         Y_LIM["astro"]),
        (axes[1], micro, "activate_micro", "Micro activation score",
         Y_LIM["micro"]),
    ]:
        sub = sub[sub.region.isin(REGION_ORDER)].copy()
        sub["region"] = pd.Categorical(sub.region, categories=REGION_ORDER,
                                      ordered=True)
        # dodge=True 与 boxplot 共用同一套 dodge 几何 → 点带正对各自彩色箱体中心
        sns.stripplot(data=sub, x="region", y=col, hue="Age",
                      order=REGION_ORDER, hue_order=AGE_LABELS,
                      dodge=True, jitter=0.25, palette="dark:black",
                      size=1.0, alpha=0.6, edgecolor="black", linewidth=0,
                      zorder=5, ax=ax)
        sns.boxplot(data=sub, x="region", y=col, hue="Age",
                    order=REGION_ORDER, hue_order=AGE_LABELS,
                    palette=[C_JUV, C_YA, C_OLD], saturation=1.0,
                    ax=ax, width=0.8, linewidth=0.9, whis=1.0,
                    showfliers=False,
                    boxprops=dict(linewidth=0.9, edgecolor="#3a3a3a",
                                  zorder=3),
                    whiskerprops=dict(linewidth=0.9, color="#3a3a3a",
                                      zorder=3),
                    capprops=dict(linewidth=0.9, color="#3a3a3a", zorder=3),
                    medianprops=dict(linewidth=1.2, color="#111", zorder=4))
        ax.set_ylim(*yl)
        ax.set_ylabel(ylab, fontsize=10)
        ax.set_xlabel("")
        ax.grid(False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="both", length=0, labelsize=8.5)
        ax.set_xticks(range(len(REGION_ORDER)))
        ax.set_xticklabels(REGION_ORDER, fontsize=8.5)

    handles = [Line2D([], [], marker="o", linestyle="none", markersize=6,
                      markerfacecolor=c, markeredgecolor="#4a4a4a",
                      markeredgewidth=0.6, label=lab)
               for c, lab in zip([C_JUV, C_YA, C_OLD], AGE_LABELS)]
    axes[0].legend(handles=handles, fontsize=8, frameon=True,
                   facecolor="white", edgecolor="#bbbbbb", framealpha=0.9,
                   loc="upper left", ncol=1, handlelength=0.9,
                   borderpad=0.25, handletextpad=0.4, labelspacing=0.35)
    if axes[1].get_legend():
        axes[1].get_legend().remove()
    fig.tight_layout(h_pad=1.6)
    p = os.path.join(OUT, "5D_region_boxplot.png")
    fig.savefig(p, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("      5D ->", p, flush=True)


# =============================================================================
# 第三部分：统计输出
# =============================================================================
def report(astro, micro):
    lines = []

    def add(t, tab):
        lines.extend(["", "=" * 78, t, "=" * 78, tab.round(4).to_string()])

    lines.append("=" * 78)
    lines.append("校验：幼年组(4wk)均值应≈ 0")
    lines.append("=" * 78)
    lines.append("  activate_astro @ 4wk mean = "
                 f"{astro[astro.age == '4wk'].activate_astro.mean():.8f}")
    lines.append("  activate_micro @ 4wk mean = "
                 f"{micro[micro.age == '4wk'].activate_micro.mean():.8f}")

    add("activate_astro by age", astro.groupby("Age", observed=True)
        ["activate_astro"].agg(["mean", "median", "count"]).reset_index())
    add("activate_micro by age", micro.groupby("Age", observed=True)
        ["activate_micro"].agg(["mean", "median", "count"]).reset_index())
    add("Astro subcluster x age (median)", astro.groupby(
        ["clust_annot", "Age"], observed=True)["activate_astro"].median().unstack())
    add("Micro subcluster x age (median)", micro.groupby(
        ["clust_annot", "Age"], observed=True)["activate_micro"].median().unstack())
    for nm, sub, col in [("ASTRO", astro, "activate_astro"),
                         ("MICRO", micro, "activate_micro")]:
        s = sub[sub.region.isin(REGION_ORDER)]
        add(f"{nm} {col} by region x age", s.groupby(
            ["region", "Age"], observed=True)[col].agg(
            ["mean", "median", "count"]).reset_index())

    txt = "\n".join(lines)
    print("\n" + txt)
    p = os.path.join(OUT, "summary_tables.txt")
    with open(p, "w", encoding="utf-8") as f:
        f.write(txt)


# =============================================================================
def main():
    adata_annot = compute_scores(IN)

    obs = adata_annot.obs.copy()
    obs["region"] = obs["tissue"].astype(str).map(TISSUE_MAP)
    clust = obs["clust_annot"].astype(str)
    astro = obs[clust.str.startswith("Astro")].copy()
    micro = obs[clust.str.startswith("Micro")].copy()
    amap = dict(zip(["4wk", "24wk", "90wk"], AGE_LABELS))
    for s in (astro, micro):
        s["Age"] = pd.Categorical(s["age"].map(amap),
                                  categories=AGE_LABELS, ordered=True)
    print(f"      Astro n={len(astro)}  Micro n={len(micro)}", flush=True)

    plot_5A(astro, micro)
    plot_5B(astro, micro)
    plot_5D(astro, micro)
    report(astro, micro)

    adata_annot.write(os.path.join(OUT, "astro_micro_scores.h5ad"))
    print("\n全部完成 ->", OUT, flush=True)


if __name__ == "__main__":
    main()