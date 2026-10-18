# Figure 5A / 5B / 5D 复现

复现对象：庄小威实验室 Cell 2022
**Molecular and spatial signatures of mouse brain aging at single-cell resolution**
（Allen, Blosser, Sullivan, Dulac, Zhuang）

代码来源：<https://github.com/ZhuangLab/SpatialBrainAgingCell22>
作者 notebook：`notebooks/merfish_spatial_celltype_org.ipynb`（CELL 41–43）

---

## 一、运行

```bash
python reproduce_figure5.py
```

依赖：`scanpy` / `anndata` / `seaborn` / `matplotlib` / `numpy` / `pandas`

输入数据（默认路径，可在脚本顶部 `IN` 修改）：
```
C:\Users\yyj\Downloads\c93d78c2-ee17-4504-8d1c-17cf093ad7b5.h5ad
378,918 细胞 × 374 基因
```

输出到脚本同级 `figure5_out/`：

| 文件 | 内容 |
|---|---|
| `5A_activation_violin.png` | Astro / Micro × 三年龄小提琴图 |
| `5B_cluster_violin.png` | Astro-1/2、Micro-1/2/3 小提琴图 |
| `5D_region_boxplot.png` | 8 脑区 × 三年龄箱线图 |
| `astro_micro_scores.h5ad` | 含`activate_astro`/`activate_micro` 的全量数据 |
| `summary_tables.txt` | 全部统计汇总 |

---

## 二、作者原始代码（复现的唯一依据）

**CELL 41 — 算分**
```python
sc.tl.score_genes(adata_annot,
    gene_list=['B2m','Trem2','Ccl2','Apoe','Axl','Itgax',
               'Cd9','C1qa','C1qc','Lyz2','Ctss'],
    score_name='activate_micro', use_raw=False)
sc.tl.score_genes(adata_annot,
    gene_list=['C4b','C3','Serpina3n','Cxcl10','Gfap','Vim','Il18','Hif3a'],
    score_name='activate_astro', use_raw=False)
```

**CELL 42 / 43 — 以 4wk 为基准归零**
```python
adata_micro = adata_annot[adata_annot.obs.cell_type=="Micro"]
adata_annot.obs.activate_micro -= np.mean(
    adata_micro[adata_micro.obs.age=='4wk'].obs.activate_micro)

adata_astro = adata_annot[adata_annot.obs.cell_type=="Astro"]
adata_annot.obs.activate_astro -= np.mean(
    adata_astro[adata_astro.obs.age=='4wk'].obs.activate_astro)
```

### 三个必须注意的细节

1. **不做 Z-score。** 论文 STAR Methods 文字提到 "Z-scored gene expression values"，
   但作者代码里没有 `sc.pp.scale`，直接 `use_raw=False` 在 log 数据上算分。**以代码为准。**
2. **不设 `random_state` / `ctrl_size`**，用 scanpy 默认（`ctrl_size=50`）。
3. **归零作用于 `obs` 全部细胞**（不置 NaN），基准是各自细胞类型中 4wk 批次的均值。

---

## 三、踩坑记录

### 1. `cell_type` 标签与作者不一致 → 全 NaN

作者代码用 `cell_type=="Micro"` / `"Astro"`，但本数据集 `obs.cell_type` 用的是**全称**：

```
astrocyte         37,480
microglial cell   16,514
```

直接照抄会取到空集。解决：建语义映射（脚本中 `CT_MAP`）。
细胞数校验一致 —— `astrocyte` = Astro-1 + Astro-2，`microglial cell` = Micro-1/2/3。

### 2. matplotlib `violinplot` 不支持 `cut`

`plt.violinplot()` 没有 `cut` 参数（seaborn 才有）。
若要严格裁剪 KDE 尾部，改用 `scipy.stats.gaussian_kde` 手绘；
本复现用 seaborn 的 `cut=0`。

### 3. `stripplot` 参数冲突

该版本 seaborn 中 `stripplot(linewidths=)` 与内部传入的 `linewidth` 别名冲突报
`TypeError: Got both 'linewidth' and 'linewidths'`；
`stripplot(width=...)` 会被透传到 `PathCollection.set()` 报
`AttributeError: unexpected keyword argument 'width'`。
解决：用 `linewidth=0`，并让 stripplot 与 boxplot 共用 dodge 几何自动对齐。

### 4. 5D 散点 jitter 不能为 0

实测论文原图（500 dpi 裁剪 + 逐列统计暗像素）：

| | 宽度 |
|---|---|
| 黑色点带 | **19–23 px** |
| 彩色箱体 | **11–13 px** |

点带比箱体**更宽**（约 1.7×）且极其密集。`jitter=0` 会让所有点重叠成1 px 细线，
完全丢掉"密密麻麻"的效果。最终采用 `jitter=0.25, size=1.0, alpha=0.6`。

### 5. zorder 必须显式设定

`stripplot` 未传 zorder 时其 `PathCollection` 会拿到 zorder=3，与 `boxplot` 相同，
导致层级不可控。目标效果是**黑点压在箱体之上**，故显式设`zorder=5`（箱体 3、中位线 4）。

---

## 四、绘图参数

| 位置 | 参数 |
|---|---|
| 5A/5B 小提琴 | `sns.violinplot(cut=0, bw_adjust=1)` |
| 5A/5B 散点 | 手绘 `scatter(s=1.5, alpha=0.4, jitter=±0.1)` — 中央窄竖带 |
| 5D 散点 | `sns.stripplot(dodge=True, jitter=0.25, s=1.0, alpha=0.6, zorder=5)` |
| 纵轴 | astro `(-1.3, 3.1)`、micro `(-1.3, 2.1)`（下限留白，小提琴底端悬空） |
| 图例 | 左上角，竖排，无标题，彩色圆点 |

### 配色（采样自论文原图 PDF，非目测）

| 用途 | 色值 |
|---|---|
| Juvenile | `#739ADA` 蓝 |
| Young Adult | `#D7C2D6` 淡紫 |
| Old | `#E08E8F` 鲑红 |
| Astro-1/2 | `#367F5A` 深绿（同色） |
| Micro-1/2/3 | `#892B5D` → `#B42E77` → `#E23290` 品红渐变 |

提取脚本思路：400 dpi 裁剪 PDF → 按颜色连通块自动定位各 violin 中心 → 取众数色。

---

## 五、复现结果

校验（幼年组均值应 ≈ 0）：
```
activate_astro @ 4wk mean = 0.00000000
activate_micro @ 4wk mean = 0.00000000
```

| | 4wk | 24wk | 90wk |
|---|---|---|---|
| activate_astro 均值 | 0.0000 | +0.1288 | **+0.7123** |
| activate_astro 中位数 | −0.0917 | +0.0046 | +0.6982 |
| activate_micro 均值 | 0.0000 | −0.1113 | **+0.4697** |
| activate_micro 中位数 | +0.0161 | −0.1028 | +0.4499 |

**亚群中位数**

| | 4wk | 24wk | 90wk |
|---|---|---|---|
| Astro-1 | −0.2413 | −0.1276 | +0.3642 |
| Astro-2 | +0.7701 | +1.0758 | **+1.4494** |
| Micro-1 | +0.0405 | +0.0141 | +0.3122 |
| Micro-2 | −0.2193 | −0.2745 | +0.3010 |
| Micro-3 | +0.2612 | +0.1906 | **+0.6017** |

范围：astro [−1.99, 4.13]、micro [−2.19, 2.94]
细胞数：Astro 37,480（12,831 / 8,352 / 16,297）；Micro 16,514（6,317 / 3,496 / 6,701）

**核心结论**：astro 激活随年龄单调升高；micro 在 24wk 略低于 4wk、到 90wk 才大幅跃升
—— "青年期短暂回落 + 老年激增" 的模式在两类胶质细胞中一致。

---

## 六、制作过程中不成功的案例

以下中间版本均已废弃，保留作过程记录：

| 问题 | 原因 |
|---|---|
| 分数整体偏高、5D 须线过长 | 漏了作者的 4wk 归零步骤 |
| 5D 箱体与黑点错位 | 手写 dodge 偏移公式与 seaborn 0.13 几何不一致 → 改用原生 `dodge=True` |
| 5D 黑点稀疏不像原图 | `jitter=0` 退化成 1px 细线 → 改 `jitter=0.25` |
| 图层反复翻转 | `stripplot` 未显式设 zorder，与箱体同为 3 → 显式 `zorder=5` |
| 小提琴底部碰轴线 | 未预留底部空白 → ylim 下限改 −1.3 |
