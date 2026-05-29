# Enrichment Dot Plot

对 **GO** 和 **KEGG** 富集分析结果进行可视化，绘制气泡图（dot plot）。

---

## 功能特点

- **气泡大小**：代表富集到的基因数量（Count）
- **气泡颜色**：代表校正后的 p 值（p.adjust），红 → 蓝表示显著性从高到低
- **GO 分面**：按 BP / CC / MF 三大分类自动分面展示
- **高亮通路**：支持将感兴趣的通路标红加粗
- **自动换行**：长通路名自动换行，避免标签重叠

---

## 依赖包

```r
library(dplyr)
library(stringr)
library(ggplot2)
library(ggtext)    # 用于支持富文本标签（markdown/html）
library(rlang)
library(ggpubr)    # 提供 theme_pubclean()
```

---

## 函数说明

### 1. `plot_go_dotplot()` — GO 富集分析气泡图

```r
plot_go_dotplot(
  data,
  top_n = 10,
  x_var = "FoldEnrichment",
  wrap_length = 50,
  highlight_pathways = NULL,
  facet = TRUE
)
```

| 参数 | 说明 | 默认值 |
|:---|:---|:---|
| `data` | GO 富集结果数据框，需包含：`Description`, `p.adjust`, `Count`, `ONTOLOGY`, 以及 `x_var` 指定的列 | — |
| `top_n` | 每个 ONTOLOGY 分类下展示的 top N 条通路 | `10` |
| `x_var` | X 轴使用的变量名 | `"FoldEnrichment"` |
| `wrap_length` | 通路描述文字的最大宽度，超出自动换行 | `50` |
| `highlight_pathways` | 需要高亮的通路名字符向量，标红加粗 | `NULL` |
| `facet` | 是否按 ONTOLOGY（BP/CC/MF）分面展示 | `TRUE` |

**返回值**：`ggplot` 对象，可继续用 `+` 叠加主题或保存。

---

### 2. `plot_kegg_dotplot()` — KEGG 富集分析气泡图

```r
plot_kegg_dotplot(
  data,
  top_n = 10,
  x_var = "FoldEnrichment",
  wrap_length = 50,
  highlight_pathways = NULL
)
```

| 参数 | 说明 | 默认值 |
|:---|:---|:---|
| `data` | KEGG 富集结果数据框，需包含：`Description`, `p.adjust`, `Count`, 以及 `x_var` 指定的列 | — |
| `top_n` | 全局展示 top N 条通路（KEGG 无分类，不按组） | `10` |
| `x_var` | X 轴使用的变量名 | `"FoldEnrichment"` |
| `wrap_length` | 通路描述文字的最大宽度，超出自动换行 | `50` |
| `highlight_pathways` | 需要高亮的通路名字符向量，标红加粗 | `NULL` |

**返回值**：`ggplot` 对象。

---

## 使用示例

### GO 富集分析可视化

```r
# 基础用法
plot_go_dotplot(go_enrichments)

# 每个分类展示 top 15，高亮指定通路
plot_go_dotplot(
  go_enrichments,
  top_n = 15,
  highlight_pathways = c("cell cycle", "DNA replication", "apoptotic process")
)

# 使用 GeneRatio 作为 X 轴，取消分面
plot_go_dotplot(
  go_enrichments,
  x_var = "GeneRatio",
  facet = FALSE
)
```

### KEGG 富集分析可视化

```r
# 基础用法
plot_kegg_dotplot(kegg_enrichments)

# 展示 top 20，高亮 MAPK 通路
plot_kegg_dotplot(
  kegg_enrichments,
  top_n = 20,
  highlight_pathways = c("MAPK signaling pathway")
)

# 使用 GeneRatio 作为 X 轴
plot_kegg_dotplot(
  kegg_enrichments,
  x_var = "GeneRatio"
)
```

### 保存图片

```r
p <- plot_go_dotplot(go_enrichments, top_n = 15)
ggsave("go_enrichment.png", p, width = 10, height = 8, dpi = 300)
```

---

## 输入数据格式

脚本兼容常见的富集分析工具输出格式，如 `clusterProfiler` 的 `enrichGO()` / `enrichKEGG()` 结果。

### GO 输入数据必备列

| 列名 | 说明 |
|:---|:---|
| `Description` | GO term 描述 |
| `ONTOLOGY` | 分类：BP / CC / MF |
| `p.adjust` | 校正后的 p 值 |
| `Count` | 富集到的基因数 |
| `FoldEnrichment` 或 `GeneRatio` | X 轴变量（取决于 `x_var` 参数） |

### KEGG 输入数据必备列

| 列名 | 说明 |
|:---|:---|
| `Description` | KEGG 通路描述 |
| `p.adjust` | 校正后的 p 值 |
| `Count` | 富集到的基因数 |
| `FoldEnrichment` 或 `GeneRatio` | X 轴变量（取决于 `x_var` 参数） |

> **提示**：若使用 `clusterProfiler`，其输出结果通常已包含上述所有列，可直接传入函数。

---

## 作者 & 维护

如有问题或功能建议，欢迎反馈。
