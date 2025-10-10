seu_obj <- readRDS("./CD8_2021Nature+2025CC+2023eLife+2023JITC_genefiltered_fastMNN_Celltype.rds")
unique(seu_obj$patient_id)

counts <- GetAssayData(seu_obj, slot = "counts", assay = "RNA")
meta <- seu_obj@meta.data
meta$group <- paste(meta$patient_id, meta$cell_names, sep = "_")
group_levels <- unique(meta$group)

avg_expr_list <- lapply(group_levels, function(g) {
  cells_in_group <- rownames(meta)[meta$group == g]
  if (length(cells_in_group) > 100) {  # 只保留细胞数足够的样本
    rowMeans(counts[, cells_in_group, drop = FALSE])
  } else {
    NULL
  }
})

names(avg_expr_list) <- group_levels

avg_expr_list <- avg_expr_list[!sapply(avg_expr_list, is.null)]

avg_expr_matrix <- do.call(cbind, avg_expr_list)

gene_mean <- rowMeans(avg_expr_matrix)
gene_sd <- apply(avg_expr_matrix, 1, sd)
gene_cv <- gene_sd / gene_mean

gene_stats <- data.frame(
  gene = rownames(avg_expr_matrix),
  mean_expr = gene_mean,
  cv = gene_cv,
  stringsAsFactors = FALSE
)

library(ggplot2)
library(ggrepel)

ggplot(gene_stats, aes(x = mean_expr, y = cv)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +
  labs(title = "Reference Gene Selection", x = "Mean Expression (log10)", y = "Coefficient of Variation (CV)")

genes_to_label <- c("MALAT1",  "RPL28", "TMSB4X", "B2M", "RPL10", "RPLP1")
gene_stats$label <- ifelse(gene_stats$gene %in% genes_to_label, gene_stats$gene, NA)

point_color <- "#6A8CAF"     
label_color <- "#D1495B"     
label_point_color <- "#D1495B" 

ggplot(gene_stats, aes(x = mean_expr, y = cv)) +
  geom_point(color = point_color, alpha = 0.4, size = 1.5) +
  geom_point(
    data = subset(gene_stats, gene %in% genes_to_label),
    color = label_point_color,
    size = 2.5
  ) +
  geom_text_repel(
    aes(label = label),
    size = 4,
    fontface = "bold",
    color = label_color,
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2, 
    segment.size = 0.3
  ) +
  scale_x_log10() +
  labs(
    title = "Reference Gene Selection",
    x = "Mean Expression (log10)",
    y = "Coefficient of Variation (CV)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    panel.grid.major = element_line(size = 0.3, color = "grey90"),
    panel.grid.minor = element_blank()
  )

ref_genes <- gene_stats %>%
  filter(mean_expr > 25, cv < 0.4) %>%
  arrange(cv)

write.table(data.frame(ID=rownames(ref_genes),ref_genes), file = "ref_genes.txt", sep = "\t", col.names = T, row.names = F, quote = F)

library(tidyverse)
ref_gene_names <- c("MALAT1",  "RPL28", "TMSB4X", "B2M", "RPL10", "RPLP1")
ref_gene_names <- c("TCF7",  "FOS", "TOX", "RBPJ", "ENTPD1", "GZMK")

counts <- GetAssayData(seu_obj, assay = "RNA", slot = "counts")
ref_gene_names <- ref_gene_names[ref_gene_names %in% rownames(counts)]

counts_cpm <- t(t(counts) / colSums(counts)) * 1e6
log2_cpm <- log2(counts_cpm + 1)
log2_cpm_sub <- log2_cpm[ref_gene_names, ]

expr_long <- as.data.frame(log2_cpm_sub) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "expression")

p <- expr_long %>%
  ggplot(aes(x = expression)) +
  geom_histogram(aes(y = ..density..), bins = 60, fill = "grey70", color = "black", alpha = 0.7, size = 0.2) +
  geom_density(color = "red", size = 1) +
  facet_wrap(~ gene, scales = "free", ncol = 4) +
  theme_bw(base_size = 14) +
  labs(
    x = "log2(CPM+1) expression in single cells",
    y = "Density",
    title = "Distribution of Stable Genes in Single Cells"
  )

print(p)
p <- expr_long %>%
  ggplot(aes(x = expression)) +
  geom_histogram(aes(y = ..density..), bins = 60, 
                 fill = "#94bfe6", color = "#5073b8", alpha = 0.7, size = 0.18) +  # CNS淡蓝主色
  geom_density(color = "#253494", size = 0.5) +                                 # 深蓝密度线
  facet_wrap(~ gene, scales = "free", ncol = 4) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(size = 0.6, color = "#253494"),
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold", color = "#253494"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13)
  ) +
  labs(
    x = "log2(CPM+1) expression in single cells",
    y = "Density",
    title = "Distribution of Stable Genes in Single Cells"
  )

print(p)




























