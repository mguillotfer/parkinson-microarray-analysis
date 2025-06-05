###############################################################################
# Script: scripts/03_differential_expression_analysis.R
# Author: Marina Guillot Fernández
# Last updated: 04-Jun-2025
# Description: Perform differential gene expression analysis using linear
#              modeling (Condition + Sex), generate volcano plot, and save
#              annotated results. Also prepares annotation file and heatmap.
###############################################################################

# ---------------------- Load required libraries ----------------------
required_pkgs <- c("limma", "ggplot2", "dplyr", "ComplexHeatmap", "circlize", "ggrepel", "dendextend")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ---------------------- Load normalized and filtered data ----------------------
message("Loading normalized & filtered expression data...")
data <- readRDS("results/normalized_filtered_data.rds")
metadata <- pData(data)


# ---------------------- Prepare annotation table ----------------------
message(" Preparing annotation table...")

annotation_path <- "data/annotation/Clariom_S_Human.r1.na36.hg38.a1.transcript.csv"
annotation_file <- read.csv(annotation_path, comment.char = "#")

anno <- annotation_file[, c("probeset_id", "gene_assignment")]
anno$gene_name <- sapply(strsplit(as.character(anno$gene_assignment), " // "), function(x) x[2])
colnames(anno) <- c("PROBEID", "GENE_ASSIGNMENT", "GENE_NAME")
anno <- anno[, c("PROBEID", "GENE_NAME", "GENE_ASSIGNMENT")]

write.csv(anno, file = "results/annotation_table.csv", row.names = FALSE)


# ---------------------- Differential Expression Analysis ----------------------
message("Running linear model for Condition + Sex...")

design <- model.matrix(~ Condition + Sex, data = metadata)
fit <- lmFit(data, design)
fit <- eBayes(fit)

results_ConditionSex <- topTable(fit, coef = "ConditionParkinson", adjust = "fdr", number = Inf)

# Merge with annotation
results_annotated <- merge(results_ConditionSex, anno, by.x = "row.names", by.y = "PROBEID")
DE_results <- results_annotated[order(results_annotated$adj.P.Val), ]

# Save results
dir.create("results/Differential_expression", recursive = TRUE, showWarnings = FALSE)
save(DE_results, file = "results/Differential_expression/DE_results_Condition_PDvsCT_Sex.RData")
write.csv(DE_results, file = "results/Differential_expression/DE_resultsCondition_PDvsCT_Sex.csv", row.names = FALSE)

# Significant genes
sig_0.05 <- DE_results[DE_results$adj.P.Val < 0.05, ]
sig_0.1 <- DE_results[DE_results$adj.P.Val < 0.1, ]

message("Significant genes (FDR < 0.05): ", nrow(sig_0.05))
message("Significant genes (FDR < 0.1): ", nrow(sig_0.1))

# ---------------------- Volcano Plot ----------------------
message("Generating volcano plot...")

# Initial table
sortedResultsTable <- DE_results %>%
  mutate(
    # Define colors directly
    Color = case_when(
      adj.P.Val < 0.05 & logFC > 0.5 ~ "red",
      adj.P.Val < 0.05 & logFC < -0.5 ~ "blue",
      adj.P.Val >= 0.05 & adj.P.Val < 0.1 & logFC > 0.5 ~ "orange",
      adj.P.Val >= 0.05 & adj.P.Val < 0.1 & logFC < -0.5 ~ "green",
      TRUE ~ "gray"
    ),
    # Assign point size
    point_size = ifelse(Color == "gray", 0.8, 0.6)
  )


# Filter significant genes for labeling (p-value < 0.1 and |logFC| > 0.5)
significant_genes <- sortedResultsTable %>%
  filter(Color != "gray") %>%
  mutate(
    label_expr = paste0("bolditalic('", GENE_NAME, "')")
  )

# Color palette
colores_volcano <- c(
  "orange" = "#f26926",
  "green" = "#379062",
  "red" = "#be5567",
  "blue" = "#5567be",
  "gray" = "gray50"
)

pdf(file = "plots/VolcanoPlot.pdf", width = 6, height = 6)

ggplot(sortedResultsTable, aes(x = logFC, y = -log10(adj.P.Val))) +
  
  # Significant points
  geom_point(data = subset(sortedResultsTable, Color != "gray"),
             aes(color = Color, size = point_size), alpha = 0.9, show.legend = FALSE) +
  
  # Non-significant points
  geom_point(data = subset(sortedResultsTable, Color == "gray"),
             aes(color = Color, size = point_size), alpha = 0.5, show.legend = FALSE) +
  
  scale_color_manual(values = colores_volcano) +
  scale_size_continuous(range = c(0.6, 0.8)) +
  
  labs(
    x = expression(Log[2]~"(Fold Change)"),
    y = expression(-Log[10]~"(Adjusted p-value)")
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22, face = "bold"),
    plot.margin = margin(4.5, 4.5, .5, .5, "pt"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  ) +
  
  xlim(min(sortedResultsTable$logFC), max(sortedResultsTable$logFC)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), color = "grey", linetype = "dashed") +
  
  scale_x_continuous(
    breaks = seq(floor(min(sortedResultsTable$logFC)), ceiling(max(sortedResultsTable$logFC)), by = 0.5)
  ) +
  scale_y_continuous(
    breaks = seq(0, 5, by = 1),
    limits = c(0, 5) 
  ) +
  
  # Text color = Point color
  geom_text_repel(
    data = significant_genes,
    aes(label = label_expr, color = Color),  # <- Color viene de sortedColors
    parse = TRUE,
    size = 3.5,
    box.padding = 0.3,
    max.overlaps = 20,
    show.legend = FALSE
  )

dev.off()

# ---------------------- Heatmap (Z-score for top DEGs) ----------------------
message("Generating heatmap...")

# Filter genes with p-value < 0.1 and |logFC| > 0.5
filtered_probes <- sortedResultsTable %>%
  filter(adj.P.Val < 0.1 & abs(logFC) > 0.5)

# Retrieve only the corresponding expressed rows
heatmap_data <- data[filtered_probes$Row.names, ]

# Calculate Z-scores
z_scores <- t(apply(heatmap_data, 1, function(x) {
  (x - mean(x)) / sd(x)
}))

# Create styles: all italicized, and significant ones also bolded
gene_labels <- filtered_probes$GENE_NAME
font_styles <- ifelse(filtered_probes$adj.P.Val < 0.05, "bold.italic", "italic")

# Assign row names to the heatmap
rownames(z_scores) <- filtered_probes$GENE_NAME

# Create a vector of customized gpar (graphical parameters)
row_font_styles <- mapply(function(gene, pval) {
  if (pval < 0.05) {
    gpar(fontface = "bold.italic", fontsize = 10)
  } else {
    gpar(fontface = "italic", fontsize = 10)
  }
}, gene = filtered_probes$GENE_NAME, pval = filtered_probes$adj.P.Val, SIMPLIFY = FALSE)

# Original Clustering 
column_dend <- hclust(dist(t(z_scores)))  # Transpuesto: columnas = muestras
dend <- as.dendrogram(column_dend)

#  Reorder dendrogram labels to group CT (if possible)
labels_ordered <- labels(dend)

# Desired order: CT first if possible
desired_order <- c(grep("CT", labels_ordered, value = TRUE),
                   grep("PD", labels_ordered, value = TRUE))

# Reorder dendrogram
dend <- rotate(dend, order = desired_order)

# Create group factor (e.g., CT and PD)
group_labels <- ifelse(grepl("^CT", colnames(data)), "CT", "PD")
group_labels_factor <- factor(group_labels, levels = c("CT", "PD"))

# Create the heatmap with ComplexHeatmap
pdf("plots/Heatmap_DE.pdf", width = 7, height = 6)
Heatmap(
  z_scores,
  name = "Z-score",
  cluster_columns = dend,  # <- dendrograma personalizado
  cluster_rows = TRUE,
  top_annotation = HeatmapAnnotation(
    Group = group_labels,
    col = list(Group = c("CT" = "gray", "PD" = "#D48CF8")),
    annotation_name_side = "left"
  ),
  col = colorRamp2(c(-1.5, 0, 1.5), c("#4558b5", "white", "#b44558")),
  show_row_names = TRUE,
  row_names_gp = gpar(fontface = font_styles, fontsize = 10),
  column_names_gp = gpar(fontsize = 7),
  row_title = "Genes",
  column_title = "Samples",
  heatmap_legend_param = list(title = "Z-score", legend_height = unit(3, "cm"))
)

dev.off()

