###############################################################################
# Script: scripts/05_Comparison_with_other_studies.R
# Author: Marina Guillot Fernández
# Last updated: 04-Jun-2025
# Description: Volcano and Scatter Plots Comparing our Differential Expression 
#              Results with Public Datasets from Craig et al. (2021)
#              and Calligaris et al. (2015)
#
#
# References: Craig, D. W. et al. (2021). RNA sequencing of whole blood reveals
#             early alterations in immune cells and gene expression in
#             Parkinson’s disease. Nature Aging, 1(8), 734–747.
#             Calligaris, R., Banica, M., Roncaglia, P., Robotti, E., 
#             Finaurini, S., Vlachouli, C., ... & Gustincich, S. (2015). Blood
#             transcriptomics of drug-naive sporadic Parkinson’s disease patients.
#             BMC genomics, 16, 1-14.
# 
###############################################################################
# PACKAGES
library(ggplot2)
library(ggrepel)
library(dplyr)

#------------Comparison with Craig et al. (2021) -----------------------------

craig_df <- read.csv("data/Supplementary_Table3.csv", sep = "\t")

# Ensure numeric values
craig_df$logFC <- as.numeric(gsub(",", ".", craig_df$logFC))
craig_df$adj.P.Val <- as.numeric(gsub(",", ".", craig_df$adj.P.Val))
craig_df$negLog10AdjP <- -log10(craig_df$adj.P.Val)

# Rename columns for consistency
colnames(craig_df) <- c("gene_id_Craig", "Gene_name", "gene_type_Craig", "logFC_Craig",
                        "AveExpr_Craig", "t_Craig", "P.Value_Craig", "adj.P.Val_Craig",
                        "B_Craig", "AvgExpr_PD_Craig", "AvgExpr_CT_Craig", "negLog10AdjP_Craig")

# ------------------------------
# Load and filter our differential expression results
# ------------------------------

our_results <- read.csv("results/Differential_expression/DE_resultsCondition_PDvsCT_Sex.csv")
our_results <- our_results[, !colnames(our_results) %in% "GENE_ASSIGNMENT"]

colnames(our_results) <- c("PROBEID_Navarrete", "logFC_Navarrete", "AveExpr_Navarrete", 
                           "t_Navarrete", "P.Value_Navarrete", "adj.P.Val_Navarrete", 
                           "B_Navarrete", "Gene_name")

highlight_genes <- our_results$Gene_name[our_results$adj.P.Val_Navarrete < 0.1]

# Merge both datasets
merged_data <- merge(craig_df, our_results, by = "Gene_name")

# ------------------------------
# Pearson Correlation
# ------------------------------

cor_test <- cor.test(merged_data$logFC_Craig, merged_data$logFC_Navarrete, method = "pearson")

cat("Pearson Correlation:\n")
cat("r =", round(cor_test$estimate, 3), "\n")
cat("p-value =", format.pval(cor_test$p.value, digits = 3), "\n")
cat("n =", cor_test$parameter + 2, "\n")


# ------------------------------
# Scatter plot highlighting overlap
# ------------------------------

merged_data <- merged_data %>%
  mutate(
    category = case_when(
      adj.P.Val_Navarrete < 0.1 & logFC_Craig > 0 & logFC_Navarrete > 0.5 ~ "up_both",
      adj.P.Val_Navarrete < 0.1 & logFC_Craig < 0 & logFC_Navarrete < -0.5 ~ "down_both",
      adj.P.Val_Navarrete < 0.1 & sign(logFC_Craig) != sign(logFC_Navarrete) ~ "discordant",
      TRUE ~ "normal"
    )
  )

colors <- c(
  "up_both" = "#be5567",
  "down_both" = "#5567be",
  "discordant" = "#379062",
  "normal" = "grey70"
)

pdf("Plots/ScatterPlot_Comparison_Craig.pdf", width = 6, height = 6)

ggplot(merged_data, aes(x = logFC_Craig, y = logFC_Navarrete)) +
  geom_point(data = filter(merged_data, category == "normal"),
             aes(color = category), size = 0.6, alpha = 0.5, show.legend = FALSE) +
  geom_point(data = filter(merged_data, category != "normal"),
             aes(color = category), size = 0.8, alpha = 0.9) +
  geom_text_repel(data = filter(merged_data, category != "normal"),
                  aes(label = Gene_name, color = category),
                  size = 4, fontface = "bold", box.padding = 0.4, 
                  point.padding = 0.3, segment.color = "grey50", show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#11111190") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#11111190") +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("Pearson r = ", round(cor_test$estimate, 3)),
           hjust = -0.1, vjust = 1.1, size = 4.5, fontface = "italic", color = "black") +
  labs(x = expression(Log[2]~"(Fold Change, Craig)"),
       y = expression(Log[2]~"(Fold Change, Navarrete)"),
       color = NULL) +
  scale_color_manual(values = colors, labels = c(
    "up_both" = "Upregulated in both",
    "down_both" = "Downregulated in both",
    "discordant" = "Opposite directions"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

dev.off()


# ---------------- Volcano Plot (Craig) ----------------
merged_data <- merged_data %>%
  mutate(direction = case_when(
    adj.P.Val_Navarrete < 0.1 & logFC_Craig > 0 & logFC_Navarrete > 0.5 ~ "up_both",
    adj.P.Val_Navarrete < 0.1 & logFC_Craig < 0 & logFC_Navarrete < -0.5 ~ "down_both",
    adj.P.Val_Navarrete < 0.1 & sign(logFC_Craig) != sign(logFC_Navarrete) ~ "discordant",
    TRUE ~ NA_character_
  ))

highlighted <- filter(merged_data, !is.na(direction))
others <- filter(merged_data, is.na(direction))
y_lim <- ceiling(max(merged_data$negLog10AdjP_Craig, na.rm = TRUE))

for (size in c(5, 6)) {
  pdf(sprintf("Plots/Metaanalysis_CraigPPMI/VolcanoPlot_Craig_%dx%d.pdf", size, size), width = size, height = size)
  
  ggplot() +
    geom_point(data = others, aes(x = logFC_Craig, y = negLog10AdjP_Craig),
               color = "grey70", size = 0.6, alpha = 0.5) +
    geom_point(data = highlighted, aes(x = logFC_Craig, y = negLog10AdjP_Craig, color = direction),
               size = 0.8, alpha = 0.9) +
    geom_text_repel(data = highlighted, aes(x = logFC_Craig, y = negLog10AdjP_Craig,
                                            label = Gene_name, color = direction),
                    size = 3.5, fontface = "bold.italic", box.padding = 0.5, 
                    point.padding = 0.3, segment.color = "grey50", show.legend = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#11111190") +
    labs(x = expression(Log[2]~"(Fold Change)"),
         y = expression(-Log[10]~"(Adjusted P-value)")) +
    scale_color_manual(values = scatter_colors, name = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14),
      legend.position = "bottom",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      plot.margin = margin(5, 15, 5, 5, "pt")
    ) +
    scale_x_continuous(limits = c(-0.75, 0.6), breaks = seq(-0.75, 0.5, 0.25)) +
    scale_y_continuous(limits = c(0, y_lim), breaks = seq(0, y_lim, 2.5))
  
  dev.off()
}

#------------Comparison with Calligaris et al. (2015) -----------------------------

# Load Calligaris PUMA table
puma <- read.csv("data/12864_2015_2058_MOESM6_ESM_PUMA.csv", sep = "\t")
puma$Gene.symbol <- as.character(puma$Gene.symbol)

# Load our differential expression results again
manzanares <- read.csv("Differential_expression/Condition_sex/DE_resultsCondition_PDvsCT_Sex.csv")
manzanares <- manzanares[, -ncol(manzanares)]
colnames(manzanares) <- paste0(colnames(manzanares), "_manzanares")
colnames(manzanares)[colnames(manzanares) == "GENE_NAME_manzanares"] <- "Gene.symbol"

colnames(puma) <- paste0(colnames(puma), "_calligaris_puma")
colnames(puma)[colnames(puma) == "Gene.symbol_calligaris_puma"] <- "Gene.symbol"
colnames(puma)[colnames(puma) == "Fold.change..LOG2._calligaris_puma"] <- "logFC_calligaris"

# Merge and process
merged_puma <- merge(puma, manzanares, by = "Gene.symbol", all = FALSE)
merged_puma$logFC_calligaris <- as.numeric(gsub(",", ".", merged_puma$logFC_calligaris))
merged_puma$logFC_manzanares <- as.numeric(gsub(",", ".", merged_puma$logFC_manzanares))

# Correlation
cor_test2 <- cor.test(merged_puma$logFC_calligaris, merged_puma$logFC_manzanares, method = "pearson")
print(cor_test2)

# Plot
highlighted_genes <- c("AQP3")
merged_puma$highlight <- ifelse(
  merged_puma$Gene.symbol %in% highlighted_genes, "Highlighted",
  ifelse(
    (merged_puma$logFC_calligaris * merged_puma$logFC_manzanares > 0) &
      (abs(merged_puma$logFC_calligaris) > 0.1 & abs(merged_puma$logFC_manzanares) > 0.1),
    "Concordant", "Normal"
  )
)

ggplot(merged_puma, aes(x = logFC_calligaris, y = logFC_manzanares)) +
  geom_point(data = subset(merged_puma, highlight == "Normal"), color = "grey70", size = 0.6, alpha = 0.4) +
  geom_point(data = subset(merged_puma, highlight == "Concordant"), color = "#be5567", size = 0.8) +
  geom_point(data = subset(merged_puma, highlight == "Highlighted"), color = "#7e303d", size = 1.0) +
  geom_text_repel(data = subset(merged_puma, highlight == "Highlighted"), 
                  aes(label = Gene.symbol), color = "#b9475a", size = 4, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#11111190") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#11111190") +
  annotate("text", x = -Inf, y = Inf,
           label = paste("Pearson r =", round(cor_test2$estimate, 3)),
           hjust = -0.1, vjust = 1.1, size = 4.5, fontface = "italic") +
  labs(x = expression(Log[2]~"(Fold Change, Calligaris)"),
       y = expression(Log[2]~"(Fold Change, Navarrete)")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
    axis.ticks = element_line(color = "black", size = 0.5)
  )

ggsave("Plots/Metaanalisis_Calligaris/ScatterPlot_Calligaris_Navarrete.pdf", width = 5, height = 3.5, dpi = 300)
