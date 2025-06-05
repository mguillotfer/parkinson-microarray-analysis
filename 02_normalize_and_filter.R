###############################################################################
# Script: scripts/02_normalize_and_filter.R
# Author: Marina Guillot Fernández
# Last updated: 04-Jun-2025
# Description: Normalize raw microarray data using RMA, filter out non-annotated
#              and low intensity probes, and perform expression 
#              intensity filtering.
###############################################################################

# ---------------------- Load required libraries ----------------------
required_pkgs <- c("oligo", "Biobase", "matrixStats")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ---------------------PCA plotting function---------------------------

plot_PCA_ggplot_v2 <- function(data, metadata, color_by, top = 500, show_labels = FALSE,
                               save_data = TRUE, save_plot = TRUE, file_name = NULL,
                               n_breaks = 5) {
  
  library(ggplot2)
  library(scales)
  library(limma)
  
  expr_data <- exprs(data)
  meta <- pData(data)
  color_vector <- meta[[color_by]]
  
  # Force numeric if applicable
  if (color_by %in% c("Age", "RIN")) {
    color_vector <- as.numeric(as.character(color_vector))
  }
  is_numeric_color <- is.numeric(color_vector)
  
  # Select probes
  top_probes <- names(sort(apply(expr_data, 1, var), decreasing = TRUE))[1:top]
  expr_data_top <- expr_data[top_probes, ]
  
  # MDS
  mds <- plotMDS(expr_data_top, gene.selection = "common", plot = FALSE)
  var_explained <- round(mds$var.explained * 100, 1)
  
  # DataFrame
  mds_df <- data.frame(
    Sample = colnames(expr_data_top),
    MDS1 = mds$x,
    MDS2 = mds$y,
    Color = color_vector
  )
  
  # Categorical variables
  color_map <- NULL
  if (color_by == "Stage") {
    mds_df$Color <- as.character(mds_df$Color)
    mds_df$Color[is.na(mds_df$Color)] <- "N.D."
    mds_df$Color <- factor(mds_df$Color, levels = c("I", "II", "III", "N.D."))
  }
  
  if (color_by %in% c("Condition", "Stage", "Hospital", "Sex")) {
    color_map <- switch(color_by,
                        "Condition" = c("Control" = "grey50", "Parkinson" = "#D48CF8"),
                        "Stage" = c("I" = "#D48CF8", "II" = "#34b7b9", "III" = "#85E582", "N.D." = "grey50"),
                        "Hospital" = c("12OUH" = "grey50", "GUHA" = "#D48CF8"),
                        "Sex" = c("Male" = "grey50", "Female" = "#D48CF8"))
    mds_df$Color <- factor(mds_df$Color)
  }
  
  # Define color scale
  if (is_numeric_color) {
    color_range <- range(mds_df$Color, na.rm = TRUE)
    breaks_vals <- seq(from = color_range[1], to = color_range[2], length.out = n_breaks)
    
    color_scale <- scale_color_gradient(
      low = "#D48CF8", 
      high = "#34b7b9",
      breaks = breaks_vals,
      labels = label_number(accuracy = 1),
      guide = guide_colorbar(
        title = NULL,
        barwidth = unit(4, "cm"),
        barheight = unit(0.4, "cm"),
        label.theme = element_text(size = 18)
      )
    )
  } else if (!is.null(color_map)) {
    color_scale <- scale_color_manual(values = color_map)
  } else {
    color_scale <- NULL
  }
  
  # File name
  label_tag <- if (show_labels) "_labelsTRUE" else ""
  base_name <- if (is.null(file_name)) {
    paste0("Plots/PCA/PCA_data_", color_by, "_", top, "probes", label_tag, "_", Sys.Date())
  } else {
    tools::file_path_sans_ext(file_name)
  }
  
  csv_file <- paste0(base_name, ".csv")
  probe_file <- paste0("Plots/PCA/PCA_probes_", top, "_probes_", Sys.Date(), ".txt")
  pdf_file <- paste0(base_name, ".pdf")
  
  if (save_data) {
    write.csv(mds_df, file = csv_file, row.names = FALSE)
    write.table(top_probes, file = probe_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  # Plot
  p <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = Color, label = Sample)) +
    geom_point(size = 4, alpha = 0.9) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 19),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20, hjust = 0.5)
    ) +
    labs(
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)"),
      color = NULL
    ) +
    color_scale
  
  if (show_labels) {
    p <- p + geom_text(vjust = -1, size = 3)
  }
  
  if (save_plot) {
    dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(pdf_file, plot = p, width = 3.62, height = 4, units = "in")
    message(paste("✅ Plot saved as:", pdf_file))
  }
  
  return(p)
}

# ---------------------- Define paths ----------------------
results_dir <- "results"
plots_dir <- "Plots"



# ---------------------- Load processed data ----------------------
message("Loading filtered raw data...")
data_reordered <- readRDS(file.path(results_dir, "raw_data.rds"))

# ---------------------- RMA normalization ----------------------
message("Performing RMA normalization...")
data <- oligo::rma(data_reordered)

# ---------------------- Filter probes with no annotation ----------------------
message("Filtering unannotated probes...")

annotation_path <- "data/annotation/Clariom_S_Human.r1.na36.hg38.a1.transcript.csv"
annotation_file <- read.csv(annotation_path, comment.char = "#")

probe_ids <- featureNames(data)
matching_values <- probe_ids %in% annotation_file$probeset_id
matching_probe_ids <- probe_ids[matching_values]
non_matching_probe_ids <- setdiff(probe_ids, annotation_file$probeset_id)

message(length(non_matching_probe_ids), " probes removed (no gene annotation).")

write.csv(non_matching_probe_ids, file = file.path(matrix_dir, "non_matching_probe_ids.csv"))

# Filter expression set
filtered_data <- data[matching_probe_ids, ]
message("Retained ", nrow(filtered_data), " annotated probes.")

write.csv(exprs(filtered_data), file = file.path(matrix_dir, "filteredProbes_norm_data.csv"))

# ---------------------- Intensity-based filtering ----------------------
message("Filtering probes based on expression intensity...")

# Based on: Palmieri et al., Bioconductor workflow for Affymetrix microarrays
palmieri_medians <- rowMedians(exprs(filtered_data))
samples_cutoff <- 16
intensity_threshold <- 4

pdf(file = file.path(plots_dir, "Histogram.pdf"), width = 7, height = 6)
hist(palmieri_medians, 100, col = "cornsilk", freq = FALSE,
     main = "Histogram of the median intensities",
     border = "antiquewhite4",
     xlab = "Median intensities")
abline(v = intensity_threshold, col = "coral4", lwd = 2)
dev.off()

idx_intensity <- apply(exprs(filtered_data), 1, function(x) {
  sum(x > intensity_threshold) >= samples_cutoff
})

table(idx_intensity)

data_filtered <- filtered_data[idx_intensity, ]
message("Retained ", nrow(data_filtered), " probes after intensity filtering.")

message("Normalization and filtering complete.")

# ---------------------- Quality Control: Normalized & Filtered Data ----------------------
message("Generating QC plots for normalized & filtered data...")

qc_dir <- file.path(plots_dir, "QC")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# --- Boxplot of normalized & filtered intensities ---
pdf(file = file.path(qc_dir, "boxplot_intensities_norm_filtered.pdf"), width = 6, height = 6)
par(mar = c(10, 4, 4, 2) + 0.1)
boxplot(data, main = "Boxplot of Intensities (Norm + Filtered)", las = 2, cex.axis = 0.7)
dev.off()

# --- MA plot of normalized & filtered data ---
pdf(file = file.path(qc_dir, "MAplot_norm_filtered.pdf"), width = 6, height = 6)
oligo::MAplot(data)
dev.off()

# --- Density plot across all samples ---
pdf(file = file.path(qc_dir, "density_plot_norm_filtered.pdf"), width = 6, height = 6)
plot(density(exprs(data)[, 1]), col = "steelblue", lwd = 2,
     main = "Density Plot (Norm + Filtered)", xlab = "Expression Value")
for (i in 2:ncol(data)) {
  lines(density(exprs(data)[, i]), col = "steelblue", lwd = 2)
}
dev.off()

# --- PCA Using Custom Function ----------
dir.create("Plots/PCA", recursive = TRUE, showWarnings = FALSE)

# Call custom PCA function
plot_PCA_ggplot_v2(
  data = data,
  metadata = pData(data),
  color_by = "Condition",       # ← Change this to another variable if needed
  top = 500,
  show_labels = TRUE,
  save_data = TRUE,
  save_plot = TRUE
)

message("QC + PCA plots for normalized & filtered data complete.")

# --- Save normalized & filtered data for downstream analysis (Step 3) ---
saveRDS(data_filtered, file = file.path(results_dir, "normalized_filtered_data.rds"))

message("QC + PCA plots for normalized & filtered data complete.")
