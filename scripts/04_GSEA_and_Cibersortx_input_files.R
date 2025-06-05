###############################################################################
# Script: scripts/04_GSEA_and_Cibersortx_input_files.R
# Author: Marina Guillot Fernández
# Last updated: 04-Jun-2025
# Description: Perform differential gene expression analysis using linear
#              modeling (Condition + Sex), generate volcano plot, and save
#              annotated results. Also prepares annotation file and heatmap.
###############################################################################

#-----------------------------#
# Load required libraries     #
#-----------------------------#
library(Biobase)

#-----------------------------#
# Load ExpressionSet object   #
#-----------------------------#
data <- readRDS("results/normalized_filtered_data.rds")  # The variable 'data' should be loaded here

# ---------------------- Generate GSEA files ----------------------
# Extract expression matrix
expr_matrix <- exprs(data)

# Get probe and sample names
probe_names <- featureNames(data)
sample_names <- colnames(expr_matrix)

# Create GCT-like dataframe: Name | Description | Expression values
gct_data <- data.frame(
  Name = probe_names,
  Description = rep("NA", length(probe_names)),
  expr_matrix,
  check.names = FALSE
)

# Write GCT function
write_gct <- function(df, file) {
  num_rows <- nrow(df)
  num_cols <- ncol(df) - 2  # Exclude 'Name' and 'Description'
  
  cat("#1.2\n", file = file)
  cat(num_rows, num_cols, "\n", sep = "\t", file = file, append = TRUE)
  cat("Name", "Description", paste(colnames(df)[-(1:2)], collapse = "\t"), "\n", sep = "\t", file = file, append = TRUE)
  write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

# Export GCT file
write_gct(gct_data, file = "results/expression_set_PD_ClariomS_human.gct")

#-----------------------------------------#
# Generate CHIP annotation file for GSEA
#-----------------------------------------#

# Load Clariom S Human chip annotation file
annotation_file <- read.csv("data/annotation/Clariom_S_Human.r1.na36.hg38.a1.transcript.csv", comment.char = "#")

# Extract relevant annotation columns
chip_file <- annotation_file[, c("probeset_id", "gene_assignment")]

# Parse gene symbols and titles
chip_file$gene_symbol <- sapply(strsplit(chip_file$gene_assignment, " // "), function(x) x[2])
chip_file$gene_title <- sapply(strsplit(chip_file$gene_assignment, " // "), function(x) x[3])

# Create .chip file format: Probe Set ID | Gene Symbol | Gene Title
chip_output <- data.frame(
  `Probe Set ID` = chip_file$probeset_id,
  `Gene Symbol` = chip_file$gene_symbol,
  `Gene Title` = chip_file$gene_title,
  stringsAsFactors = FALSE
)

# Export .chip file
write.table(chip_output, file = "results/PD_ClariomS_human.chip", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

################################################################################
##                     Generate files for CIBERSORTx                         ##
################################################################################

# Select relevant annotation columns
anno <- annotation_file[, c("probeset_id", "gene_assignment")]

# Parse gene symbols
anno$gene_name <- sapply(strsplit(as.character(anno$gene_assignment), " // "), function(x) x[2])

# Rename columns
colnames(anno) <- c("PROBEID", "GENE_ASSIGNMENT", "GENE_NAME")
anno <- anno[, c("PROBEID", "GENE_NAME", "GENE_ASSIGNMENT")]

# Short version: only probe ID and gene name
anno_short <- anno[, c("PROBEID", "GENE_NAME")]

# Merge annotation with expression matrix
data_exprs <- exprs(data)
data_exprs_annotated <- merge(data_exprs, anno_short, by.x = 0, by.y = "PROBEID")  # Row.names become a column

# Reorder columns to place gene name first
data_exprs_cibersort <- data_exprs_annotated[, c("GENE_NAME", setdiff(names(data_exprs_annotated), c("GENE_NAME", "Row.names")))]

# Export annotated matrix to CSV and TSV
write.table(data_exprs_cibersort,
            file = "results/dataexprs_filtered_annotated_CIBERSORTx.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)



