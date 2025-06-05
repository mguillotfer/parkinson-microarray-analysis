###############################################################################
# Script: 01_load_and_prepare_data.R
# Author: Marina Guillot Fernández
# Date: 11-Nov-2024 | Last updated: 09-May-2025
# Description: Load raw CEL files, process metadata, 
#               perform quality control and exclude failed sample (PD16)
###############################################################################

# ---------------------- Load required libraries ----------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_pkgs <- c("oligo", "Biobase", "pd.clariom.s.human", "affy", "limma", "ggplot2")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("/mnt/DATOS/JorgeManzanares_FranNavarrete_Sept2024_MFG/Nuevo_analisis_Sexo_cambiado/Codigo_Github")

# ---------------------- Define paths ----------------------
data_dir <- "data"
results_dir <- "results"
qc_dir <- file.path("Plots", "QC")
sample_density_dir <- file.path(qc_dir, "Density_intensities_sample")

dir.create(results_dir, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sample_density_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------- Load CEL files ----------------------
cel_files <- list.celfiles(data_dir, full.names = TRUE)
raw_data <- read.celfiles(cel_files)

# Clean and reorder sample names
sample_names <- colnames(exprs(raw_data))
sample_names_clean <- gsub("_\\(Clariom_S_Human\\).CEL", "", sample_names)
colnames(raw_data) <- sample_names_clean

manual_order <- c("JM-CT1", "JM-CT2", "JM-CT5", "JM-CT6", "JM-CT13", "JM-CT14", "JM-CT15", 
                  "JM-CT16", "JM-CT17", "JM-CT18", "JM-CT19", "JM-CT20", "JM-CT21", 
                  "JM-CT23", "JM-CT25", "JM-CT26", "JM-PK1", "JM-PK2", "JM-PK3", 
                  "JM-PK4", "JM-PK5", "JM-PK6", "JM-PK7", "JM-PK8", "JM-PK9", 
                  "JM-PK16", "JM-PK17", "JM-PK18", "JM-PK19", "JM-PK20", 
                  "JM-PK21", "JM-PK22", "JM-PK23", "JM-PK25", "JM-PK26", 
                  "JM-PK27", "JM-PK28", "JM-PK29", "JM-PK30")

order_index <- match(manual_order, sample_names_clean)
if (anyNA(order_index)) stop("Some sample names in 'manual_order' not found in CEL files.")

data_reordered <- raw_data[, order_index]
message("CEL files loaded and reordered")

# ---------------------- Add metadata ----------------------
message("Adding metadata...")


z <- pData(data_reordered)

# Actualización de metadatos
otherdat <- data.frame(
  SampleName = c("CT1", "CT2", "CT3", "CT4", "CT5", 
                 "CT6", "CT7", "CT8", "CT9", "CT10", 
                 "CT11", "CT12", "CT13", "CT14", "CT15",
                 "CT16", "PD1", "PD2", "PD3", "PD4", 
                 "PD5", "PD6", "PD7", "PD8", "PD9", 
                 "PD10", "PD11", "PD12", "PD13", "PD14", 
                 "PD15", "PD16", "PD17", "PD18", "PD19", 
                 "PD20", "PD21", "PD22", "PD23"),
  Condition = c(rep("Control", 16), rep("Parkinson", 23)),
  
  Hospital = c("12OUH", "12OUH", "GUHA", "GUHA", "12OUH",
               "12OUH", "12OUH", "12OUH", "12OUH", "12OUH",
               "12OUH", "12OUH", "12OUH", "12OUH", "12OUH",
               "12OUH", "12OUH", "12OUH", "12OUH", "12OUH",
               "12OUH", "12OUH", "12OUH", "12OUH", "12OUH",
               "GUHA", "GUHA", "GUHA", "12OUH", "12OUH",
               "12OUH", "12OUH", "12OUH", "12OUH", "GUHA",
               "GUHA", "GUHA", "12OUH", "12OUH"),
  Sex = c("Female", "Female", "Female", "Male", "Female", 
          "Male", "Male", "Male", "Female", "Female", 
          "Female", "Female", "Male", "Male", "Male", 
          "Male", # Hasta aquí van los controles (16)
          "Male", "Female", "Female", "Male", "Male", 
          "Male", "Male", "Female", "Male", "Male", 
          "Male", "Female", "Male", "Female", "Female", 
          "Female", "Female", "Female", "Female", "Female", 
          "Male", "Male", "Male"), 
  Age = c(66, 60, 61, 58, 66, 69, 68, 49, 47, 49, 62, 50, 55, 45, 52, 66,
          71, 72, 62, 65, 82, 75, 106, 81, 74, 71, 69, 70, 77, 63, 71, 84,
          82, 78, 42, 55, 56, 47, 71),
  Stage = c(rep(NA, 16), 
            "I", "II", "I", "I", "II", "I", "II", "II", "I", "I", "II", "I",
            "II", "I", "III", "II", "II", "I", "I", "I", "II", "I", "II"),
  RIN = c(7.3, 8.6, 7.4, 8.4, 8.7, 8.7, 8.7, 7.8, 8.8, 6.4, 6.5, 6.3, 6.3, 6.4, 6.3, 7.5,
          7.6, 7.7, 7.4, 6.7, 8.2, 7.8, 8.2, 8.4, 8.1, 9.7, 9.6, 9.7, 9.3, 5.1, 7.3, 8.1,
          8.9, 6.5, 9.4, 9.1, 9.1, 8.0, 6.4)
  )


z <- data.frame(z, otherdat)
validObject(z)
pData(data_reordered) <- z
colnames(data_reordered) <- pData(data_reordered)$SampleName
message("Metadata added")


# ---------------------- Quality Control Plots ----------------------
message("Generating quality control plots...")

pdf(file = file.path(qc_dir, "boxplot_intensities_raw.pdf"), width = 10, height = 10)
oligo::boxplot(data_reordered, main = "Boxplot of Intensities (Raw Data)", las = 2, cex.axis = 0.7)
dev.off()

pdf(file = file.path(qc_dir, "MAplot_rawdata.pdf"), width = 10, height = 10)
oligo::MAplot(data_reordered)
dev.off()

pdf(file = file.path(qc_dir, "density_intensities_raw.pdf"), width = 10, height = 10)
affy::hist(data_reordered)
dev.off()

# Expression density histograms per sample
for (x in seq_len(ncol(data_reordered))) {
  sample_name <- colnames(data_reordered)[x]
  png_file <- file.path(sample_density_dir, paste0("density_sample_", x, "_", sample_name, ".png"))
  pdf_file <- file.path(sample_density_dir, paste0("density_sample_", x, "_", sample_name, ".pdf"))
  
  png(png_file, width = 800, height = 600)
  affy::hist(data_reordered[, x],
             main = paste("Histogram -", sample_name),
             xlab = "Intensity", ylab = "Frequency")
  dev.off()
  
  pdf(pdf_file, width = 8, height = 6)
  affy::hist(data_reordered[, x],
             main = paste("Histogram -", sample_name),
             xlab = "Intensity", ylab = "Frequency")
  dev.off()
}

# ---------------------- Filter failed sample ----------------------
message(" Excluding sample PD16 due to QC failure...")

data <- data_reordered[, !(sampleNames(data_reordered) %in% "PD16")]

# ---------------------- Save processed data ----------------------
message("Saving processed data...")
saveRDS(data, file = file.path(results_dir, "raw_data.rds"))
write.csv(pData(data), file = file.path(results_dir, "sample_metadata.csv"))
write.csv(exprs(data), file = file.path(results_dir, "raw_expression_matrix.csv"))

# ---------------------- Save Session Info ----------------------
sink(file.path(results_dir, "sessionInfo.txt"))
sessionInfo()
sink()

message("Done! Data prepared and QC plots saved.")
