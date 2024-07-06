#!/usr/bin/env Rscript
###
 # @author Andi Liu
 # @email andi.liu@uth.tmc.edu
 # @create date 2024-07-05 23:20:22
 # @modify date 2024-07-05 23:20:22
 # @desc [description]
###

args = commandArgs(trailingOnly=TRUE)
setwd(args)
# Seurat processing in R
library(Seurat)
library(dplyr)
library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
library(future)
library(stringr)
library(tidydr)
library(DoubletFinder)
library(sctransform)
ldir <- args
csvlist <- list.files(path = ldir, pattern = 'csv')
for (x in csvlist) {
  counts_data <- read.csv(x, header = TRUE, row.names = 1)
  counts_data <- t(counts_data)
  obj <- CreateSeuratObject(counts = counts_data, min.cells = 3, min.features = 200)
  sob_re <- paste(x, '.rds', sep = '')
  saveRDS(obj, file = sob_re)
}
# Automatically generate the list of .rds files
files_to_load <- list.files(path = ldir, pattern = 'csv.rds')


#Create a list to hold the loaded objects
object_list <- list()

for (file in files_to_load) {
  object_name <- sub("\\..*$", "", file)
  object_list[[object_name]] <- readRDS(file)
}

sample_ids <- sub(".csv.rds", "", basename(files_to_load))
obj <- merge(
  x = object_list[[1]],
  y = object_list[-1],
  add.cell.ids = sample_ids
)

# Further processing
obj <- PercentageFeatureSet(obj, '^MT.', col.name = 'percent.mt')
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
meta <- obj@meta.data

# Load the processed Seurat object
obj <- JoinLayers(obj)
# Normalize and find variable features
#obj <- SCTransform(obj)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = 'vst', nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:10, verbose = FALSE)

# Add individual and cell barcode information
obj$individualID <- str_split_fixed(rownames(obj@meta.data), pattern = '_', n = 2)[,1]
obj$cell_barcode <- str_split_fixed(rownames(obj@meta.data), pattern = '_', n = 2)[,2]
obj$individual_cell_id <- rownames(obj@meta.data)

# DoubletFinder processing
sweep.res.list <- paramSweep(obj, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

# Estimate homotypic doublet proportion
nExp_poi <- round(0.05 * nrow(obj@meta.data))  # Assuming 5% doublet formation rate
pK_optimal <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric),]$pK))
obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = pK_optimal, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

# Add doublet information
meta <- obj@meta.data
meta$doublet <- ifelse(meta[,9] == 'Doublet', TRUE, FALSE)
obj@meta.data <- meta
# Load reference data
load(paste0('../Reference/', '02_reference_SEA_AD_1234_subclass.RData'))
load(paste0('../Reference/','01_Reference_ROSMAP_Control.RData'))

# Update meta data for references
meta_ref_seaad <- reference_1234_subclass@meta.data
meta_ref_methy <- reference@meta.data

# Change broad cell type names to major cell type names
meta_ref_methy$major.cell.type <- NA
meta_ref_methy[meta_ref_methy$broad.cell.type == 'Ex',]$major.cell.type <- 'Excitatory'
meta_ref_methy[meta_ref_methy$broad.cell.type == 'In',]$major.cell.type <- 'Inhibitory'
meta_ref_methy[meta_ref_methy$broad.cell.type == 'Ast',]$major.cell.type <- 'Astrocyte'
meta_ref_methy[meta_ref_methy$broad.cell.type == 'End',]$major.cell.type <- 'Endothelial'
meta_ref_methy[meta_ref_methy$broad.cell.type == 'Mic',]$major.cell.type <- 'Microglia'
meta_ref_methy[meta_ref_methy$broad.cell.type == 'Oli',]$major.cell.type <- 'Oligodendrocyte'
meta_ref_methy[meta_ref_methy$broad.cell.type == 'Opc',]$major.cell.type <- 'OPC'
meta_ref_methy[meta_ref_methy$broad.cell.type == 'Per',]$major.cell.type <- 'VLMC/Per'

# Update SEA-AD meta data
meta_ref_seaad[meta_ref_seaad$major.cell.type == 'Microglia-PVM',]$major.cell.type <- 'Microglia'
meta_ref_seaad[meta_ref_seaad$major.cell.type == 'VLMC',]$major.cell.type <- 'VLMC/Per'

# Assign updated meta data back to references
reference_1234_subclass@meta.data <- meta_ref_seaad
reference@meta.data <- meta_ref_methy

# Transfer cell type labels from SEA-AD reference to query

obj <- SCTransform(obj)

DefaultAssay(obj) <- 'SCT'
transfer_anchors_seaad <- FindTransferAnchors(
  reference = reference_1234_subclass,
  query = obj,
  normalization.method = 'SCT',
  query.assay = 'SCT',
  recompute.residuals = FALSE,
  dims = 1:30
)
Idents(reference_1234_subclass) <- 'major.cell.type'
predictions_seaad <- TransferData(
  anchorset = transfer_anchors_seaad, 
  reference = reference_1234_subclass,
  refdata = reference_1234_subclass$major.cell.type,
  query = obj,
  dims = 1:30
)
obj$major.cell.type.prediction.seaad <- predictions_seaad$predicted.id
obj$major.cell.type.prediction.seaad.score <- predictions_seaad$predicted.id.score
# Transfer cell type labels from Mathy et al. reference to query
transfer_anchors_mathy <- FindTransferAnchors(
  reference = reference,
  query = obj,
  normalization.method = 'SCT',
  query.assay = 'SCT',
  recompute.residuals = FALSE,
  dims = 1:30
)
Idents(reference) <- 'major.cell.type'
predictions_mathy <- TransferData(
  anchorset = transfer_anchors_mathy, 
  reference = reference,
  refdata = reference$major.cell.type,
  query = obj,
  dims = 1:30
)
obj$major.cell.type.prediction.mathy <- predictions_mathy$predicted.id
obj$major.cell.type.prediction.mathy.score <- predictions_mathy$predicted.id.score

# Visualize annotations
p1 <- DimPlot(obj, label = TRUE, repel = TRUE, group.by = 'major.cell.type.prediction.seaad')
p2 <- DimPlot(obj, label = TRUE, repel = TRUE, group.by = 'major.cell.type.prediction.mathy')

ggsave("cell_annotation_plot.pdf", plot = patchwork::wrap_plots(p1, p2, ncol=2), scale = 1, width = 13, height = 6)

# Filter based on prediction scores
obj <- subset(obj, major.cell.type.prediction.mathy.score > 0.95)
obj <- subset(obj, major.cell.type.prediction.seaad.score > 0.95)

# Save the final processed object

diagnosis_mapping <- read.table("../donor_disease.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("individualID", "Diagnosis"))
obj@meta.data <- obj@meta.data %>%
  left_join(diagnosis_mapping, by = "individualID") %>%
  mutate(Diagnosis = ifelse(is.na(Diagnosis), "Unknown", Diagnosis))  # Assign "Unknown" if no matching individualID
saveRDS(obj, file='final_processed_seurat_obj.rds')
