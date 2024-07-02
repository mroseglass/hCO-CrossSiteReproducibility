# integration analysis
# perform SCTransform on each sample individually
# then integrate using UNC R3 (D14 and D84) as reference

library(Seurat)
library(glmGamPoi)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(tidytext)
library(here)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# integrated dataset
seurat.integrated.rds <- here("results/seurat/20230314_PGP1_UNCR3_reference_integrated_seurat_object.rds")

# INPUT FILES ##########################################################################################################
# merged seurat filtered file
seurat.merged.filtered.rds <- here("results/seurat/20230313_PGP1_QC_DF_filtered_seurat_object.rds")

# gencode gtf file
#gencode.gtf <- here("data/refgenome/gencode/gencode.v40.annotation.gtf.gz")
# gencode gene names, types, and ensgid
df.gencode.csv <- here("data/refgenome/gencode/gencode.v40.genes.csv.gz")

# GLOBALS ##############################################################################################################
# reference datasets
ref.datasets <- c("UNC_D14_R3", "UNC_D84_R3")



# Load Seurat Object ########
printMessage("Loading Filtered Seurat Object...")
seur.filtered <- readRDS(seurat.merged.filtered.rds)
dim(seur.filtered)


# split object into list by sample
seur.list <- SplitObject(seur.filtered, split.by = "Sample")

# reference indexes
ref.idx <- which(names(seur.list) %in% ref.datasets)

# sample to x cells per samples
# numCells <- 400
# printMessage(paste("Sampling to", numCells, "cells per sample"))
# seur.list <- lapply(X = seur.list, FUN = function(x) {
#     if (length(colnames(x)) < numCells) {
#         x <- x
#     } else {
#         x <- x[, sample(colnames(x), size = numCells, replace = FALSE)]
#     }
# })

# sample to 10% cells per samples (unless thats under 200)
# percentCells <- 0.1
# printMessage(paste("Sampling to", percentCells*100, "percent cells per sample"))
# seur.list <- lapply(X = seur.list, FUN = function(x) {
#     numCells <- length(colnames(x))
#     targetCells <- numCells * percentCells
#     if (targetCells < 200) {
#         x <- x[, sample(colnames(x), size = 200, replace = FALSE)]
#     } else {
#         x <- x[, sample(colnames(x), size = targetCells, replace = FALSE)]
#     }
# })

# SCTransform each object in list
printMessage("SCTransform")
seur.list <- lapply(X = seur.list, FUN = function(x) {
    x <- SCTransform(x,
                     vars.to.regress = c("percent.mt"),
                     return.only.var.genes = FALSE,
                     vst.flavor = "v2",
                     method = 'glmGamPoi',
                     verbose = FALSE)
})

# select features that are variable across datasets
printMessage("Selecting integration features")
features <- SelectIntegrationFeatures(object.list = seur.list,
                                      nfeatures = 5000)

# prep integration
printMessage("Prepping integration")
seur.list <- PrepSCTIntegration(object.list = seur.list,
                                anchor.features = features,
                                verbose = TRUE)

# # run PCA
# printMessage("Running PCA")
# seur.list <- lapply(X = seur.list, FUN = function(x) {
#     x <- RunPCA(x, features = features, verbose = FALSE)
# })

# find anchors
printMessage("Finding integration anchors")
pgp1.anchors <- FindIntegrationAnchors(object.list = seur.list,
                                       reference = ref.idx,
                                       #reduction = "cca",
                                       #dims = 1:50,
                                       normalization.method = "SCT",
                                       anchor.features = features,
                                       verbose = TRUE)

# integrate data
printMessage("Integrating data")
seur.integrated <- IntegrateData(anchorset = pgp1.anchors,
                                 normalization.method = "SCT",
                                 k.weight = 100,
                                 verbose = TRUE)

# save integrated data
printMessage("Saving integrated data")
saveRDS(seur.integrated, seurat.integrated.rds)
printMessage("Done saving data")



















