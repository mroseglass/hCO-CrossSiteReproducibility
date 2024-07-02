library(Seurat)
library(glmGamPoi)
library(stringr)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(tidytext)
library(here)
library(plyranges)
# OUTPUT FILES #########################################################################################################
dir.pdf <- ("/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/")
dir.create(dir.pdf, showWarnings = FALSE, recursive = TRUE)
MK.dir <- "/work/users/r/o/roseg/single-cell_reproducibility/"

# output pdf
output.pdf <- paste0(dir.pdf, "20230406_PGP1_UNCR3_reference_integrated_Key_Marker_Genes_expression_on_TSNE.pdf")

# INPUT FILES ##########################################################################################################
# integrated PGP1 seurat object
seurat.pgp1.rds <- paste0(MK.dir,"results/seurat/20230320_PGP1_UNCR3_reference_integrated_clustered_seurat_object.rds")

#ASSIGN CELL CYCLE ##########################################################################
seur.pgp1 <- readRDS(seurat.pgp1.rds)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

iddrc.cc <- CellCycleScoring(seur.pgp1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


p <- FeaturePlot(iddrc.cc, reduction = "tsne", features = "G2M.Score")

print(p)
