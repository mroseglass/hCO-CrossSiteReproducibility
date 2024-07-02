# cluster integrated dataset for PGP1
# cluster at different resolutions
# save marker genes for each cluster

library(Seurat)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(here)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# seurat objects at various cluster resolutions
output.rds <- paste0(dir.cluster.output, format(Sys.time(), "%Y%m%d"), "_PGP1_UNCR3_reference_integrated_clustered_seurat_object.rds")

# INPUT FILES ##########################################################################################################
# integrated dataset
seurat.integrated.rds <- here("results/seurat/20230314_PGP1_UNCR3_reference_integrated_seurat_object.rds")

# GLOBALS ##############################################################################################################
# cluster resolution (could have been looped, but ran individually to brute force it quickly)
CLUSTER.RESOLUTION <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0)

# get min PCs function
get_min_pc <- function(stdev) {
    sum.stdev <- sum(stdev)
    percent.stdev <- (stdev / sum.stdev) * 100
    cumulative <- cumsum(percent.stdev)

    co1 <- which(cumulative > 90 & percent.stdev < 5)[1]
    co2 <- sort(which((percent.stdev[1:length(percent.stdev) - 1] - percent.stdev[2:length(percent.stdev)]) > 0.1), decreasing = TRUE)[1] + 1

    min.pc <- min(co1, co2)

    return(min.pc)
}


# Load PGP1 Integrated Data ############
printMessage("Loading Seurat Object...")
seur.pgp1 <- readRDS(seurat.integrated.rds)
printMessage("Finished Loading.")


# Cluster Integrated Data #########
printMessage("Running PCA.")

DefaultAssay(seur.pgp1) <- "integrated"
seur.pgp1 <- RunPCA(seur.pgp1,
                    npcs = 20,
                    verbose = FALSE)


min.pc <- get_min_pc(seur.pgp1[["pca"]]@stdev)

print(paste("Number of PCs for UMAP and Neighbors:", min.pc))

printMessage("Running UMAP.")
seur.pgp1 <- RunUMAP(seur.pgp1,
                     dims = 1:min.pc,
                     #dims = 1:10,
                     verbose = TRUE)

printMessage("Running TSNE.")
seur.pgp1 <-RunTSNE(seur.pgp1,
                    dims = 1:min.pc)

printMessage("Running FindNeighbors.")
seur.pgp1 <- FindNeighbors(seur.pgp1,
                           dims = 1:min.pc,
                           #dims = 1:10,
                           verbose = TRUE)

printMessage(paste("Running FindClusters at", paste0(CLUSTER.RESOLUTION, collapse = ","), "resolution."))
seur.pgp1 <- FindClusters(seur.pgp1,
                          verbose = TRUE,
                          algorithm = 2,
                          resolution = CLUSTER.RESOLUTION)

# Save Clustered Data ##################
printMessage("Saving Clustered Seurat Object...")
saveRDS(seur.pgp1, output.rds)
printMessage("Finished Saving Clustered Seurat Object.")


