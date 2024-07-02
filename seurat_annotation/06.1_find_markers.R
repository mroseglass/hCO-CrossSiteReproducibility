# find marker genes of the clustered data at different resolutions

library(Seurat)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(writexl)
library(here)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# directory for markers at various cluster resolutions
dir.marker.output <- here("results/seurat/cluster_markers/")
dir.create(dir.marker.output, recursive = TRUE, showWarnings = FALSE)

# rescaled seurat object
seurat.rescale.rds <- here("results/seurat/20230321_PGP1_UNCR3_reference_integrated_SCTRescale_seurat_object.rds")

# INPUT FILES ##########################################################################################################
# integrated and clustered dataset
seurat.integrated.rds <- here("results/seurat/20230320_PGP1_UNCR3_reference_integrated_clustered_seurat_object.rds")

# GLOBALS ##############################################################################################################


# Load PGP1 Integrated Data ############
printMessage("Loading Seurat Object...")
seur.pgp1 <- readRDS(seurat.integrated.rds)
printMessage("Finished Loading.")

# get cluster resolution colunn names
cluster.colnames <- colnames(seur.pgp1@meta.data)[grepl("integrated_snn_res", colnames(seur.pgp1@meta.data))]

# Find All Marker Genes ###############

# re-scale data in the SCT slot
seur.pgp1 <- PrepSCTFindMarkers(seur.pgp1,
                                assay = "SCT",
                                verbose = TRUE)

printMessage("Saving Rescaled Data Seurat Object...")
saveRDS(seur.pgp1, seurat.rescale.rds)
printMessage("Finished Saving Rescaled Data Seurat Object.")


# loop over cluster resolutions
for (i in 1:length(cluster.colnames)) {

    printMessage(paste("Finding marker genes for resolution:", cluster.colnames[i]))

    # set clusters by selecting a resolution and setting idents
    Idents(seur.pgp1) <- seur.pgp1@meta.data[,match(cluster.colnames[i], colnames(seur.pgp1@meta.data))]
    print("Cluster idents:")
    print(summary(Idents(seur.pgp1)))

    markers <- FindAllMarkers(seur.pgp1,
                              assay = "SCT",
                              slot = "data",
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              verbose = TRUE)

    # save marker genes
    printMessage(paste("Saving marker genes for resolution:", cluster.colnames[i]))
    saveRDS(markers, paste0(dir.marker.output, format(Sys.time(), "%Y%m%d"), "_cluster_markers_for_", cluster.colnames[i], ".rds"))

}


printMessage("Finished")







