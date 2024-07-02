# Doublet Filter pgp1 samples

library(Seurat)
library(DoubletFinder)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(here)

library(mikelaffr)

HERE <- "/work/users/r/o/roseg/single-cell_reproducibility/"
# OUTPUT FILES #########################################################################################################
# qc and doublet filtered seurat file
seurat.qc.df.rds <- paste0(HERE,"results/seurat/20230313_PGP1_QC_DF_filtered_seurat_object.rds")

# qc and doublet LABELED (doublets are still in the object) seurat file
seurat.qc.dfLabeled.rds <- paste0(HERE,"results/seurat/20230313_PGP1_QC_DF_labeled_seurat_object.rds")

# INPUT FILES ##########################################################################################################
# merged seurat filtered file, filtered for mitochondrial read percentage, gene count, and read count
seurat.merged.filtered.rds <- paste0(HERE,"results/seurat/20230313_PGP1_QC_filtered_seurat_object.rds")


# GLOBALS ##############################################################################################################
# get min PCs
get_min_pc <- function(stdev) {
    sum.stdev <- sum(stdev)
    percent.stdev <- (stdev / sum.stdev) * 100
    cumulative <- cumsum(percent.stdev)

    co1 <- which(cumulative > 90 & percent.stdev < 5)[1]
    co2 <- sort(which((percent.stdev[1:length(percent.stdev) - 1] - percent.stdev[2:length(percent.stdev)]) > 0.1), decreasing = TRUE)[1] + 1

    min.pc <- min(co1, co2)

    return(min.pc)
}


# Load Seurat Object ########
printMessage("Loading Seurat Object...")
seur.object <- readRDS(seurat.merged.filtered.rds)
printMessage("Finished Loading Seurat Object.")

print(paste("Number of cells:", dim(seur.object)[2]))

# Doublet Filter ##############
printMessage("Running DoubletFilter")
# split object into list by sample
seur.list <- SplitObject(seur.object, split.by = "Sample")

pdf(here("doc/seurat/pdf/20230313_doublet_filter.pdf"))
# loop over each sample and run doublet filter
for (i in 1:length(seur.list)) {
    printMessage(paste(i, "of", length(seur.list), ":", names(seur.list)[i]))

    wtk.id <- seur.list[[i]]@meta.data$WTK_ID[1]

    if (wtk.id %in% c("WTK1", "WTK2", "WTK3")) {
        exp.D.ratio <- 0.034
    } else {
        exp.D.ratio <- 0.0425
    }

    seur.list[[i]] <- SCTransform(seur.list[[i]],
                                  vars.to.regress = c("percent.mt"),
                                  verbose = FALSE,
                                  return.only.var.genes = FALSE,
                                  vst.flavor = "v2",
                                  method = "glmGamPoi")

    seur.list[[i]] <- RunPCA(seur.list[[i]],
                             npcs = 20,
                             verbose = FALSE)

    min.pc <- get_min_pc(seur.list[[i]][["pca"]]@stdev)

    seur.list[[i]] <- RunUMAP(seur.list[[i]],
                              dims = 1:min.pc,
                              verbose = FALSE)

    seur.list[[i]] <- FindNeighbors(seur.list[[i]],
                                    dims = 1:min.pc,
                                    verbose = FALSE)

    seur.list[[i]] <- FindClusters(seur.list[[i]],
                                   resolution = 0.1,
                                   verbose = FALSE)

    # pK identification (no ground-truth)
    sweep.res.list <- paramSweep_v3(seur.list[[i]], PCs = 1:min.pc, sct = TRUE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)

    # Optimal pK is the max of bomdality coefficent (BCmvn) distribution
    pK <- as.numeric(as.character(bcmvn$pK))
    BCmetric <- bcmvn$BCmetric
    optimal.pk <- pK[which(BCmetric %in% max(BCmetric))]

    # Homotypic doublet proportion estimate
    nExp.poi <- round(exp.D.ratio * nrow(seur.list[[i]]@meta.data))

    # DoubletFinder
    seur.list[[i]] <- doubletFinder_v3(seur.list[[i]],
                                       pN = 0.25,
                                       pK = optimal.pk,
                                       nExp = nExp.poi,
                                       PCs = 1:min.pc,
                                       sct = TRUE)

    colnames(seur.list[[i]]@meta.data)[grepl("DF.classifications", colnames(seur.list[[i]]@meta.data))] <- "DF.name"

    p <- VlnPlot(seur.list[[i]], features = "nFeature_RNA", group.by = "DF.name", pt.size = 0.4) + ggtitle(names(seur.list)[i])

    print(p)

}

dev.off()

printMessage("Saving seurat list...")
saveRDS(seur.list, here("results/seurat/20230313_TEMP_DF_seur_list.rds"))
printMessage("Finished Saving.")

# seur.list <- readRDS(here("results/seurat/TEMP_seur_list_DF.rds"))

# get doublets
printMessage("Getting Doublets")
df.doublets <- data.frame()

for (i in 1:length(seur.list)) {

    df.doublets <- bind_rows(df.doublets, seur.list[[i]]@meta.data[,c("cell_name", "DF.name")])

}

# remove doublet list
rm(seur.list)

# add doublet data to seurat object

printMessage("Adding doublet data to metadata")
head(df.doublets)
head(df.doublets[rownames(seur.object@meta.data),])

df.doublets <- df.doublets[rownames(seur.object@meta.data),]

print("all rownames of metadata == rownames of doublet table")
all(rownames(seur.object@meta.data) == rownames(df.doublets))

seur.object@meta.data$DF.name <- df.doublets$DF.name

# save doublet labeled qc data

printMessage("Saving doublet labeled data")
saveRDS(seur.object, seurat.qc.dfLabeled.rds)

# filter out doublets
printMessage("Filtering doublets")
dim(seur.object)
seur.object <- subset(seur.object, subset = DF.name == "Singlet")
dim(seur.object)

# save doublet filtered qc data
printMessage("Saving doublet filtered QC filtered seurat object")
saveRDS(seur.object, seurat.qc.df.rds)

printMessage("Finished")
