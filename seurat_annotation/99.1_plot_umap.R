# plot UMAP

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

#library(mikelaffr)

# OUTPUT FILES #########################################################################################################
dir.pdf <- ("/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/")
dir.create(dir.pdf, showWarnings = FALSE, recursive = TRUE)
MK.dir <- "/work/users/r/o/roseg/single-cell_reproducibility/"

# output pdf
output.pdf <- paste0(dir.pdf, "20230406_PGP1_UNCR3_reference_integrated_Key_Marker_Genes_expression_on_TSNE.pdf")

# INPUT FILES ##########################################################################################################
# integrated PGP1 seurat object
seurat.pgp1.rds <- paste0(MK.dir,"results/seurat/20230320_PGP1_UNCR3_reference_integrated_clustered_seurat_object.rds")

# gencode gene names, types, and ensgid
df.gencode.csv <- paste0(MK.dir,"data/refgenome/gencode/gencode.v40.genes.csv.gz")

# marker genes
df.markers.csv <- paste0(MK.dir,"data/marker_genes/HumanOrganoidThalamus.csv")

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


# Load Marker Genes ###########
# marker genes
df.markers <- read_csv(df.markers.csv)
df.markers <- dplyr::rename(df.markers, gene_name = gene)
#try for cluster makers of "NSC", "hindbrain", "dorsal", "neural stem cells", "ventral", "upper layer"
#df.markers <- df.markers %>% filter(str_detect(cluster, "neural stem cells"))


# gencode genes
df.gencode <- read_csv(df.gencode.csv)

print("Gene names not matching in dataset:")
df.markers$gene_name[!df.markers$gene_name %in% df.gencode$gene_name]

# gencode genes
df.gencode <- read_csv(df.gencode.csv)

# filter out markers without matching name
df.markers <- df.markers %>%
    filter(gene_name %in% df.gencode$gene_name)

# Load Seurat Object #############
printMessage("Loading Seurat Object...")
seur.pgp1 <- readRDS(seurat.pgp1.rds)
printMessage("Finished Loading.")

pdf(output.pdf, height = 5, width = 7)

DefaultAssay(seur.pgp1) <- "SCT"
printMessage("Plotting")

#DimPlot(seur.pgp1,
#        group.by = "Day",
#        label = FALSE,
#        reduction = "umap")

#DimPlot(seur.pgp1,
#        group.by = "Site",
#        label = FALSE,
#        reduction = "umap")

#DimPlot(seur.pgp1,
#        group.by = "WTK_ID",
#        label = FALSE,
#        reduction = "umap")

#DimPlot(seur.pgp1,
#        group.by = "Sample",
#        label = FALSE,
#        reduction = "umap")

DimPlot(seur.pgp1,
        group.by = "Matrigel",
        label = FALSE,
        reduction = "tsne")

DimPlot(seur.pgp1,
        group.by = "Site",
        label = FALSE,
        reduction = "tsne")

DimPlot(seur.pgp1,
        group.by = "WTK_ID",
        label = FALSE,
        reduction = "tsne")

DimPlot(seur.pgp1,
        group.by = "Sample",
        label = FALSE,
        reduction = "tsne")

#for (i in 1:5) {
for (i in 1:nrow(df.markers)) {

    if (!df.markers$gene_name[i] %in% rownames(seur.pgp1)) {
        print(paste("Gene:", df.markers$gene_name[i], "not in seurat object. Skipping."))
        next()
    }

    #p <- FeaturePlot(seur.pgp1, reduction = "umap", features = df.markers$gene_name[i])

    #print(p)

    p <- FeaturePlot(seur.pgp1, reduction = "tsne", features = df.markers$gene_name[i])

    print(p)

}

dev.off()


printMessage("Finished")

