# Get average expression for each cluster from Kreigsten

library(Seurat)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(here)
library(ComplexHeatmap)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################

dir.pdf <- here("doc/seurat/pdf/")


# INPUT FILES ##########################################################################################################
# integrated and clustered dataset
seurat.integrated.rds <- here("results/seurat/20230320_PGP1_UNCR3_reference_integrated_clustered_seurat_object.rds")

# integrated and clustered dataset, rescaled seurat object
seurat.rescale.rds <- here("results/seurat/20230321_PGP1_UNCR3_reference_integrated_SCTRescale_seurat_object.rds")


# gencode gene names, types, and ensgid
df.gencode.csv <- here("data/refgenome/gencode/gencode.v40.genes.csv.gz")

# marker genes
df.markers.csv <- here("data/marker_genes/marker_genes.csv")

# GLOBALS ##############################################################################################################



# Load PGP1 Integrated Data ############
printMessage("Loading Seurat Object...")
seur.pgp1 <- readRDS(seurat.integrated.rds)
printMessage("Finished Loading.")

# get cluster resolution column names
cluster.colnames <- colnames(seur.pgp1@meta.data)[grepl("integrated_snn_res", colnames(seur.pgp1@meta.data))]


# Load Marker Genes ###########
# marker genes
df.markers <- read_csv(df.markers.csv)

# gencode genes
df.gencode <- read_csv(df.gencode.csv)

print("Gene names not matching in dataset:")
df.markers$gene_name[!df.markers$gene_name %in% df.gencode$gene_name]

# filter out markers without matching name
df.markers %<>%
    filter(gene_name %in% df.gencode$gene_name)



# loop over cluster resolutions
for (i in 1:length(cluster.colnames)) {


    printMessage(paste("Plotting PGP1 violin plots for resolution:", cluster.colnames[i]))

    # set clusters by selecting a resolution and setting idents
    Idents(seur.pgp1) <- seur.pgp1@meta.data[,match(cluster.colnames[i], colnames(seur.pgp1@meta.data))]
    print("Cluster idents:")
    print(summary(Idents(seur.pgp1)))

    DefaultAssay(seur.pgp1) <- "SCT"

    pdf(paste0(dir.pdf, "20230321_PGP1_UNCR3_reference_", cluster.colnames[i], "_violin_plots.pdf"), height = 6, width = 8)

    for (j in 1:nrow(df.markers)) {

        if (!df.markers$gene_name[j] %in% rownames(seur.pgp1)) {
            print(paste("Gene:", df.markers$gene_name[j], "not in seurat object. Skipping."))
            next()
        }

        p <- VlnPlot(seur.pgp1,
                features = df.markers$gene_name[j],
                raster = NULL,
                assay = "RNA",
                slot = "counts",
                pt.size = 0.05) + labs(y = "Counts")


        print(p)

    }

    dev.off()

}

