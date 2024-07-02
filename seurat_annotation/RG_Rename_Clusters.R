#rename IDDRC clusters
library(Seurat)
library(glmGamPoi)
library(stringr)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(tidytext)
library(RColorBrewer)

library(plyranges)

# OUTPUT FILES #########################################################################################################
dir.obj <- ("/work/users/r/o/roseg/single-cell_reproducibility/results/seurat")
MK.dir <- "/work/users/r/o/roseg/single-cell_reproducibility/"

# INPUT FILES ##########################################################################
# integrated PGP1 seurat object
seurat.pgp1.rds <- paste0(MK.dir,"results/seurat/20230320_PGP1_UNCR3_reference_integrated_clustered_seurat_object.rds")

S.obj <- readRDS(seurat.pgp1.rds)
###################################

#need to set active identity to correct resolution, 0.6 for IDDRC
S.obj <- SetIdent(S.obj, value = S.obj@meta.data$integrated_snn_res.0.6)

mixing_pallete <-c('#5F95B2',
                   '#BCC6E5',
                   '#B09AB1',
                   '#F1C4DC',
                   '#EA9F8B',
                   '#F8C893',
                   '#89A48C',
                   '#369F48',
                   '#DB7B87',
                   '#E12228',
                   '#B177B3',
                   '#2179B4',
                   '#F47B20',
                   '#F89B40',
                   '#F15A29')


S.obj <- RenameIdents(S.obj, `0` = "Neuroepithelial stem cells",
                    `6` = "Dividing neural progenitor cells, S",
                    `10` = "Dividing neural progenitor cells, G2",
                    `12` = "Medial Pallium/Marginal Zone",
                    `14` = "Dividing intermediate progenitors, S",
                    `4` = "Intermediate progenitors",
                    `11` = "Radial glia",
                    `5` = "Outer radial glia",
                    `2` = "Lower layer neuron a",
                    `3` = "Lower layer neuron b",
                    `13` = "Upper layer neuron a, Layer 2/3/4",
                    `1` = "Upper layer neuron b, layer 2/3",
                    `7` = "Pan cortical neuron",
                    `8` = "Pan neuron/Cajal - Retzius",
                    `9` = "Unspecified neuron")

#Save new R object
#saveRDS(obj, paste0(dir.obj, "May12023_IDDRC_PGP1_UNCR3_reference_integrated_clustered_seurat_object.rds"))

#looks correct
DimPlot(S.obj, reduction = "tsne", raster = FALSE, cols = mixing_pallete, label=FALSE)+ NoLegend()
ggsave("June21_ColoredNamesTSNE.png",plot = last_plot(), width = 250, height = 250,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")

# add passage at thaw to metadata and would need to add sample, easier to do in excel sheets
#db <- readRDS("/work/users/r/o/roseg/IDDRC/IDDRCDatabase/IDDRC_culture.rds")
#db <- rename(db, Sample = Experiment.name)
#S.obj@meta.data <- inner_join(S.obj@meta.data, db, by = "Sample")

S.obj@meta.data <- mutate(S.obj@meta.data, PassageatThaw = paste0(S.obj@meta.data$Site,S.obj@meta.data$Rep))
check.meta.data <- S.obj@meta.data

DimPlot(S.obj, reduction = "tsne", group.by = "PassageatThaw", label = FALSE, raster = FALSE)
ggsave("June21_Colored_PassageatSeet_TSNE.png",plot = last_plot(), width = 175, height = 175,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")

Day.color <- c("#F3756D","white","#19BDC2")
DimPlot(S.obj, reduction = "tsne", group.by = "Day", cols = Day.color, label = FALSE, raster = FALSE)
ggsave("June21_Colored_Day_TSNE.png",plot = last_plot(), width = 175, height = 175,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")

Site.color <-c("#F3756D", "#32B44A", "#6F94CD", "#6F94CD")
DimPlot(S.obj, reduction = "tsne", group.by = "WTK_ID", cols = Site.color, label = FALSE, raster = FALSE)
ggsave("June21_Colored_WTK_TSNE.png",plot = last_plot(), width = 250, height = 150,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

S.obj <- CellCycleScoring(S.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(S.obj, reduction = "tsne", group.by = "Phase", cols = Site.color, label = FALSE, raster = FALSE)
check.meta.data <- S.obj@meta.data
ggsave("June21_Colored_CellCycle_TSNE.png",plot = last_plot(), width = 250, height = 150,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")


#markers
DefaultAssay(S.obj) <- "SCT"
FeaturePlot(S.obj, reduction = "tsne", features = "SATB2", raster = FALSE)
ggsave("June21_Colored_SATB2_TSNE.png",plot = last_plot(), width = 175, height = 175,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")

FeaturePlot(S.obj, reduction = "tsne", features = "TBR1", raster = FALSE)
ggsave("June21_Colored_TBR1_TSNE.png",plot = last_plot(), width = 175, height = 175,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")

FeaturePlot(S.obj, reduction = "tsne", features = "EOMES", raster = FALSE)
ggsave("June21_Colored_EOMES_TSNE.png",plot = last_plot(), width = 175, height = 175,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")

FeaturePlot(S.obj, reduction = "tsne", features = "FOXG1", raster = FALSE)
ggsave("June21_Colored_FOXG1_TSNE.png",plot = last_plot(), width = 175, height = 175,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")

FeaturePlot(S.obj, reduction = "tsne", features = "PAX6", raster = FALSE)
ggsave("June21_Colored_PAX6_TSNE.png",plot = last_plot(), width = 175, height = 175,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")

FeaturePlot(S.obj, reduction = "tsne", features = "NES", raster = FALSE)
ggsave("June21_Colored_NES_TSNE.png",plot = last_plot(), width = 175, height = 175,
       units = "mm",dpi = 300, device = "png",
       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/ColorbyTSNE")


