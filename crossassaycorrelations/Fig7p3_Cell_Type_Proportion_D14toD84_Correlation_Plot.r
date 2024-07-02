#plot cell type proportion to cell type proportion
#Just for graphing
library(tidyverse)
library(readxl)
library(plotrix)
library(devtools)  

#
setwd("/work/users/r/o/roseg")
HERE <- "/work/users/r/o/roseg/single-cell_reproducibility/"
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
results.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/CellTypeProportions/"

celltocell <- read_csv(paste0(results.here,"Cell_Type_Proportion_CorrelationtoEachOther_withGUFFI_Jun22.csv"))
#.x is day 84 and .y is day14
celltocell$cluster <- gsub(".x","", celltocell$cluster)
celltocell$Day14Cluster <- gsub(".y","", celltocell$Day14Cluster)
celltocell$Day14Cluster <- gsub("ler","layer", celltocell$Day14Cluster)
celltocell$Day14Cluster <- gsub("Ler","layer", celltocell$Day14Cluster)

celltocell <- mutate(celltocell, Nominal = case_when(FDR < 0.05 ~ "#",
                                                   pval < 0.05 ~ "*",
                                                   pval > 0.05 ~ ""))
celltocell$cluster <- factor(celltocell$cluster, 
                            levels=c("Neuroepithelial stem cells","Dividing neural progenitor cells, S","Dividing neural progenitor cells, G2",
                                     "Medial Pallium/Marginal Zone","Dividing intermediate progenitors, S","Intermediate progenitors",
                                     "Radial glia","Stressed radial glia","Outer radial glia","Lower layer neuron I",
                                     "Lower layer neuron II","Upper layer neuron I, Layer 2/3/4",
                                     "Upper layer neuron II, layer 2/3","Pan cortical neuron","Pan neuron/Cajal - Retzius",
                                     "Unspecified neuron","Stressed unspecified neuron"))

celltocell$Day14Cluster <- factor(celltocell$Day14Cluster, 
                                  levels=c("Neuroepithelial stem cells","Dividing neural progenitor cells, S","Dividing neural progenitor cells, G2",
                                           "Medial Pallium/Marginal Zone","Dividing intermediate progenitors, S","Intermediate progenitors",
                                           "Radial glia","Stressed radial glia","Outer radial glia","Lower layer neuron I",
                                           "Lower layer neuron II","Upper layer neuron I, layer 2/3/4",
                                           "Upper layer neuron II, layer 2/3","Pan cortical neuron","Pan neuron/Cajal - Retzius",
                                           "Unspecified neuron","Stressed unspecified neuron"))
#for D84 as dotplot with size and color 
celltocell <- mutate(celltocell, negativelog10 = -log10(pval))
#for D84 as dotplot with size and color
cell.hm <- ggplot(celltocell, aes(x = Day14Cluster, y=cluster)) + 
  scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+
  scale_size(limits = c(0,5.6)) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_y_discrete(limits=rev)
cell.hm

cell.hm.r <- cell.hm + 
  geom_text(aes(x = Day14Cluster, 
                y=cluster,
                label = Nominal), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) +
  labs(title = "Day 14 Cell Type Proportions Correlates to Day 84 Cell Type Proportions",
       x = "Day 14 Cell Type", y = "Day 84 Cell Type")
cell.hm.r
ggsave("Jun24_CellTypeProportion_D14_to_D84_reorders.pdf",plot = last_plot(), width = 150, height = 100,
       units = "mm",dpi = 300, device = "pdf",
       path = "/work/users/r/o/roseg/IDDRC/IDDRCPlots/CrossMeasurements")
