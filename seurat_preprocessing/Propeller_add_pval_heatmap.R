library(tidyverse)

setwd("/work/users/r/o/roseg")
HERE <- "/work/users/r/o/roseg/single-cell_reproducibility/"
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
results.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/CellTypeProportions"

# INPUTS
CellTypeID <- read_csv(paste0(db.here, "ClustertoCellType.csv"))
CellTypeID$Cluster <- as.character(CellTypeID$Cluster)
CellTypeID <- rename(CellTypeID, clusters = Cluster)

D84prop <- read_csv(paste0(db.here,"CellTypeProportions/5823_IDDRCD84Prop.csv"))
D84prop <- rename(D84prop, D84_FDR = FDR)
D84prop <- dplyr::select(D84prop, c(clusters,D84_FDR))
D14prop <- read_csv(paste0(db.here,"CellTypeProportions/5823_IDDRCD14Prop.csv"))
D14prop <- rename(D14prop, D14_FDR = FDR)
D14prop <- dplyr::select(D14prop, c(clusters,D14_FDR))

FDR <- inner_join(D14prop,D84prop)
FDR$clusters <- as.character(FDR$clusters)
FDR <- inner_join(FDR, CellTypeID, by = "clusters")
FDR <- dplyr::select(FDR, !"clusters")
FDR <- FDR %>%
  pivot_longer(!CellType, names_to = "Day", values_to = "padj")

FDR$CellType <- factor(FDR$CellType, 
                            levels=c("Neuroepithelial stem cells","Dividing neural progenitor cells, S","Dividing neural progenitor cells, G2",
                                     "Medial Pallium/Marginal Zone","Dividing intermediate progenitors, S","Intermediate progenitors",
                                     "Radial glia","Outer radial glia","Lower layer neuron a",
                                     "Lower layer neuron b","Upper layer neuron a, Layer 2/3/4",
                                     "Upper layer neuron b, layer 2/3","Pan cortical neuron","Pan neuron/Cajal - Retzius",
                                     "Unspecified neuron"))




D84P.hm <- ggplot(FDR, aes(x=Day, y=CellType, fill=padj)) + 
  geom_tile() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_gradient2(low = "#2179B4", mid = "white",high = "#DDCB77", midpoint = 0.1) +
  labs(title = "FDR", x = "Day")+
  scale_y_discrete(limits=rev)
D84P.hm
