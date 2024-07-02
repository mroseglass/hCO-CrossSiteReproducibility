#Just for graphing
library(tidyverse)
library(readxl)
library(plotrix)
library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
setwd("/work/users/r/o/roseg")

# INPUTS ##########################################################################
HERE <- "/work/users/r/o/roseg/single-cell_reproducibility/"
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"

#just one object for all assays
D84Prop.r<- read_csv(paste0(db.here,"Revisions.D84.Proportions.GRUFFI.to.Assays.Technical.Dec28.r.Pearsons.csv"))
D14Prop.r<- read_csv(paste0(db.here,"Revisions.D14.Proportions.GRUFFI.to.Assays.Technical.Dec28.r.Pearsons.csv"))

#remove data if df <4
D14Prop.r <- filter(D14Prop.r, df >3)
D84Prop.r <- filter(D84Prop.r, df >3)

#add labels to object###########################################################
D84Prop.r <- mutate(D84Prop.r, Nominal = case_when(FDR < 0.05 ~ "#",
                                      pval < 0.05 ~ "*",
                                      pval > 0.05 ~ ""))
D84Prop.r$cluster <- factor(D84Prop.r$cluster, 
                            levels=c("Neuroepithelial stem cells","Dividing neural progenitor cells, S","Dividing neural progenitor cells, G2",
                                     "Medial Pallium/Marginal Zone","Dividing intermediate progenitors, S","Intermediate progenitors",
                                     "Radial glia","Stressed radial glia","Outer radial glia","Lower layer neuron I",
                                     "Lower layer neuron II","Upper layer neuron I, Layer 2/3/4",
                                     "Upper layer neuron II, layer 2/3","Pan cortical neuron","Pan neuron/Cajal - Retzius",
                                     "Unspecified neuron","Stressed unspecified neuron"))

D84Prop.r$Assay <- gsub("_deltadeltaCT"," qPCR",D84Prop.r$Assay)

#D84Prop.r$Assay <- gsub("SSEA3+/SSEA4+ -1","iPSC SSEA3+/SSEA4+ Flow",D84Prop.r$Assay)
#D84Prop.r$Assay <- gsub("TRA Backbone -1","iPSC TRA Backbone Flow",D84Prop.r$Assay) 
#D84Prop.r$Assay <- gsub("PercentBud","Visible Budding Day 14",D84Prop.r$Assay)
#D84Prop.r$Assay <- gsub("PercentNest","Growth in Matrigel Day 14",D84Prop.r$Assay)
#D84Prop.r$Assay <- gsub("PercentCystsD35","Expected Morphology Day 35",D84Prop.r$Assay)
#D84Prop.r$Assay <- gsub("PercentCystsD56","Expected Morphology Day 56",D84Prop.r$Assay)

D84Prop.r$Assay <- factor(D84Prop.r$Assay, 
                          levels=c("PassageatThaw","PassageatSeed","SeedViability","SeedingTime","TimeinAccutase",        
                                   "TotalCellCount","NumberWellsforSeeding","TRA backbone -1","SSEA3+SSEA4+ -1",      
                                   "Nanog qPCR","OCT4 qPCR","DUSP6 qPCR","ZIC2 qPCR","OTX2 qPCR","TFCP2L1 qPCR",  
                                   "KLF4 qPCR","TFAP2C qPCR", "Day 07","NumberEBsinCookie","EmbeddingTimes","PercentNest", 
                                   "PercentBud","Volume","PerPAX6","NCADperPAX6","NCADperToPRO","Day 14","Growth714",
                                   "CounthCOatStart", "Day 35","Growth1435",
                                   "PercentCystsD35","CountD35","SOX2+PAX6+ 35","FOXG1 35","TBR2 35","TBR1 35","GSX2 35", 
                                   "Day 56","Growth3556",
                                   "CountD56","PercentCystsD56","scRNAseqAreaD84","DissoTotalCell","DissoViability","CellperuLprefreeze"))

D14Prop.r <- mutate(D14Prop.r, Nominal = case_when(FDR < 0.05 ~ "#",
                                               pval < 0.05 ~ "*"))

D14Prop.r$cluster <- factor(D14Prop.r$cluster, 
                            levels=c("Neuroepithelial stem cells","Dividing neural progenitor cells, S","Dividing neural progenitor cells, G2",
                                     "Medial Pallium/Marginal Zone","Dividing intermediate progenitors, S","Intermediate progenitors",
                                     "Radial glia","Stressed radial glia","Outer radial glia","Lower layer neuron I",
                                     "Lower layer neuron II","Upper layer neuron I, Layer 2/3/4",
                                     "Upper layer neuron II, layer 2/3","Pan cortical neuron","Pan neuron/Cajal - Retzius",
                                     "Unspecified neuron","Stressed unspecified neuron"))

D14Prop.r$Assay <- gsub("_deltadeltaCT"," qPCR",D14Prop.r$Assay)

D14Prop.r$Assay <- factor(D14Prop.r$Assay, 
                          levels=c("PassageatThaw","PassageatSeed","SeedViability","SeedingTime","TimeinAccutase",        
                                   "TotalCellCount","NumberWellsforSeeding","TRA backbone -1","SSEA3+SSEA4+ -1",      
                                   "Nanog qPCR","OCT4 qPCR","DUSP6 qPCR","ZIC2 qPCR","OTX2 qPCR","TFCP2L1 qPCR",  
                                   "KLF4 qPCR","TFAP2C qPCR", "Day 07","NumberEBsinCookie","EmbeddingTimes","PercentNest", 
                                   "PercentBud","Volume","PerPAX6","NCADperPAX6","NCADperToPRO","Day 14","Growth714", "scRNAseqAreaD14",
                                   "CounthCOatStart", "Day 35","Growth1435",
                                   "PercentCystsD35","CountD35","SOX2+PAX6+ 35","FOXG1 35","TBR2 35","TBR1 35","GSX2 35", 
                                   "Day 56","Growth3556",
                                   "CountD56","PercentCystsD56","DissoTotalCell","DissoViability","CellperuLprefreeze"))


#for D84 as dotplot with size and color 
D84Prop.r <- mutate(D84Prop.r, negativelog10 = -log10(pval))

Assay <- as_tibble_col(c("PassageatThaw","PassageatSeed","SeedViability","SeedingTime","TimeinAccutase",        
           "TotalCellCount","NumberWellsforSeeding"), column_name = "Assay")
D84Prop.ipsc <- inner_join(Assay,D84Prop.r)

Assay <- as_tibble_col(c("TRA backbone -1","SSEA3+SSEA4+ -1",      
                         "Nanog qPCR","OCT4 qPCR","DUSP6 qPCR","ZIC2 qPCR","OTX2 qPCR","TFCP2L1 qPCR",  
                         "KLF4 qPCR","TFAP2C qPCR"), column_name = "Assay")
D84Prop.ipsc.gene <- inner_join(Assay,D84Prop.r)

Assay <- as_tibble_col(c("Day 07","PercentNest", 
                         "PercentBud","Day 14","Growth714", "Day 35","Growth1435",
                         "PercentCystsD35","Day 56","Growth3556","PercentCystsD56"), column_name = "Assay")
D84Prop.phasecontrast <- inner_join(Assay,D84Prop.r)

Assay <- as_tibble_col(c("NumberEBsinCookie","EmbeddingTimes","Volume","PerPAX6","NCADperPAX6","NCADperToPRO"), column_name = "Assay")
D84Prop.embedding <- inner_join(Assay,D84Prop.r)

Assay <- as_tibble_col(c("scRNAseqAreaD84","DissoTotalCell","DissoViability","CellperuLprefreeze"), column_name = "Assay")
D84Prop.scRNAseq <- inner_join(Assay,D84Prop.r)

Assay <- as_tibble_col(c("CounthCOatStart","CountD35","SOX2+PAX6+ 35","FOXG1 35","TBR2 35","TBR1 35","GSX2 35"), column_name = "Assay")
D84Prop.1month <- inner_join(Assay,D84Prop.r)

#for D84 as dotplot with size and color
D84P.hm <- ggplot(D84Prop.ipsc, aes(x = Assay, y=cluster)) + 
  scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "iPSC Culture Factor Correlates", x = "Assay", y = "Cell Type")
D84P.ipsc <- D84P.hm + 
  geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()) 

D84P.hm <- ggplot(D84Prop.ipsc.gene, aes(x = Assay, y=cluster)) + scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "iPSC Gene Expression Correlates", x = "Assay", y = "Cell Type")
D84P.ipsc.gene <- D84P.hm + geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()) 

D84P.hm <- ggplot(D84Prop.phasecontrast, aes(x = Assay, y=cluster)) + scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "Phase Contrast Quantification Correlates", x = "Assay", y = "Cell Type")
D84P.phasecontrast <- D84P.hm + geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()) 

D84P.hm <- ggplot(D84Prop.embedding, aes(x = Assay, y=cluster)) + scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "Day 14 Marker Gene and Technical Factor Correlate", x = "Assay", y = "Cell Type")
D84P.embedding <- D84P.hm + geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()) 

D84P.hm <- ggplot(D84Prop.1month, aes(x = Assay, y=cluster)) + scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "Day 35 Marker Gene and Technical Factor Correlates", x = "Assay", y = "Cell Type")
D84P.1month <- D84P.hm + geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()) 

D84P.hm <- ggplot(D84Prop.scRNAseq, aes(x = Assay, y=cluster)) + scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "scRNAseq Factor Correlates", x = "Assay", y = "Cell Type")
D84P.scRNAseq <- D84P.hm + geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()) 


figure <- ggarrange(D84P.ipsc,D84P.ipsc.gene,D84P.phasecontrast,D84P.embedding,
                    D84P.1month,D84P.scRNAseq,
                    labels = c("A", "B"),common.legend = TRUE,
                    ncol = 3, nrow = 2)
figure
ggsave("AprilRevisions_AllCluster_AllAssay_D84_reordered_long.pdf",plot = last_plot(), 
       width = 350, height = 350,
       units = "mm",dpi = 300, device = "pdf",
       path = "/work/users/r/o/roseg/IDDRC/IDDRCPlots/CrossMeasurements")

  

#for D14
D14Prop.r <- mutate(D14Prop.r, negativelog10 = -log10(pval))
D84Prop.r <- D14Prop.r
#for D14 as dotplot with size and color 
Assay <- as_tibble_col(c("PassageatThaw","PassageatSeed","SeedViability","SeedingTime","TimeinAccutase",        
                         "TotalCellCount","NumberWellsforSeeding"), column_name = "Assay")
D84Prop.ipsc <- inner_join(Assay,D84Prop.r)

Assay <- as_tibble_col(c("TRA backbone -1","SSEA3+SSEA4+ -1",      
                         "Nanog qPCR","OCT4 qPCR","DUSP6 qPCR","ZIC2 qPCR","OTX2 qPCR","TFCP2L1 qPCR",  
                         "KLF4 qPCR","TFAP2C qPCR"), column_name = "Assay")
D84Prop.ipsc.gene <- inner_join(Assay,D84Prop.r)

Assay <- as_tibble_col(c("Day 07","PercentNest", 
                         "PercentBud","Day 14","Growth714", "Day 35","Growth1435",
                         "PercentCystsD35","Day 56","Growth3556","PercentCystsD56"), column_name = "Assay")
D84Prop.phasecontrast <- inner_join(Assay,D84Prop.r)

Assay <- as_tibble_col(c("NumberEBsinCookie","EmbeddingTimes","Volume","PerPAX6","NCADperPAX6","NCADperToPRO"), column_name = "Assay")
D84Prop.embedding <- inner_join(Assay,D84Prop.r)

Assay <- as_tibble_col(c("scRNAseqAreaD84","DissoTotalCell","DissoViability","CellperuLprefreeze"), column_name = "Assay")
D84Prop.scRNAseq <- inner_join(Assay,D84Prop.r)

Assay <- as_tibble_col(c("CounthCOatStart","CountD35","SOX2+PAX6+ 35","FOXG1 35","TBR2 35","TBR1 35","GSX2 35"), column_name = "Assay")
D84Prop.1month <- inner_join(Assay,D84Prop.r)

#for D84 as dotplot with size and color
D84P.hm <- ggplot(D84Prop.ipsc, aes(x = Assay, y=cluster)) + 
  scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "iPSC Culture Factor Correlates", x = "Assay", y = "Cell Type")
D84P.ipsc <- D84P.hm + 
  geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) 

D84P.hm <- ggplot(D84Prop.ipsc.gene, aes(x = Assay, y=cluster)) + scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "iPSC Gene Expression Correlates", x = "Assay", y = "Cell Type")
D84P.ipsc.gene <- D84P.hm + geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()) 

D84P.hm <- ggplot(D84Prop.phasecontrast, aes(x = Assay, y=cluster)) + scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "Phase Contrast Quantification Correlates", x = "Assay", y = "Cell Type")
D84P.phasecontrast <- D84P.hm + geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()) 

D84P.hm <- ggplot(D84Prop.embedding, aes(x = Assay, y=cluster)) + scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "Day 14 Marker Gene and Technical Factor Correlate", x = "Assay", y = "Cell Type")
D84P.embedding <- D84P.hm + geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()) 

D84P.hm <- ggplot(D84Prop.1month, aes(x = Assay, y=cluster)) + scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "Day 35 Marker Gene and Technical Factor Correlates", x = "Assay", y = "Cell Type")
D84P.1month <- D84P.hm + geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()) 

D84P.hm <- ggplot(D84Prop.scRNAseq, aes(x = Assay, y=cluster)) + scale_color_gradient2(low = "#2179B4", mid = "white", high = "#DDCB77", midpoint = 0) +
  geom_point(aes(color=r,size=negativelog10))+theme_minimal()+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text.y=element_blank()) +
  scale_y_discrete(limits=rev)+ scale_size(limits = c(0,5.6)) +
  labs(title = "scRNAseq Factor Correlates", x = "Assay", y = "Cell Type")
D84P.scRNAseq <- D84P.hm + geom_text(aes(x = Assay,y=cluster,label = Nominal), color = "black", size = 4) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        panel.border = element_blank(),panel.background = element_blank(),axis.ticks = element_blank()) 


figure <- ggarrange(D84P.ipsc,D84P.ipsc.gene,D84P.phasecontrast,D84P.embedding,
                    D84P.1month,D84P.scRNAseq,
                    labels = c("A", "B"),common.legend = TRUE,
                    ncol = 6, nrow = 1)
figure
#ggsave("Revisions_AllCluster_AllAssay_D14_reordered_long.pdf",plot = last_plot(), width = 250, height = 100,
#       units = "mm",dpi = 300, device = "pdf",
#       path = "/work/users/r/o/roseg/IDDRC/IDDRCPlots/CrossMeasurements")


