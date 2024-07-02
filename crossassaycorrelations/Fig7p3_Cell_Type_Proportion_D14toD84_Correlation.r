#Cell type proportions at D14 to cell type proportions at D84 
library(tidyverse)
library(readxl)
library(plotrix)
library(lme4)
library(lmtest)
library(ppcor)
library("ggplot2")
library("GGally") 

setwd("/work/users/r/o/roseg")
HERE <- "/work/users/r/o/roseg/single-cell_reproducibility/"
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
results.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/CellTypeProportions"

# INPUTS
D84prop <- read_csv(paste0(db.here,"CellTypeProportions/Day84CellTypePropwithGRUFFI.csv"))
D14prop <- read_csv(paste0(db.here,"CellTypeProportions/Day14CellTypePropwithGRUFFI.csv"))

D84prop$Experiment.scRNAseq <- gsub("_D84_","",D84prop$Experiment.scRNAseq)
D14prop$Experiment.scRNAseq <- gsub("_D14_","",D14prop$Experiment.scRNAseq)

Props <- inner_join(D84prop,D14prop, by = "Experiment.scRNAseq")
#.x is day 84 and .y is day14

list.d14c <- colnames(Props)
list.d84c <- list.d14c[1:18]
list.d84c <- list.d84c[-16]
list.d14c <- list.d14c[19:35]

# pull one cluster to run p-cor with 
y=1
Test <- Props %>% dplyr::select(all_of(c(list.d14c,list.d84c[y])))
Test <- as.data.frame(Test)

#run loops of filtered test object through the assays
i=1
Prop.Result <- cor.test(unlist(Test[ ,length(Test)]), unlist(Test[ ,i]), method="pearson")
r <- Prop.Result$estimate
pval <- Prop.Result$p.value
df <- Prop.Result$parameter
ready <- cbind(r,pval,df,list.d14c[i])

for(i in 2:length(Test)-1) {
  Prop.Result <- cor.test(unlist(Test[ ,length(Test)]), unlist(Test[ ,i]), method="pearson")
  r <- Prop.Result$estimate
  pval <- Prop.Result$p.value
  df <- Prop.Result$parameter
  ready2 <- cbind(r,pval,df,list.d14c[i])
   ready <- rbind(ready2,ready)
}
#make readable and add FDR
ready <- as_tibble(ready)
ready <- dplyr::rename(ready, Day14Cluster = V4)
ready <- ready %>% mutate_at(c('pval', 'r'), as.numeric)
ready <- mutate(ready, FDR = p.adjust(pval, method="BH"))
cluster <- c(rep(list.d84c[y],length(ready$FDR)))
ready <- cbind(cluster,ready)
Prop.Cor.sum <- ready

#Then loop for all other cluster ############################################################
for(y in 2:length(list.d84c)){
Test <- Props %>% dplyr::select(all_of(c(list.d14c,list.d84c[y])))
Test <- as.data.frame(Test)

#run loops of filtered test object through the assays
i=1
Prop.Result <- cor.test(unlist(Test[ ,length(Test)]), unlist(Test[ ,i]), method="pearson")
r <- Prop.Result$estimate
pval <- Prop.Result$p.value
df <- Prop.Result$parameter
ready <- cbind(r,pval,df,list.d14c[i])

for(i in 2:length(Test)-1) {
  Prop.Result <- cor.test(unlist(Test[ ,length(Test)]), unlist(Test[ ,i]), method="pearson")
  r <- Prop.Result$estimate
  pval <- Prop.Result$p.value
  df <- Prop.Result$parameter
  ready2 <- cbind(r,pval,df,list.d14c[i])
  ready <- rbind(ready2,ready)
}
#make readable and add FDR
ready <- as_tibble(ready)
ready <- dplyr::rename(ready, Day14Cluster = V4)
ready <- ready %>% mutate_at(c('pval', 'r'), as.numeric)
ready <- mutate(ready, FDR = p.adjust(pval, method="BH"))
cluster <- c(rep(list.d84c[y],length(ready$FDR)))
ready <- cbind(cluster,ready)
Prop.Cor.sum2 <- ready
Prop.Cor.sum <- rbind(Prop.Cor.sum2,Prop.Cor.sum)
}

#write_csv(Prop.Cor.sum, paste0(results.here,"Cell_Type_Proportion_CorrelationtoEachOther_withGUFFI_Jun22.csv"))

# Check Correlations ####################################################################################
#Pairs plot to look at spread - is huge
#cut down just to one D14 cell type nd
#.x is day 84 and .y is day14

Props <- mutate(Props, Site = gsub("R.*","",Props$Experiment.scRNAseq))
Site.color <-c("#12783D", "#882155", "#332F85")

Testcors <- dplyr::select(Props,c(`Upper layer neuron II, layer 2/3.x` ,`Unspecified neuron.x`,
                                 `Dividing neural progenitor cells, G2.y`))
NPCG2 <- ggpairs(Testcors)

Testcors <- dplyr::select(Props,c(`Upper layer neuron II, layer 2/3.x` ,`Radial glia.x`,`Unspecified neuron.x`,
                                  `Radial glia.y`,))
RG <- ggpairs(Testcors)

Testcors <- dplyr::select(Props,c(`Lower layer neuron I.x`,`Lower layer neuron II.x`,
                                  `Stressed unspecified neuron.y`,`Outer radial glia.y`,`Intermediate progenitors.y`))
LLNs <- ggpairs(Testcors)

pdf("/work/users/r/o/roseg/IDDRC/IDDRCPlots/CrossMeasurements/ExampleCellTypeProportionsAcrossDays.pdf",
    onefile = TRUE, width = 12, height = 8)
print(NPCG2)
print(RG)
print(LLNs)
dev.off()


