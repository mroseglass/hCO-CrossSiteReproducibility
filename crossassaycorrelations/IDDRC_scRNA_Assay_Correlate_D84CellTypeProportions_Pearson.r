#Other assays & technical variables to cell type proportions at D14 including GRUFFI
library(tidyverse)
library(readxl)
library(plotrix)
library(lme4)
library(lmtest)
library(ppcor)
setwd("/work/users/r/o/roseg")
HERE <- "/work/users/r/o/roseg/single-cell_reproducibility/"
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"

# perform lm correlations for all assays ###########################################################
#load other assays 
D84.Prop.FC.Rank <- read_csv(paste(db.here,"D84PropFCRank.csv", sep = "")) %>% dplyr::select(!Experiment.name)
#add average area and pull out specific scRNAseq
area.scrnaseq <- read_csv("/work/users/r/o/roseg/IDDRC/IDDRCDatabase/CrossSectionalAreaRNAseq.csv") %>%
  filter(FileName != "Copy of 20211008-EW38R2-Day14-C4-scRNAIntensity.jpeg") %>%
  filter(FileName != "Copy of 20211008-EW38R2-Day14-C5-scRNAIntensity.jpeg") %>%
  mutate(Experiment = paste0(Site,Replicate)) %>%
  filter(Site != "CN_First") %>%
  dplyr::select(c(Experiment,Day,AreaMM2))

area.scrnaseq <- pivot_wider(area.scrnaseq, names_from = Day, values_from = AreaMM2)%>%
  rename(scRNAseqAreaD14 = D14) %>%
  rename(scRNAseqAreaD84 = D84)

AreaGrowthMeasures <- read_csv("/work/users/r/o/roseg/IDDRC/IDDRCDatabase/AverageAreaMeasures.csv") %>% 
  dplyr::select(!c(Site,Replicate))
D84.Prop.FC.Rank <- left_join(D84.Prop.FC.Rank,AreaGrowthMeasures)
D84.Prop.FC.Rank <- left_join(D84.Prop.FC.Rank,area.scrnaseq)

#unique(D14.Prop.FC.Rank$Experiment)

#use lm to see if any biological variables correlated to D84 proportion
Prop.Area <- matrix()
Prop.Area.Result <- list()

# FIRST EXAMPLE FOR CLUSTER TO ASSAY ####################################################
#partial correlation because only testing 2 variables
#create new data fram with one cluster plus all assays 
#add 2 cluster lists
D84.Prop.FC.Rank <- dplyr::select(D84.Prop.FC.Rank, !c("hCOAge","scRNAseqAreaD14"))
cluster.list <- D84.Prop.FC.Rank[,1:15]
cluster.list <- cbind(cluster.list,D84.Prop.FC.Rank[,32:33])
cluster.list <- colnames(cluster.list)
assay.list <- colnames(D84.Prop.FC.Rank[,17:31])
assay.list2 <- colnames(D84.Prop.FC.Rank[,34:64])
assay.list <- c(assay.list,assay.list2)


# pull one cluster to run p-cor with 
y=1
Test <- D84.Prop.FC.Rank %>% dplyr::select(all_of(c(assay.list,cluster.list[y])))
Test <- as.data.frame(Test)
#Test <- column_to_rownames(Test,"Experiment")

#run loops of filtered test object through the assays
i=1
Prop.Area.Result <- cor.test(unlist(Test[ ,length(Test)]), unlist(Test[ ,i]), method="pearson")
r <- Prop.Area.Result$estimate
pval <- Prop.Area.Result$p.value
df <- Prop.Area.Result$parameter
ready <- cbind(r,pval,df,assay.list[i])

  for(i in 2:length(Test)-1) {
    Prop.Area.Result <- cor.test(unlist(Test[ ,length(Test)]), unlist(Test[ ,i]), method="pearson")
    r <- Prop.Area.Result$estimate
    pval <- Prop.Area.Result$p.value
    df <- Prop.Area.Result$parameter
    ready2 <- cbind(r,pval,df,assay.list[i])
    ready <- rbind(ready2,ready)
  }
#make readable and add FDR
ready <- as_tibble(ready)
ready <- dplyr::rename(ready, Assay = V4)
ready <- ready %>% mutate_at(c('pval', 'r'), as.numeric)
ready <- mutate(ready, FDR = p.adjust(pval, method="BH"))
cluster <- c(rep(cluster.list[y],length(ready$FDR)))
count <- 
ready <- cbind(cluster,ready)
assay.D14.prop.sum <- ready

#Then loop for all other cluster ############################################################
for(y in 2:length(cluster.list)){
Test <- D84.Prop.FC.Rank %>% dplyr::select(all_of(c(assay.list,cluster.list[y])))
Test <- as.data.frame(Test)
#Test <- column_to_rownames(Test,"Experiment")

#run loops of filtered test object through the assays
i=1
Prop.Area.Result <- cor.test(unlist(Test[ ,length(Test)]), unlist(Test[ ,i]), method="pearson")
r <- Prop.Area.Result$estimate
pval <- Prop.Area.Result$p.value
df <- Prop.Area.Result$parameter
ready <- cbind(r,pval,df,assay.list[i])

for(i in 2:length(Test)-1) {
  Prop.Area.Result <- cor.test(unlist(Test[ ,length(Test)]), unlist(Test[ ,i]), method="pearson")
  r <- Prop.Area.Result$estimate
  pval <- Prop.Area.Result$p.value
  df <- Prop.Area.Result$parameter
  ready2 <- cbind(r,pval,df,assay.list[i])
  ready <- rbind(ready2,ready)
}
#make readable and add FDR
ready <- as_tibble(ready)
ready <- dplyr::rename(ready, Assay = V4)
ready <- ready %>% mutate_at(c('pval', 'r'), as.numeric)
ready <- mutate(ready, FDR = p.adjust(pval, method="BH"))
cluster <- c(rep(cluster.list[y],length(ready$FDR)))
ready <- cbind(cluster,ready)
assay.D14.prop.sum2 <- ready
assay.D14.prop.sum <- rbind(assay.D14.prop.sum2,assay.D14.prop.sum)
}

#write.csv(assay.D14.prop.sum, paste0(db.here,"Revisions.D84.Proportions.GRUFFI.to.Assays.Technical.Dec28.r.Pearsons.csv"))

# CORRELATION FOLLOW UP ########################################################################################
D14.check <- mutate(D14.Prop.FC.Rank, Replicate = paste0("Rep",str_extract(D14.Prop.FC.Rank$Experiment,"[0-9]")))
D14.check <- mutate(D14.check, Site = gsub("*Rep[0-9]","",D14.check$Experiment))
D14.check$Site <- gsub("P1","",D14.check$Site)
# GLOBALS ##########################################################################################
Site.color <-c("#12783D", "#882155", "#332F85")
Replicate.shapes <- c(1,2,3,4,15,16,17)

#Check if SSEAand ZIC2 correlate
Prop.Area.Result <- cor.test(D14.check$`SSEA3+SSEA4+ -1`, D14.check$`Dividing neural progenitor cells, S`, method="pearson")
mainpoint <- paste("r =",round(Prop.Area.Result$estimate, digits = 2),
                   "p-value = ",round(Prop.Area.Result$p.value, digits = 6))

D84toPRbyExpt.p <- (ggplot(D14.check, aes(x=`SSEA3+SSEA4+ -1`, y=`Dividing neural progenitor cells, S`))) +
  scale_color_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  geom_smooth(method = lm, color = "black") +
  geom_point(size = 2, aes(col = Site, shape = Replicate)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw()+
  labs(title = "SSEA Backbone and NPC in S Phase in iPSCs",
       subtitle = mainpoint)+
  xlab("SSEA Backbone") +
  ylab("Proportion NPC in S Phase")
D84toPRbyExpt.p
ggsave("Revisions_AllCluster_AllAssay_D84_reordered_SSEA_NPCinS.pdf",plot = last_plot(), width = 100, height = 100,
       units = "mm",dpi = 300, device = "pdf",
       path = "/work/users/r/o/roseg/IDDRC/IDDRCPlots/CrossMeasurements")

#Check if LLNb and ZIC2 correlate
Prop.Area.Result <- cor.test(D14.check$Nanog_deltadeltaCT, D14.check$`Dividing intermediate progenitors, S`, method="pearson")
mainpoint <- paste("r =",round(Prop.Area.Result$estimate, digits = 2),
                   "p-value = ",round(Prop.Area.Result$p.value, digits = 6))

D84toPRbyExpt.p <- (ggplot(D14.check, aes(x=Nanog_deltadeltaCT, y=`Dividing intermediate progenitors, S`))) +
  scale_color_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  geom_smooth(method = lm, color = "black") +
  geom_point(size = 2, aes(col = Site, shape = Replicate)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw()+
  labs(title = "iPSC Nanog and IP in S Phase at Day 84",
       subtitle = mainpoint)+
  xlab("Nanog Delta CT") +
  ylab("IP in S Phase")
D84toPRbyExpt.p
ggsave("Revisions_AllCluster_AllAssay_D84_reordered_Nanog_IPinS.pdf",plot = last_plot(), width = 100, height = 100,
       units = "mm",dpi = 300, device = "pdf",
       path = "/work/users/r/o/roseg/IDDRC/IDDRCPlots/CrossMeasurements")

#Check if DUSP6 and NPCG2 correlate
Prop.Area.Result <- cor.test(D14.check$PercentBud, D14.check$`Pan cortical neuron`, method="pearson")
mainpoint <- paste("r =",round(Prop.Area.Result$estimate, digits = 2),
                   "p-value = ",round(Prop.Area.Result$p.value, digits = 6))

D84toPRbyExpt.p <- (ggplot(D14.check, aes(x=PercentBud, y=`Pan cortical neuron`))) +
  scale_color_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  geom_smooth(method = lm, color = "black") +
  geom_point(size = 2, aes(col = Site, shape = Replicate)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw()+
  labs(title = "Day 14 Visual Quality and Pan Neurons at Day 84",
       subtitle = mainpoint)+
  xlab("Percent Budding Day 14") +
  ylab("Pan Cortical Neuron")
D84toPRbyExpt.p
ggsave("Revisions_AllCluster_AllAssay_D84_reordered_Budding_PanNeuron.pdf",plot = last_plot(), width = 100, height = 100,
       units = "mm",dpi = 300, device = "pdf",
       path = "/work/users/r/o/roseg/IDDRC/IDDRCPlots/CrossMeasurements")

#Check if DUSP6 and NPCG2 correlate
Prop.Area.Result <- cor.test(D14.check$`Day 35`, D14.check$`Dividing neural progenitor cells, S`, method="pearson")
mainpoint <- paste("r =",round(Prop.Area.Result$estimate, digits = 2),
                   "p-value = ",round(Prop.Area.Result$p.value, digits = 6))

D84toPRbyExpt.p <- (ggplot(D14.check, aes(x=`Day 35`, y=`Dividing neural progenitor cells, S`))) +
  scale_color_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  geom_smooth(method = lm, color = "black") +
  geom_point(size = 2, aes(col = Site, shape = Replicate)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw()+
  labs(title = "Day 35 Average Cross-sectional Area and NPC in S Phase at Day 84",
       subtitle = mainpoint)+
  xlab("Day 35 Average Cross-sectional Area") +
  ylab("Proportion NPC in S Phase")
D84toPRbyExpt.p
ggsave("Revisions_AllCluster_AllAssay_D84_reordered_D35Area_NPCinS.pdf",plot = last_plot(), width = 100, height = 100,
       units = "mm",dpi = 300, device = "pdf",
       path = "/work/users/r/o/roseg/IDDRC/IDDRCPlots/CrossMeasurements")

#Check if DUSP6 and NPCG2 correlate
Prop.Area.Result <- cor.test(D14.check$Growth1435, D14.check$`Dividing neural progenitor cells, S`, method="pearson")
mainpoint <- paste("r =",round(Prop.Area.Result$estimate, digits = 2),
                   "p-value = ",round(Prop.Area.Result$p.value, digits = 6))

D84toPRbyExpt.p <- (ggplot(D14.check, aes(x=Growth1435, y=`Dividing neural progenitor cells, S`))) +
  scale_color_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  geom_smooth(method = lm, color = "black") +
  geom_point(size = 2, aes(col = Site, shape = Replicate)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw()+
  labs(title = "Day 14 to 35 Average Cross-sectional Area Growth and NPC in S Phase at Day 84",
       subtitle = mainpoint)+
  xlab("Day 14 to 35 Growth in Average Cross-sectional Area") +
  ylab("Proportion NPC in S Phase")
D84toPRbyExpt.p
ggsave("Revisions_AllCluster_AllAssay_D84_reordered_D1435Area_NPCinS.pdf",plot = last_plot(), width = 100, height = 100,
       units = "mm",dpi = 300, device = "pdf",
       path = "/work/users/r/o/roseg/IDDRC/IDDRCPlots/CrossMeasurements")
