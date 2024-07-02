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
D14Prop <- read_csv("/work/users/r/o/roseg/IDDRC/IDDRCDatabase/CellTypeProportions/D14ProportionsbyExperiment.csv")
IDconverter <- read_csv(paste0(db.here,"ML_D14scRNAseq_IDconverter.csv"))
D14Prop <- dplyr::rename(D14Prop, "Experiment.scRNAseq" = "Experiment")
D14Prop <- inner_join(D14Prop,IDconverter)
cluster.list <- as_tibble(colnames(D14Prop))

#add technical with left joins
technical <- read_csv(paste0(db.here, "scRNAseqD84_all_technical.csv"))
technical <- mutate(technical, Experiment = gsub("_","",technical$Experiment.name))
D14Prop <- left_join(D14Prop, technical)

#D14Prop <- dplyr::rename(D14Prop, "Experiment.scRNAseq" = "Experiment")
D14Gruffi <- read_csv("/work/users/r/o/roseg/single-cell_reproducibility/data/June16.Gruffi.Suggested.Threshold.IDDRC.Stress.Experiment.Day.Cluster.csv")
D14Gruffi <- mutate(D14Gruffi, Site = gsub("_.*","",D14Gruffi$ExpDayCluster))
D14Gruffi <- mutate(D14Gruffi, Day = sub('.*D','D',D14Gruffi$ExpDayCluster))
D14Gruffi <- mutate(D14Gruffi, Day = gsub("_.*","",D14Gruffi$Day))
D14Gruffi <- mutate(D14Gruffi, Rep = sub('.*R','R',D14Gruffi$ExpDayCluster))
D14Gruffi <- mutate(D14Gruffi, Cluster = sub('.*!','',D14Gruffi$Rep))
D14Gruffi <- mutate(D14Gruffi, Rep = sub('!.*','',D14Gruffi$Rep))
D14Gruffi <- mutate(D14Gruffi, Experiment.scRNAseq = paste(Site,Day,Rep, sep = "_"))
D14Gruffi$Experiment.scRNAseq <- gsub("UNC_D14_P1","P1", D14Gruffi$Experiment.scRNAseq)
D14Gruffi$Experiment.scRNAseq <- gsub("UNC_D84_P1","P1", D14Gruffi$Experiment.scRNAseq)
D14Gruffi <- mutate(D14Gruffi, PcStressed = as.numeric(sub('%','',D14Gruffi$PcStressed))/100)
#pick day
D14Gruffi <- filter(D14Gruffi, Day == "D14")
Stressed <- dplyr::select(D14Gruffi, c("Experiment.scRNAseq","PcStressed","Cluster"))
Stressed <- Stressed %>%
  pivot_wider(names_from = Cluster, values_from = "PcStressed", values_fill = 0)
Stressed <- dplyr::select(Stressed, "Experiment.scRNAseq","9","11")
Stressed <- dplyr::rename(Stressed, `Stressed unspecified neuron`=`9`)
Stressed <- dplyr::rename(Stressed, `Stressed radial glia`=`11`)
D14Prop <- inner_join(D14Prop,Stressed, by = "Experiment.scRNAseq")
#write_csv(D14Prop, paste0(db.here,"CellTypeProportions/Day14CellTypePropwithGRUFFI.csv"))

FC <- read_csv(paste0(db.here,"FC.Assay.by.Experiment.csv"))
FC <- dplyr::rename(FC, "Experiment" = "Experiment.name")
FC$Experiment <- gsub("_","",FC$Experiment,)

#remove D7 and D14 FC
FC <- dplyr::select(FC, c(Experiment,`TBR1 35`,`TBR2 35`,`SOX2+PAX6+ 35`,`SSEA3+SSEA4+ -1`,`TRA backbone -1`,`FOXG1 35`,`GSX2 35`))

Ranks <- read_csv(paste0(db.here,"Percent.Rank.by.Experiment.csv"))
Ranks <- dplyr::rename(Ranks, "Experiment" = "Experiment.name")
Ranks$Experiment <- gsub("_","",Ranks$Experiment,)
qPCR <- read_csv(paste0(db.here, "SY_IDDRC_qPCR.csv"))
qPCR$Experiment <- gsub("_","",qPCR$Experiment)

D14iDISCO <- read_csv(paste0(db.here,"iDISCO/D14_IMARIS_Combined.csv"))
D14iDISCO <- mutate(D14iDISCO, Experiment = paste0(Site,Rep))
D14iDISCO <- D14iDISCO %>%
  group_by(Experiment) %>%
  dplyr::summarise(PerPAX6 = mean(percentpax6, na.rm = TRUE),
                   Volume = mean(ToPROVolume, na.rm = TRUE),
                   NCADperToPRO = mean(NCADperTotalVol, nam.rm=TRUE),
                   NCADperPAX6 = mean(NCADperPAX6, na.rm=TRUE))


D14.Prop.FC.Rank <- full_join(D14Prop,IDconverter)
D14.Prop.FC.Rank <- left_join(D14.Prop.FC.Rank,qPCR, by = "Experiment")
D14.Prop.FC.Rank <- left_join(D14.Prop.FC.Rank,Ranks, by = "Experiment")
D14.Prop.FC.Rank <- left_join(D14.Prop.FC.Rank,FC, by = "Experiment")
D14.Prop.FC.Rank <- left_join(D14.Prop.FC.Rank,D14iDISCO, by = "Experiment")

#add it area measurements
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
D14.Prop.FC.Rank <- left_join(D14.Prop.FC.Rank,AreaGrowthMeasures)
D14.Prop.FC.Rank <- left_join(D14.Prop.FC.Rank,area.scrnaseq)
D14.Prop.FC.Rank <- dplyr::select(D14.Prop.FC.Rank, !c(scRNAseqAreaD84))

#unique(D14.Prop.FC.Rank$Experiment)

#use lm to see if any biological variables correlated to D84 proportion
Prop.Area <- matrix()
Prop.Area.Result <- list()

# FIRST EXAMPLE FOR CLUSTER TO ASSAY ####################################################
#partial correlation because only testing 2 variables
#create new data fram with one cluster plus all assays 
#create new data fram with one cluster plus all assays 
cluster.list2 <- as_tibble(colnames(Stressed))
#add 2 cluster lists
cluster.list <- rbind(cluster.list,cluster.list2)
assay.list <- as_tibble(colnames(D14.Prop.FC.Rank))
assay.list <- anti_join(assay.list,cluster.list)
assay.list <- assay.list$value
assay.list <- assay.list[-1] #to remove NameofCookie
assay.list <- assay.list[-13] #to remove Experiment.name
assay.list <- assay.list[-13] #to remove hCOAge
cluster.list <- cluster.list$value
cluster.list <- cluster.list[-16] #to remove Experiment.scRNAseq
cluster.list <- cluster.list[-16] #to remove Experiment.scRNAseq
cluster.list <- cluster.list[-16] #to remove Experiment

# pull one cluster to run p-cor with 
y=1
Test <- D14.Prop.FC.Rank %>% dplyr::select(all_of(c(assay.list,cluster.list[y])))
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
assay.D14.prop.sum <- ready

#Then loop for all other cluster ############################################################
for(y in 2:length(cluster.list)){
Test <- D14.Prop.FC.Rank %>% dplyr::select(all_of(c(assay.list,cluster.list[y])))
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

#write.csv(assay.D14.prop.sum, paste0(db.here,"Revisions.D14.Proportions.GRUFFI.to.Assays.Technical.Dec28.r.Pearsons.csv"))

# CORRELATION FOLLOW UP ########################################################################################
D14.check <- mutate(D14.Prop.FC.Rank, Replicate = gsub(".*_","",D14.Prop.FC.Rank$Experiment.name))
D14.check <- mutate(D14.check, Site = gsub("_.*","",D14.check$Experiment.name))

# GLOBALS ##########################################################################################
Site.color <-c("#12783D", "#882155", "#332F85")
Replicate.shapes <- c(1,2,3,4,15,16,17)

#Check if DUSP6 and ZIC2 correlate
Prop.Area.Result <- cor.test(D14.check$DUSP6_deltadeltaCT, D14.check$ZIC2_deltadeltaCT, method="pearson")
mainpoint <- paste("r =",round(Prop.Area.Result$estimate, digits = 2),
                   "p-value = ",round(Prop.Area.Result$p.value, digits = 6))

D84toPRbyExpt.p <- (ggplot(D14.check, aes(x=DUSP6_deltadeltaCT, y=ZIC2_deltadeltaCT))) +
  scale_color_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  geom_smooth(method = lm, color = "black") +
  geom_point(size = 2, aes(col = Site, shape = Replicate)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw()+
  labs(title = "DUSP6 and ZIC2 in iPSCs",
       subtitle = mainpoint)+
  xlab("DUSP6 Delta CT") +
  ylab("ZIC2 Delta CT")
D84toPRbyExpt.p

#Check if LLNb and ZIC2 correlate
Prop.Area.Result <- cor.test(D14.check$ZIC2_deltadeltaCT, D14.check$`Lower layer neuron II`, method="pearson")
mainpoint <- paste("r =",round(Prop.Area.Result$estimate, digits = 2),
                   "p-value = ",round(Prop.Area.Result$p.value, digits = 6))

D84toPRbyExpt.p <- (ggplot(D14.check, aes(x=ZIC2_deltadeltaCT, y=`Lower layer neuron II`))) +
  scale_color_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  geom_smooth(method = lm, color = "black") +
  geom_point(size = 2, aes(col = Site, shape = Replicate)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw()+
  labs(title = "iPSC ZIC2 and Lower layer neuron b at Day 14",
       subtitle = mainpoint)+
  xlab("ZIC2 Delta CT") +
  ylab("Lower layer neuron b")
D84toPRbyExpt.p

#Check if DUSP6 and NPCG2 correlate
Prop.Area.Result <- cor.test(D14.check$DUSP6_deltadeltaCT, D14.check$`Dividing neural progenitor cells, G2`, method="pearson")
mainpoint <- paste("r =",round(Prop.Area.Result$estimate, digits = 2),
                   "p-value = ",round(Prop.Area.Result$p.value, digits = 6))

D84toPRbyExpt.p <- (ggplot(D14.check, aes(x=DUSP6_deltadeltaCT, y=`Dividing neural progenitor cells, G2`))) +
  scale_color_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  geom_smooth(method = lm, color = "black") +
  geom_point(size = 2, aes(col = Site, shape = Replicate)) +
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw()+
  labs(title = "iPSC DSUP6 and NPC, G2 Phase at Day 14",
       subtitle = mainpoint)+
  xlab("DUSP6 Delta CT") +
  ylab("NPC, G2 Phase")
D84toPRbyExpt.p
