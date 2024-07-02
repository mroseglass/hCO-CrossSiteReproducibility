# Rose plots for iPSC gene expression
library(vctrs)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readxl)
library(cowplot)
library(viridis)
library(ggpubr)
library(forcats)

db.here <- "/work/users/r/o/roseg/IVIVData/InfoVerifiedID/"
qPCR <- "IVIV_KE_qPCR_detlaCT.csv"

qPCR <- read_csv(paste0(db.here,qPCR))
qPCR$target_name <- gsub("TOCT4","OCT4", qPCR$target_name)
qPCR$target_name <- factor(qPCR$target_name,
                           levels = c("EIF4A2", 
                                      "OCT4", "SOX2","NANOG", 
                                      "ZIC2", "DUSP6", "OTX2",
                                      "TFAP2C","KLF4","TFCP2L1",
                                      "PAX6","SOX17","SNAI2")) 
                                 
qPCR %>%
  ggplot(aes(x = target_name, y=deltaCT)) +
  geom_boxplot() +
  geom_point(position=position_dodge(width = .75),aes(col = sample_name)) +
  labs(title = "IPSC CT Means",
       x = "Target Genes",
       y = "CT Means",
       color = "condition")

qPCR.ctfilter <- filter(qPCR, ct_mean <30)
qPCR.ctfilter <- filter(qPCR.ctfilter, target_name != "EIF4A2")
qPCR.ctfilter %>%
  ggplot(aes(x = target_name, y=deltaCT)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.01),aes(col = EBID)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(title = "iPSC Gene Expression Before Seeding",
       x = "Target Genes",
       y = "Delta CT",
       color = "EBID")

qPCR.ctfilter <- dplyr::rename(qPCR.ctfilter, EBID = sample_name)

IVIV <- readRDS(paste0(db.here,"IVIV_IBIS_Culture_Technical_VerifiedIDs.rds"))
outcome.ebid <- distinct(dplyr::select(IVIV, c(EBID, DCCID,VerifiedEBIDMapping,VerifiedBatchIDMapping,VerifiedBuddingIDMapping, OUTCOME)))

qPCR.ctfilter <- inner_join(qPCR.ctfilter, outcome.ebid, by = "EBID")

EBIDyDCCID <- c("Beau","Orange","Orange2","Sara","Sodi","Green2",
"Baor","Violet","McMac","Violet2","Ginger","PA","PA2","Cost","Cher","Leo","Nina","GA","GA2","Cate","IA","Owen",
"Posh","Posh2","Dana","Dana2","FL","Yellow2","Baby","Arizona","Arizona2","Black","Dach","Adam","Telly","Telly2",
"Slol","Mary","Splinter","Splinter2","Sebe","Sebe2","mifi","White","White2","Lucky","Nile","Val","Apple","AppleW",
"Bnan","BnanW","Mich","Mich2","Sporty","Blue","Blue2","Scary","frul","Evan","WI","Raph","Raph2","RaphW","Red")

qPCR.ctfilter %>%
  ggplot(aes(x = target_name, y=deltaCT)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(dodge.width=0.01),aes(col = EBID)) +
  labs(title = "iPSC Gene Expression Before Seeding",
       x = "Target Genes",
       y = "Delta CT",
       color = "EBID")

#taqman is 100 - 189
Taqman <- c(names(IVIV)[105:184])
FoldChangeFilter <- function(x){
  dummyReplace <- NA
  x <- ifelse(x > 50, dummyReplace, x)
  return(x)
}
IVIV <- IVIV %>%
  mutate_at(Taqman, FoldChangeFilter)
Taqman.order <- read_csv("/work/users/r/o/roseg/IVIVData/InfoVerifiedID/IVIV_Taqman_Order.csv")

taqman <- IVIV[100:191]
taqman <- dplyr::select(taqman, !c(Taqman_Ectoderm,Taqman_Endoderm,Taqman_Mesoderm,Taqman_SelfRenewal,SampleSwapDetectediPSCExpansion))
taqman.psc <- filter(taqman, Taqman_Culture == "iPSC")
taqman.tri <- filter(taqman, Taqman_Culture == "Tri")

#make long for geom_tile
taqman.psc <- dplyr::select(taqman.psc, c(!Taqman_Culture))
taqman.psc <- pivot_longer(taqman.psc, cols = !iPSCLineClone, names_to = "Gene", values_to = "logfc")
taqman.psc$Gene <- gsub("Taqman_FoldChange_","",taqman.psc$Gene)
taqman.psc <- inner_join(taqman.psc, Taqman.order, by = "Gene")

#colors <- c("red","blue","green","grey","black")
colors <- top_bot_5_both$color[order(top_bot_5_both$name)]
                                      
taqman.psc.p <- ggplot(taqman.psc, aes(x = iPSCLineClone, y=Gene, fill = logfc)) + 
  geom_tile()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_gradient(name = "logfc", trans = "log",
                        breaks = c(0.1,0.5,1,5,10,20,35,50))+
  labs(title = "Day 84 Cell Type Proportions Correlates to Assays", x = "Assay", y = "Cell Type") 
#+  theme(axis.text.x = element_text(color = taqman.psc$State))
taqman.psc.p



