library(lme4)
library(lmerTest)
library(tidyverse)

setwd("/work/users/r/o/roseg")
HERE <- "/work/users/r/o/roseg/single-cell_reproducibility/"
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"

# INPUTS #####################################################################################
qPCR <- read_csv(paste0(db.here, "SY_IDDRC_qPCR.csv"))

# GLOBALS ##########################################################################################
#would need to plot each whisker separatel and then merge
Site.color <-c("#12783D", "#882155", "#4040C6")
Replicate.shapes <- c(1,2,3,4,15,16,17)

#remove samples not in full experiments
qPCR <- filter(qPCR, Experiment != "UNC_Rep1")
qPCR <- filter(qPCR, Experiment != "UNC_Rep6")
qPCR <- filter(qPCR, Experiment != "UNC_Rep7")

#Add site and replicate
qPCR <- mutate(qPCR, Site = sub("_.*", "", qPCR$Experiment))
qPCR <- mutate(qPCR, Replicate = sub(".*_", "", qPCR$Experiment))


# PLOTS ###########################################################################
#Statistical Testing for Site Differences
SiteANOVA <- aov(OCT4_deltadeltaCT~Site, data = qPCR)
SiteANOVA <- summary(SiteANOVA)
#save p-val for FDR correction
OCT4 <- 0.0862

#OCT4
qPRC.p <- (ggplot(qPCR, aes(x=Site, y=OCT4_deltadeltaCT, fill=Site))) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  geom_point(aes(shape = Replicate),size = 2, position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  ylab("Log10 Expression")+ scale_y_continuous(trans='log10') +
  ylim(0,5.5) + xlab("OCT4")
qPRC.p

#Nanog
#Statistical Testing for Site Differences
SiteANOVA <- aov(Nanog_deltadeltaCT~Site, data = qPCR)
SiteANOVA <- summary(SiteANOVA)
Nanog <- 0.0174

qPRC.p <- (ggplot(qPCR, aes(x=Site, y=OTX2_deltadeltaCT, fill=Site))) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  #geom_jitter(width = 0.1)+
  geom_point(aes(shape = Replicate),size = 2, position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(1,1,1.3))+
              ylab("Log Expression")+ scale_y_continuous(trans='log10') +
  ylim(0,5.5) + ylab("Log10 Expression") + xlab("NANOG")
qPRC.p

#Statistical Testing for Site Differences
SiteANOVA <- aov(DUSP6_deltadeltaCT~Site, data = qPCR)
SiteANOVA <- summary(SiteANOVA)
DUSP6 <- 0.0742
#DUSP6
qPRC.p <- (ggplot(qPCR, aes(x=Site, y=DUSP6_deltadeltaCT, fill=Site))) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  #geom_jitter(width = 0.1)+
  geom_point(aes(shape = Replicate),size = 2, position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  ylab("Log Expression")+ scale_y_continuous(trans='log10') +
  ylim(0,1.5) + ylab("DUSP6") 
qPRC.p

#ZIC2
#Statistical Testing for Site Differences
SiteANOVA <- aov(ZIC2_deltadeltaCT~Site, data = qPCR)
SiteANOVA <- summary(SiteANOVA)
ZIC2 <-0.832

qPRC.p <- (ggplot(qPCR, aes(x=Site, y=ZIC2_deltadeltaCT, fill=Site))) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  #geom_jitter(width = 0.1)+
  geom_point(aes(shape = Replicate),size = 2, position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  ylab("Log Expression")+ scale_y_continuous(trans='log10') +
  ylim(0,1.5) + ylab("ZIC2")
qPRC.p


# EXAMPLE PLOTS #################################################
# All Markers ################################################################################
#make one plot with all markers geom_points still group by markers and not by Site
qPCR.l <- qPCR %>% pivot_longer(!Experiment, names_to = "Marker", values_to = "DeltaDeltaCT")
qPCR.l <- mutate(qPCR.l, Site = sub("_.*", "", qPCR.l$Experiment))
qPCR.l <- mutate(qPCR.l, Replicate = sub(".*_", "", qPCR.l$Experiment))

qPCR.l$Marker <-gsub("_.*", "", qPCR.l$Marker)
qPCR.l$Marker <- gsub("Nanog", "NANOG", qPCR.l$Marker)
qPCR.l <- mutate(qPCR.l, ExptM = paste(Marker,Site))
#convert 'Marker' to factor and specify level order
qPCR.l$ExptM <- factor(qPCR.l$ExptM, 
                        levels=c('OCT4 CHOP','OCT4 CN','OCT4 UNC',
                                 'NANOG CHOP','NANOG CN','NANOG UNC', 
                                 'DUSP6 CHOP','DUSP6 CN','DUSP6 UNC',
                                 'OTX2 CHOP','OTX2 CN','OTX2 UNC',
                                 'ZIC2 CHOP','ZIC2 CN','ZIC2 UNC',
                                 'TFAP2C CHOP','TFAP2C CN','TFAP2C UNC',
                                 'KLF4 CHOP','KLF4 CN','KLF4 UNC',
                                 'TFCP2L1 CHOP','TFCP2L1 CN','TFCP2L1 UNC'))


gene.p <- (ggplot(qPCR.l, aes(x=ExptM, y=DeltaDeltaCT, fill = Site))) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+ 
  geom_point(aes(shape = Replicate),size = 2, position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  ylab("Log Expression")+ scale_y_continuous(trans='log10') 
  print(gene.p)



#Statistical Testing for Site Differences #############################################################
SiteANOVA <- aov(OTX2_deltadeltaCT~Site, data = qPCR)
SiteANOVA <- summary(SiteANOVA)
OTX2 <-0.0108

SiteANOVA <- aov(TFAP2C_deltadeltaCT~Site, data = qPCR)
SiteANOVA <- summary(SiteANOVA)
TFAP2C <-0.0545

SiteANOVA <- aov(KLF4_deltadeltaCT~Site, data = qPCR)
SiteANOVA <- summary(SiteANOVA)
KLF4 <-0.262

SiteANOVA <- aov(TFCP2L1_deltadeltaCT~Site, data = qPCR)
SiteANOVA <- summary(SiteANOVA)
TFCP2L1 <- 0.693

FDRtest <- rbind(DUSP6,KLF4,Nanog,OCT4,OTX2,TFAP2C,TFCP2L1,ZIC2)
FDRtest <- as.data.frame(FDRtest)
FDRtest <- mutate(FDRtest, FDRCheck = p.adjust(V1, method = "BH"))


#Statistical Testing for Site Differences
# FDR correct pva ~ 0.00625
SiteANOVA <- aov(TFAP2C_deltadeltaCT~Site, data = qPCR)
SiteANOVA <- summary(SiteANOVA)
SiteANOVA

#nominal sig: Nanog, OTX2,TFAP2C

#save all one-way anovas as .rds
qPCR <- mutate(qPCR, Site = sub("_.*", "", qPCR$Experiment))
qPCR <- mutate(qPCR, Replicate = sub(".*_", "", qPCR$Experiment))
output <- list()

for(i in 2:9){
marker <- cbind(qPCR[10], qPCR[i])
colnames(marker)[2] = "Marker"
SiteANOVA <- aov(Marker~Site, data = marker)
SiteANOVA <- summary(SiteANOVA)
output[[i]] <- SiteANOVA
names(output)[i] <- colnames(qPCR)[i]
}

saveRDS(output, paste0(db.here, "qPCR_OneWayANOVAs.rds"))
