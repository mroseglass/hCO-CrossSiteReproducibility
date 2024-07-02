library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(car)
library("ggsignif")

# Inputs ############################################################################
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"
FC <- read_csv(paste0(db.here, "FlowPercentsCompiled.csv"))
FC$Rep <- gsub("R","Rep",FC$Rep)
FC <- mutate(FC, Experiment.name = paste(Site, Rep, sep="_"))

# OUTPUTS ##############################################################################
plot.here <- "/work/users/r/o/roseg/IDDRC/IDDRCPlots/"

# GLOBALS ##########################################################################################
Site.color <-c("#12783D", "#882155", "#4040C6")
Replicate.shapes <- c(1,2,3,4,15,16,17)

#remove samples not in full experiments
FC <- filter(FC, Experiment.name != "UNC_Rep1")
FC <- filter(FC, Experiment.name != "UNC_Rep6")

FC.w <- pivot_wider(FC, names_from  = "Marker", values_from = "Value")

# iPSC Plot #####################################################################
#Plot Percent SSEA
FC.w.0 <- filter(FC.w, Day == "-1")
SiteANOVA <- aov(`SSEA3+SSEA4+`~Site, data = FC.w.0)
SiteANOVA <- summary(SiteANOVA)

FC0.all <- (ggplot(FC.w.0,aes(x=Site, y=`SSEA3+SSEA4+`, fill=Site)))+
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  geom_point(aes(shape = Rep),size = 2,
             position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(101,101,103)) +
    ylab("Percent Population")+ ylim(85,105)+
  ggtitle(paste("SSEA ANOVA pval=",round(SiteANOVA[[1]]$`Pr(>F)`, digits = 3)))

FC0.all

#Plot Percent TRA ###################################
FC.w.0 <- filter(FC.w, Day == "-1")
SiteANOVA <- aov(`TRA backbone`~Site, data = FC.w.0)
SiteANOVA <- summary(SiteANOVA)

FC0.all <- (ggplot(FC.w.0,aes(x=Site, y=`TRA backbone`, fill=Site)))+
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  geom_point(aes(shape = Rep),size = 2,
             position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
    geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(101,101,103)) +
  ylab("Percent Population")+ ylim(85,105)+
  ggtitle(paste("TRA ANOVA pval=",round(SiteANOVA[[1]]$`Pr(>F)`, digits = 3)))

FC0.all

# D35 ##########################################################################################
#Plot Percent `SOX2+PAX6+`
FC.w.35 <- filter(FC.w, Day == "35")
SiteANOVA <- aov(`SOX2+PAX6+`~Site, data = FC.w.35)
SiteANOVA <- summary(SiteANOVA)

FC0.all <- (ggplot(FC.w.35,aes(x=Site, y=`SOX2+PAX6+`, fill=Site)))+
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  geom_point(aes(shape = Rep),size = 2,
             position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(44,44,47)) +
  ylab("Percent Population")+ ylim(0,50)+
  ggtitle(paste("D35 SOX2+PAX6+ ANOVA pval=",round(SiteANOVA[[1]]$`Pr(>F)`, digits = 3)))

FC0.all

#Plot Percent FOXG1
FC.w.35 <- filter(FC.w, Day == "35")
SiteANOVA <- aov(FOXG1~Site, data = FC.w.35)
SiteANOVA <- summary(SiteANOVA)

FC0.all <- (ggplot(FC.w.35,aes(x=Site, y=FOXG1, fill=Site)))+
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  geom_point(aes(shape = Rep),size = 2,
             position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(30,30,35)) +
  ylab("Percent Population")+ ylim(0,50)+
  ggtitle(paste("D35 FOXG1 ANOVA pval=",round(SiteANOVA[[1]]$`Pr(>F)`, digits = 3)))

FC0.all

#Plot Percent GSX2
FC.w.35 <- filter(FC.w, Day == "35")
SiteANOVA <- aov(GSX2~Site, data = FC.w.35)
SiteANOVA <- summary(SiteANOVA)

FC0.all <- (ggplot(FC.w.35,aes(x=Site, y=GSX2, fill=Site)))+
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  geom_point(aes(shape = Rep),size = 2,
             position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(40,40,45)) +
  ylab("Percent Population")+ ylim(0,50)+
  ggtitle(paste("D35 GSX2 ANOVA pval=",round(SiteANOVA[[1]]$`Pr(>F)`, digits = 3)))

FC0.all

#Plot Percent TBR1
FC.w.35 <- filter(FC.w, Day == "35")
SiteANOVA <- aov(TBR1~Site, data = FC.w.35)
SiteANOVA <- summary(SiteANOVA)

FC0.all <- (ggplot(FC.w.35,aes(x=Site, y=TBR1, fill=Site)))+
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  geom_point(aes(shape = Rep),size = 2,
             position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(12,12,17)) +
  ylab("Percent Population")+ ylim(0,50)+
  ggtitle(paste("D35 TBR1 ANOVA pval=",round(SiteANOVA[[1]]$`Pr(>F)`, digits = 3)))

FC0.all

#Plot Percent TBR2
FC.w.35 <- filter(FC.w, Day == "35")
SiteANOVA <- aov(TBR2~Site, data = FC.w.35)
SiteANOVA <- summary(SiteANOVA)

FC0.all <- (ggplot(FC.w.35,aes(x=Site, y=TBR2, fill=Site)))+
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+
  geom_point(aes(shape = Rep),size = 2,
             position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(10,10,15)) +
  ylab("Percent Population")+ ylim(0,50)+
  ggtitle(paste("D35 TBR2 ANOVA pval=",round(SiteANOVA[[1]]$`Pr(>F)`, digits = 3)))

FC0.all

# All D35 Together ####  ##################################################
FC.35 <- filter(FC, Day =="35")
FC.35 <- mutate(FC.35, ExptM = paste(Marker,Site))
#convert 'Marker' to factor and specify level order
FC.35$ExptM <- factor(FC.35$ExptM, 
                       levels=c('FOXG1 CHOP','FOXG1 CN','FOXG1 UNC',
                                'SOX2+PAX6+ CHOP','SOX2+PAX6+ CN','SOX2+PAX6+ UNC',
                                'GSX2 CHOP','GSX2 CN','GSX2 UNC',
                                'TBR1 CHOP','TBR1 CN','TBR1 UNC',
                                'TBR2 CHOP','TBR2 CN','TBR2 UNC'))
gene.p <- (ggplot(FC.35, aes(x=ExptM, y=Value, fill = Site))) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+ 
  geom_point(aes(shape = Rep),size = 2, position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  ylab("Percent Population") 
print(gene.p)



# All iPSC Together ####  ##################################################
FC.0 <- filter(FC, Day =="-1")
FC.0 <- mutate(FC.0, ExptM = paste(Marker,Site))
#convert 'Marker' to factor and specify level order
FC.0$ExptM <- factor(FC.0$ExptM, 
                      levels=c('SSEA3+SSEA4+ CHOP','SSEA3+SSEA4+ CN','SSEA3+SSEA4+ UNC',
                               'TRA backbone CHOP','TRA backbone CN','TRA backbone UNC'))
gene.p <- (ggplot(FC.0, aes(x=ExptM, y=Value, fill = Site))) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+ 
  geom_point(aes(shape = Rep),size = 2, position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
  ylab("Percent Population") 
print(gene.p)
