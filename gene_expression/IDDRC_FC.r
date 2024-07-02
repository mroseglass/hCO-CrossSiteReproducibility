library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(car)

# Inputs ############################################################################
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"

IDDRC <- readRDS(paste0(db.here,"IDDRC_culture_FC.rds"))
Ranks <- readRDS(paste0(db.here,"Rankings_uniquehCO_random_recoded.rds"))
D7Area <- readRDS(paste0(rds.here, "IDDRCD7AreaReady.rds"))
D14Area <- readRDS(paste0(rds.here, "IDDRCD14AreaReady.rds"))
D35Area <- readRDS(paste0(rds.here, "IDDRCD35AreaReady.rds"))
D56Area <- readRDS(paste0(rds.here, "IDDRCD56AreaReady.rds"))

FC <- read_csv(paste0(db.here, "FlowPercentsCompiled.csv"))
FC$Rep <- gsub("R","Rep",FC$Rep)
FC <- mutate(FC, Experiment.name = paste(Site, Rep, sep="_"))

# OUTPUTS ##############################################################################
plot.here <- "/work/users/r/o/roseg/IDDRC/IDDRCPlots/"

# FC Overview #####################################################################
Site.color <-c("#12783D", "#882155", "#332F85")
Site.color <-c("#12783D", "#332F85")

Replicate.shapes <- c(1,2,3,4,15,16,17)

marker.list <- c(unique(FC$Marker))

# remove all CN
FC <- filter(FC, Site != "CN")
#remove late CN
FC <- filter(FC, Experiment != "CN_Rep7")
FC <- filter(FC, Experiment != "CN_Rep6")
FC <- filter(FC, Experiment != "CN_Rep5")
#

pdf(paste(plot.here,"May11_FC_Overview_NoCN_recolor",".pdf"), onefile = TRUE, width = 12, height = 8)
for(i in 1:length(marker.list)){
  FC.m <- filter(FC, Marker == marker.list[i])
  
  if (length(unique(FC.m$Day)) > 1){ 
  FC.m.longit.p <- ggplot(data=FC.m, aes(x=Day, y=Value, group=Experiment.name)) +
  geom_line(aes(color=Site, size = 1)) +
    scale_color_manual(values = Site.color) +
  scale_shape_manual(values=c(1,2,3,4,15,16,17))+
  geom_point(aes(shape = Rep),size = 5) +
  ggtitle(marker.list[i])+
    theme_bw(base_size = 18) +
  ylab(paste("Percent Cell Expressing ",marker.list[i]))
  print(FC.m.longit.p)
    
  FC.m.d <- c(unique(FC.m$Day))
  for(x in 1:length(FC.m.d)){ 
    FC.m.1d <- filter(FC.m, Day == FC.m.d[x])
   #correct 
    FC.ANOVA <- aov(FC.m.1d$Value~Site, data = FC.m.1d)
    FC.ANOVA <- summary(FC.ANOVA)
    Anova.result <- unlist(FC.ANOVA)["Pr(>F)1"]
    
    FC.m.onetime.p <- (ggplot(FC.m.1d, aes(x=Site, y=Value, fill=Site))) +
    geom_boxplot(outlier.shape = NA, lwd=0.2) + 
    scale_shape_manual(values=c(1,2,3,4,15,16,17))+
    scale_fill_manual(values = Site.color) +
    geom_point(aes(shape = Rep),size = 5, position=position_jitterdodge(),alpha=0.9) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme_bw(base_size = 18)+
    ggtitle(paste(marker.list[i],"Day",FC.m.d[x], "aov pval", Anova.result)) + 
    ylab(paste("Percent Cell Expressing ",marker.list[i])) 
  print(FC.m.onetime.p)
  }
    }  
else {
  FC.ANOVA <- aov(FC.m$Value~Site, data = FC.m)
  FC.ANOVA <- summary(FC.ANOVA)
  Anova.result <- unlist(FC.ANOVA)["Pr(>F)1"]
  
  FC.m.onetime.p <- (ggplot(FC.m, aes(x=Site, y=Value, fill=Site))) +
    geom_boxplot(outlier.shape = NA, lwd=0.2) + 
    scale_shape_manual(values=c(1,2,3,4,15,16,17))+
    scale_fill_manual(values = Site.color) +
    geom_point(aes(shape = Rep),size = 5, position=position_jitterdodge(),alpha=0.9) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggtitle(paste(marker.list[i],"aov pval", Anova.result)) + 
    theme_bw(base_size = 18) +
    ylab(paste("Percent Cell Expressing ",marker.list[i]))
    print(FC.m.onetime.p)
    }
}
dev.off()

#Add significance testing across sites for each plot
#make points larger
#label UNC Cysts replicates RED

# FC and RANK ############################################################
#Test to see if marker expression at different times regardless of site correlates with

Ranks <- mutate(Ranks, Experiment.name = gsub("\\s+", "_", Experiment)) 
#add id with Experiment and Day
Ranks <- mutate(Ranks, ExperimentDay = paste0(Experiment,Day))
FC.Rank <- inner_join(IDDRC,Ranks, by = "Experiment.name")
FC.Rank <- dplyr::select(FC.Rank, !c(Site.x))
FC.Rank <- dplyr::rename(FC.Rank, Site = Site.y)
FC.Rank <- dplyr::select(FC.Rank, !c(Experiment.x))
FC.Rank <- dplyr::rename(FC.Rank, Experiment = Experiment.y)

FC.Rank.NoCN <- filter(FC.Rank, Site != "CN")

library(lme4)
Bud <- glmer(Bud ~  `SOX2+PAX6+_7` + as.factor(Ranker) +
               (1|Experiment.name), data = FC.Rank, family = binomial)
Bud.sum <- summary(Bud)
Bud.sum$coefficients[,4] 

#Just D35 
FC.Rank.D35 <- filter(FC.Rank, Day == "D35" )
FC.Rank.D56 <- filter(FC.Rank, Day != "D35" )

Cysts <- glmer(Cysts ~ `SOX2+PAX6+_35` + as.factor(Ranker) + 
               (1|Experiment.name), data = FC.Rank.D56, family = binomial)
Cysts.sum <- summary(Cysts)

FC.Rank.D56 <- dplyr::select(FC.Rank.D56, c(Ranker,Experiment.name,Site,`SOX2+PAX6+_35`,Nest))
FC.Rank.D56 <- distinct(FC.Rank.D56)

Nest <- glmer(Nest ~ `SOX2+PAX6+_35` + as.factor(Ranker) + (1|Experiment.name),
              data = FC.Rank.D56, family = binomial)
Nest.sum <- summary(Nest)


#percents broken down by day for graphing
RankbyExperiment <- Ranks %>%
  dplyr::group_by(ExperimentDay) %>%
  dplyr::summarise(PercentBud = sum(Bud, na.rm=TRUE)/length(Bud),
                   Bud1 = sum(Bud, na.rm=TRUE),
                   BudTot = length(Bud),
                   PercentCysts = sum(Cysts, na.rm=TRUE)/length(Cysts),
                   CystsTot = length(Cysts),
                   PercentNest = sum(Nest, na.rm=TRUE)/length(Nest),
                   NestTot = length(Nest))

IDConverter <- dplyr::select(Ranks, c(Experiment.name,ExperimentDay,Day))
RankbyExperiment <- distinct(inner_join(RankbyExperiment,IDConverter))

FC <- mutate(FC, Assay = paste(Marker,Day))
FC <- dplyr::select(FC, !c(Day))
FC <- FC %>%
  pivot_wider(names_from = "Assay", values_from = "Value")

FCR <- inner_join(RankbyExperiment,FC, by = "Experiment.name")
FCR <- distinct(FCR)

FCR.test <- distinct(dplyr::select(FCR, c(ExperimentDay,`FOXG1 35`,PercentCysts,Site,Rep,Day)))
FCR.test <- filter(FCR.test, PercentCysts > 0)
FCR.test <- filter(FCR.test, Day == "D56")
FCR.test <- FCR.test %>% drop_na()
# plotting any points
FCR.p <- (ggplot(FCR.test, aes(x=`FOXG1 35`, y=PercentCysts)) +
  scale_shape_manual(values=c(1,2,3,4,15,16,17))+
  geom_point(aes(shape = Rep, col = Site, size = 5)) +
  geom_smooth(method = lm)) + theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("CystsD56 ~ `FOXG1 35` + Ranker +(1|Experiment)") +
  ylab("Percent without Cysts") +
  theme(aspect.ratio=1) +theme_bw()
FCR.p

#plot FOXG1 & SOX2PAX6
FCR.foxg <- distinct(dplyr::select(FCR, c(ExperimentDay,`FOXG1 35`,Site,Rep,Day)))
FCR.soxpax <- distinct(dplyr::select(FCR, c(ExperimentDay,`SOX2+PAX6+ 35`,Site,Rep,Day)))
FCR.together <- inner_join(FCR.foxg, FCR.soxpax)
FCR.together <- filter(FCR.together, Day == "D14")
FCR.together <- FCR.together %>% drop_na()

result <- summary(lm(FCR.together$`FOXG1 35` ~ FCR.together$`SOX2+PAX6+ 35`))

FCR.p <- (ggplot(FCR.together, aes(x=`FOXG1 35`, y=`SOX2+PAX6+ 35`)) +
            scale_shape_manual(values=c(1,2,3,4,15,16,17))+
            geom_point(aes(shape = Rep, col = Site, size = 5)) +
            geom_smooth(method = lm)) + theme(axis.text.x = element_text(angle = 45)) +
  ggtitle("Adjust r-squared = 0.4, p = 0.02") +
  #ylab("Percent without Cysts") +
  theme(aspect.ratio=1) +theme_bw()
FCR.p

# FC and Area #################################################
Area <- D35Area
Area <- mutate(Area, Experiment.name = paste(Site, Replicate, sep="_"))

FC <- mutate(FC, Assay = paste(Marker,Day))
FC <- dplyr::select(FC, !c(Day))
FC <- FC %>%
  pivot_wider(names_from = "Assay", values_from = "Value")

FCAssay <- c(unique(colnames(FC)))
#5 to 20 is assay

FCA <- inner_join(Area,FC, by = "Experiment.name")
FCA <- distinct(FCA)
FCA <- dplyr:: rename(FCA, Site = Site.x)

#strugglint to incorporate fixed and random effects in a for loop so
#make a list of functions to loop through in 1 datafram

#fit 1 including Site
fit.with.Site <- list(lmer(AreaMM2 ~  `TBR2 35`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `TBR2 14`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `TBR2 7`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `TBR1 35`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `TBR1 14`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `FOXG1 35`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `FOXG1 14`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `FOXG1 7`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `GSX2 35`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `GSX2 14`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `GSX2 7`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `SOX2+PAX6+ 35`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `SOX2+PAX6+ 14`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `SOX2+PAX6+ 7`+ as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `SSEA3+SSEA4+ -1` + as.factor(Site) + (1|Experiment.name), data = FCA),
                      lmer(AreaMM2 ~  `TRA backbone -1` + as.factor(Site) + (1|Experiment.name), data = FCA))
FwS.out <- list()

for(i in 1:length(fit.with.Site)){
  FwS.out[[i]] <- lmer(fit.with.Site[[i]], data = FCA)
  names(FwS.out)[i] <- capture.output(fit.with.Site[[i]]) # name each entry in output list for each identification
}

#fit 2 excluding site
fit.without.site <- list(lmer(AreaMM2 ~  `TBR2 35`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `TBR2 14`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `TBR2 7`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `TBR1 35`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `TBR1 14`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `FOXG1 35`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `FOXG1 14`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `FOXG1 7`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `GSX2 35`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `GSX2 14`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `GSX2 7`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `SOX2+PAX6+ 35`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `SOX2+PAX6+ 14`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `SOX2+PAX6+ 7`+ (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `SSEA3+SSEA4+ -1` + (1|Experiment.name), data = FCA),
                         lmer(AreaMM2 ~  `TRA backbone -1` + (1|Experiment.name), data = FCA))
FwoS.out <- list()

for(i in 1:length(fit.without.site)){
  FwoS.out[[i]] <- lmer(fit.without.site[[i]], data = FCA)
  names(FwoS.out)[i] <- capture.output(fit.without.site[[i]]) # name each entry in output list for each identification
}

#then compare the 2 models
# Q: is Site and FC assay predict AreaMM2 and Day x

site.FC.anova <- list()
for(i in 1:length(FwS.out)){
  site.FC.anova[[i]] <- anova(unlist(FwS.out[i]),unlist(FwoS.out[i]))
  names(site.FC.anova)[i] <- capture.output(site.FC.anova[[i]])
}

anova(fit.with.Site[[1]],fit.without.site[[1]])


#example
fit1 <- lmer(AreaMM2 ~  `TBR1 35`+ 
                     (1|Experiment.name), data = FCA)
fit2 <- lmer(AreaMM2 ~  `TBR1 35`+ as.factor(Site.x) +
                     (1|Experiment.name), data = FCA)

anova(fit1,fit2)


