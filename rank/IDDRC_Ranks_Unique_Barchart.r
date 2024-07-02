library(lme4)
library(lmerTest)
library(tidyverse)
library(patchwork)
library(plotrix)
library(ggsignif)

HERE <- "/work/users/r/o/roseg/single-cell_reproducibility/"
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"

# INPUTS #####################################################################################
Ranks <- read_csv(paste0(db.here,"Percent.Rank.by.Experiment.csv"))
Ranks <- dplyr::rename(Ranks, "Experiment" = "Experiment.name")
Ranks$Experiment <- gsub("_","",Ranks$Experiment)
Ranks <- mutate(Ranks, Site = gsub("Rep","_",Ranks$Experiment))
Ranks <- mutate(Ranks, Site = gsub("_.*","",Ranks$Site))

#remove UNC Cysts
Ranks <- filter(Ranks, Experiment != "UNCRep6")
Ranks <- filter(Ranks, Experiment != "UNCRep1")
Ranks <- filter(Ranks, Experiment != "UNCRep7")

iddrc.database <- readRDS(paste0(db.here,"IDDRC_culture_FC.rds"))

# LM testing with ranker & random effect #
Ranks.lm <- readRDS(paste0(db.here,"Rankings_uniquehCO_random_recoded.rds"))
Ranks.lm <- mutate(Ranks.lm, Experiment.name = gsub("\\s+", "_", Experiment)) 
Ranks.lm <- mutate(Ranks.lm, ExperimentDay = paste0(Experiment,Day))

# GLOBALS ##########################################################################################
Site.color <-c("#12783D", "#882155", "#4040C6")
Replicate.shapes <- c(1,2,3,4,15,16,17)

######################################################################################################
##save statistics
Bud.m1 <- glmer(Bud ~  as.factor(Site) + as.factor(Ranker) +
               (1|Experiment.name), data = Ranks.lm, family = binomial)

Bud.m2 <- glmer(Bud ~  as.factor(Ranker) +
                  (1|Experiment.name), data = Ranks.lm, family = binomial)

Bud.Site <- anova(Bud.m1, Bud.m2)
#site not significant whereas ranker is 

Nest.m1 <- glmer(Nest ~  as.factor(Site) + as.factor(Ranker) +
                (1|Experiment.name), data = Ranks.lm, family = binomial)
Nest.m2 <- glmer(Nest ~  as.factor(Ranker) +
                   (1|Experiment.name), data = Ranks.lm, family = binomial)
Nest.Site <- anova(Nest.m2, Nest.m1)


#instead dataframe filtered to remove one site
Nest.CNvUNC <-filter(Ranks.lm, Site != "CHOP")
Nest.m1 <- glmer(Nest ~  as.factor(Site) + as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CNvUNC, family = binomial)
Nest.m2 <- glmer(Nest ~  as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CNvUNC, family = binomial)
Nest.Site <- anova(Nest.m2, Nest.m1)

Nest.CHOPvUNC <-filter(Ranks.lm, Site != "CN")
Nest.m1 <- glmer(Nest ~  as.factor(Site) + as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CHOPvUNC, family = binomial)
Nest.m2 <- glmer(Nest ~  as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CHOPvUNC, family = binomial)
Nest.Site <- anova(Nest.m2, Nest.m1)

Nest.CNvCHOP <-filter(Ranks.lm, Site != "UNC")
Nest.m1 <- glmer(Nest ~  as.factor(Site) + as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CNvCHOP, family = binomial)
Nest.m2 <- glmer(Nest ~  as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CNvCHOP, family = binomial)
Nest.Site <- anova(Nest.m2, Nest.m1)

#rescale variables
Ranks35 <- filter(Ranks.lm, Day == "D35")
Cysts.m1 <- glmer(Cysts ~  as.factor(Site) + as.factor(Ranker) +
                 (1|Experiment), data = Ranks35,control = glmerControl(optimizer ="Nelder_Mead"),
                 family = binomial)
Cysts.m2 <- glmer(Cysts ~ as.factor(Ranker) +
                    (1|Experiment), data = Ranks35,control = glmerControl(optimizer ="Nelder_Mead"),
                  family = binomial)
Cysts.Site <- anova(Cysts.m1,Cysts.m2)
#0.00465 is site significant in model

#check for site v site differences
Nest.CNvUNC <-filter(Ranks35, Site != "CHOP")

Cysts.m1 <- glmer(Cysts ~  as.factor(Site) + as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CNvUNC, 
                  glmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')),
                  family = binomial)
Cysts.m2 <- glmer(Cysts ~  as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CNvUNC,
                  glmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')),
                  family = binomial)
Nest.Site <- anova(Cysts.m2, Cysts.m1)

Nest.CHOPvUNC <-filter(Ranks35, Site != "CN")
Nest.m1 <- glmer(Cysts ~  as.factor(Site) + as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CHOPvUNC,
                 #glmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')), #fails to converge
                 family = binomial)
Nest.m2 <- glmer(Cysts ~  as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CHOPvUNC,
                 glmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')),
                 family = binomial)
Nest.Site <- anova(Nest.m2, Nest.m1)

Nest.CNvCHOP <-filter(Ranks35, Site != "UNC")
Nest.m1 <- glmer(Cysts ~  as.factor(Site) + as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CNvCHOP,
                 #control = glmerControl(optimizer ="Nelder_Mead"), fails to converge
                 family = binomial)

Nest.m2 <- glmer(Cysts ~  as.factor(Ranker) +
                   (1|Experiment.name), data = Nest.CNvCHOP,control = glmerControl(optimizer ="Nelder_Mead"),
                 family = binomial)
Nest.Site <- anova(Nest.m2, Nest.m1)


#Rank d56
Ranks56 <- filter(Ranks.lm, Day == "D56")
Cysts.m1 <- glmer(Cysts ~  as.factor(Site) + as.factor(Ranker) +
                 (1|Experiment), data = Ranks56, family = binomial)
Cysts.m2 <- glmer(Cysts ~ as.factor(Ranker) +
                 (1|Experiment), data = Ranks56, family = binomial)
Cysts.Site <- anova(Cysts.m1,Cysts.m2)
#Site significant in model p= 0.377

# Make Individual Plots #########################################################################
#Percent Bud
Bud <- Ranks %>%
  ggplot(aes(fill=Site, y=PercentBud, x=Experiment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values= Site.color)+
  ylab("Percent Visible Buds") + ylim(0,1) +
    ggtitle("No Significant Site Differences Percent Bud")

#Percent Growth into matrigel
Nests <- Ranks %>%
  ggplot(aes(fill=Site, y=PercentNest, x=Experiment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values= Site.color)+
  ylab("Percent Growth in Matrigel") + ylim(0,1) +
  ggtitle("CN significantly less Growth in Matrigel")

#Percent smooth D35
#CHOPR6 not Ranked
RanksD35 <- filter(Ranks, Experiment != "CHOPRep6")
D35QC <- RanksD35 %>%
  ggplot(aes(fill=Site, y=PercentCystsD35, x=Experiment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values= Site.color)+
  ylab("Percent Ideal") + ylim(0,1)+
  ggtitle("CN UNC Significantly worse than CHOP")

#Percent Not smooth D56
D56QC <- Ranks %>%
  ggplot(aes(fill=Site, y=PercentCystsD56, x=Experiment)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values= Site.color)+
  ylab("Percent Ideal")+
  ggtitle("No Difference in QC D56")


#AverageArea per well at D56 to Bursting ##########################################
# not significant as anova or glm for D35 or D56

Bursting <- read_csv(paste0(db.here,"CHOPBurstingRank.csv"))
#get NameofCookie to be the same variable
D56AreaA <- readRDS(paste0(rds.here, "IDDRCD56AreaReady.rds"))
D56AreaA <- dplyr::rename(D56AreaA, FileName_Original = FileName_D56Original)
D56Area <- filter(D56AreaA, Site == "CHOP")
#get name of cookie for CHOP
D56Area <- mutate(D56Area, NameofCookie = gsub(".*56-","",FileName_Original))
#D35Area <- mutate(D35Area, NameofCookie = gsub(".*35-","",FileName_Original))
D56Area <- mutate(D56Area, NameofCookie = substr(NameofCookie, 1, 3))
D56Area <- mutate(D56Area, NameofCookie = gsub("-","",NameofCookie))
D56Area <- mutate(D56Area, NameofCookie = gsub("-","",NameofCookie))
D56Area <- mutate(D56Area, NameofCookie = gsub("ep","",paste0(Replicate,NameofCookie)))
D56Area <- mutate(D56Area, NameofCookie = gsub("c","C",NameofCookie))

D56AreabyNoc <- D56Area %>%
  dplyr::group_by(NameofCookie) %>%
  dplyr::summarise(D56Area = mean(AreaMM2, na.rm = TRUE),
                   D56Area.sem = std.error(AreaMM2, na.rm = TRUE),
                   Experiment = Experiment)
#for D56
D56AreabyNoc <- mutate(D56AreabyNoc, Rep = gsub(".*R","",Experiment))
D56AreabyNoc$Rep <-gsub("ep","Rep",D56AreabyNoc$Rep)

Bursting <- mutate(Bursting,NameofCookie = paste0("R",gsub("Rep","",Rep),Cookie))

AreatoBurst <- inner_join(Bursting, D56AreabyNoc, by = "NameofCookie")
AreatoBurst <- distinct(AreatoBurst)

Bursting.p <- (ggplot(AreatoBurst, aes(x=as.factor(Bursting), y=D56Area, fill = as.factor(Bursting)))) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = c("#BCBCBC","#777777")) +
  scale_shape_manual(values=Replicate.shapes)+
  theme_classic()+ 
  geom_point(aes(shape = Rep.x, colour = as.factor(MotorFail)),
             size = 2, position=position_jitterdodge(dodge.width=0.1)) +
  ylab("Area mm2 Day 56")
Bursting.p

#not significantly correlated when controlling for just replicate, need to include motor failure
BvMF.A <- glmer(Bursting ~ D56Area + as.factor(MotorFail)+(1|Rep.x), data = AreatoBurst, family = binomial)
BvMF.A <- summary(BvMF.A)

#BvMF <- glmer(Bursting ~ as.factor(MotorFail)+(1|Rep.x), data = AreatoBurst, family = binomial)
#summary(BvMF)

#BvA <- glmer(Bursting ~  D56Area + (1|Rep.x), data = AreatoBurst, family = binomial)
#summary(BvA)

#BvMF.nr <- aov(Bursting ~ as.factor(MotorFail), data = AreatoBurst)
#summary(BvMF.nr)

#BvA.nr <- aov(Bursting ~ D56Area, data = AreatoBurst)
#summary(BvA.nr)

#ANOVAbursting <- aov(BvMF.A,BvA)
#summary(ANOVAbursting)

#Stacked barchart for bursting
PercentB <- read_csv(paste0(db.here, "SiteBursting.csv"))
BurstingbyBW <- read_csv(paste0(db.here, "CHOPBurstingRank.csv"))
#add statistical test 
#ANOVAbursting <- aov(Bursting~as.factor(Site), data = BurstingbyBW)
#SiteANOVA <- summary(ANOVAbursting)
#use this below 
ANOVA.check <- glm(Bursting~as.factor(Site), family = "binomial",data = BurstingbyBW)
site <-anova(ANOVA.check, test="Chisq")
 
#instead of performing t-test try for the same glm with logistic regression
#but with datafram containing only 2 sites
Bursting.2 <- filter(BurstingbyBW, Site !="UNC")
ANOVA.check <- glm(Bursting~as.factor(Site), family = "binomial",data = Bursting.2)
CNvCHOP <- anova(ANOVA.check, test="Chisq")
Bursting.2 <- filter(BurstingbyBW, Site !="CN")
ANOVA.check <- glm(Bursting~as.factor(Site), family = "binomial",data = Bursting.2)
UNCvCHOP <- anova(ANOVA.check, test="Chisq")
Bursting.2 <- filter(BurstingbyBW, Site !="CHOP")
ANOVA.check <- glm(Bursting~as.factor(Site), family = "binomial",data = Bursting.2)
UNCvCN <- anova(ANOVA.check, test="Chisq")

#
p <- PercentB %>% ggplot(aes(fill=Site, y=PercentBurstingWells, x=Site)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values= Site.color)+
  ylab("Percent Wells Bursting")+ ylim(0,1) +
  labs(title = "Percent Wells with Bursting",
         subtitle = paste0("glm ANOVA p-value = ",round(site$`Pr(>Chi)`[2], digits = 10)),
         caption = paste0("CNvUNC",round(UNCvCN$`Pr(>Chi)`[2],digits = 3),
                          "CHOPvUNC", round(UNCvCHOP$`Pr(>Chi)`[2], digits = 10),
                          "CHOPvCN", round(CNvCHOP$`Pr(>Chi)`[2], digits = 7)))
p
ggsave("June29_BurstingbySite.pdf",plot = last_plot(), width = 66, height = 90,
       units = "mm",dpi = 300, device = "pdf",
       path = "/work/users/r/o/roseg/IDDRC/IDDRCPlots/Rank")


#Stacked barchart for bursting for motor failure
#add statistical test 
BvMF.nr <- aov(Bursting ~ as.factor(MotorFail), data = AreatoBurst)
summary(BvMF.nr)
#
bmf <- dplyr::select(AreatoBurst, c(MotorFail,Bursting))
bmf <- filter(bmf, Bursting == 1)
bmf %>% ggplot(aes(fill=as.factor(MotorFail), x=as.factor(MotorFail))) +
    geom_bar() +
  labs(title ="Motor Failure", subtitle = paste0("p-value = ", round(BvMF.A$coefficients[3,4], digits=4)))
ggsave("June29_BurstingbyCHOPmotorfail.pdf",plot = last_plot(), width = 76, height = 90,
       units = "mm",dpi = 300, device = "pdf",
       path = "/work/users/r/o/roseg/IDDRC/IDDRCPlots/Rank")

