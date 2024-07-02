library("ggsignif")
library(tidyverse)
library(dplyr)
library(readxl)
library(plotrix)
setwd("/work/users/r/o/roseg")

# INPUT FILES ############################################################################################
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"

#Removing unregistered CN images from each day ###########################################################
D7 <- readRDS(paste0(rds.here, "CNRegistered_D7Area.rds"))%>%
  dplyr::select(c(ImageName,AreaMM2,Site,Replicate,Experiment)) %>%
  filter(Site != "CN_First")
Day <- rep("Day 07", length(D7$ImageName))
D7 <- cbind(D7,Day)

D14 <- readRDS(paste0(rds.here, "CNRegistered_D14Area.rds")) %>% dplyr::rename(AreaMM2 = D14AreaMM2) %>% dplyr::rename(ImageName = FileName_D14Original) %>%
  dplyr::select(c(ImageName,AreaMM2,Site,Replicate,Experiment)) %>%
  filter(Site != "CN_First")
Day <- rep("Day 14", length(D14$ImageName))
D14 <- cbind(D14,Day)

D35 <- readRDS(paste0(rds.here, "CNRegistered_D35Area.rds")) %>% dplyr::rename(AreaMM2 = D35AreaMM2) %>% rename(ImageName = FileName_Original)%>%
  dplyr::select(c(ImageName,AreaMM2,Site,Replicate,Experiment))%>%
  filter(Site != "CN_First")
Day <- rep("Day 35", length(D35$ImageName))
D35 <- cbind(D35,Day)

D56 <- readRDS(paste0(rds.here, "CNRegistered_D56Area.rds")) %>% rename(AreaMM2 = D56AreaMM2) %>% rename(ImageName = FileName_D56Original)%>%
  dplyr::select(c(ImageName,AreaMM2,Site,Replicate,Experiment))%>%
  filter(Site != "CN_First")
Day <- rep("Day 56", length(D56$ImageName))
D56 <- cbind(D56,Day)

# OUTPUTS ##############################################################################
plot.here <- "/work/users/r/o/roseg/IDDRC/IDDRCPlots/"

# GLOBALS ##########################################################################################
Site.color <-c("#12783D", "#882155", "#4040C6")
Replicate.shapes <- c(1,2,3,4,15,16,17)

#remove samples not in complete experiments
Area <- rbind(D7,D14,D35,D56)
Area <- filter(Area, Experiment != "UNCRep1") %>% 
  filter(Experiment != "UNCRep6") %>%
  filter(Experiment != "UNCRep7") %>%
  filter(Experiment != "CHOPRep6")

AreabyExperiment <- mutate(Area, Experiment = (paste(Experiment,Day, sep=".")))

AreabyExperiment <- AreabyExperiment %>%
  dplyr::group_by(Experiment) %>%
  dplyr::summarise(Area = mean(AreaMM2, na.rm = TRUE),
                   Area.sem = std.error(AreaMM2, na.rm = TRUE))

AreabyExperiment <- mutate(AreabyExperiment, Day = str_extract(AreabyExperiment$Experiment,"Day [0-9][0-9]")) %>%
  mutate(Site = gsub( "Rep[0-9]", "", Experiment)) %>%
  mutate(Site = gsub( ".Day [0-9][0-9]", "", Site)) %>%
  mutate(Replicate = str_extract(Experiment, "Rep[0-9]")) %>%
  mutate(SiteDay = paste(Site,Day)) %>%
  mutate(SiteRep = paste(Site,Replicate))

#save object with average area per day and growth rate 7-14, 14-35, 35-56
AreaGrowthMeasures <- dplyr::select(AreabyExperiment, c(Day,Area,Site,Replicate)) %>%
  pivot_wider(names_from = Day, values_from = Area) %>%
  mutate(Growth714 = `Day 14`-`Day 07`) %>%
  mutate(Growth1435 = `Day 35`-`Day 14`) %>%
  mutate(Growth3556 = `Day 56`-`Day 35`) %>%
  mutate(Experiment = paste0(Site,Replicate))
#write_csv(AreaGrowthMeasures, "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/AverageAreaMeasures.csv")

AreabyExperiment$SiteDay <- factor(AreabyExperiment$SiteDay, levels =  c("CHOP Day 07","CN Day 07","UNC Day 07",
                                      "CHOP Day 14","CN Day 14", "UNC Day 14",
                                      "CHOP Day 35", "CN Day 35","UNC Day 35",  
                                      "CHOP Day 56", "CN Day 56", "UNC Day 56"))

pdf(paste(plot.here,"Jan10_Area_Overview_recolor",".pdf"), onefile = TRUE, width = 12, height = 8)
  ggplot(data=AreabyExperiment, aes(x=SiteDay, y=Area, group=SiteDay, fill = Site)) +
      #geom_line(aes(group = SiteRep, color="gray", size = 0.5)) +
      geom_line(aes(group = SiteRep,color = Site))+
    geom_violin()+
      #geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = TRUE) +
      geom_point(aes(shape = Replicate),size = 2)+
                 #,position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
      scale_fill_manual(values = Site.color) +
      scale_shape_manual(values=c(1,2,3,4,15,16,17))+
      ggtitle("Growth of hCOs")+
      theme_bw(base_size = 18) +
      ylab("Cross-Sectional Area (mm2)")+
      theme(axis.text.x = element_text(angle = 45))
  
dev.off()

# Statistics
#get test for all area measurements into one table.
#use replicates as biological sample
#correct for multiple comparisons 

D7ANOVA <- aov(Area~Site, data = filter(AreabyExperiment, Day == "Day 07"))
D7ANOVA <- summary(D7ANOVA)
test.area <- filter(AreabyExperiment, Day == "Day 07")
CHOPvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "UNC"))
CHOPvUNC <- t.test(Area ~ Site, data = filter(test.area, Site != "CN"))
UNCvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "CHOP"))
area.test.table <- as.data.frame(rbind(CHOPvCN,CHOPvUNC,UNCvCN))
area.test.table <- rownames_to_column(area.test.table, var = "TestName") %>%
  mutate(TestName = paste("D7", TestName))

D14ANOVA <- aov(Area~Site, data = filter(AreabyExperiment, Day == "Day 14"))
D14ANOVA <- summary(D14ANOVA)
test.area <- filter(AreabyExperiment, Day == "Day 14")
CHOPvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "UNC"))
CHOPvUNC <- t.test(Area ~ Site, data = filter(test.area, Site != "CN"))
UNCvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "CHOP"))
area.test.table.tmp <- as.data.frame(rbind(CHOPvCN,CHOPvUNC,UNCvCN))
area.test.table.tmp <- rownames_to_column(area.test.table.tmp, var = "TestName") %>%
  mutate(TestName = paste("D14", TestName))
area.test.table <-rbind(area.test.table.tmp, area.test.table)

D35ANOVA <- aov(Area~Site, data = filter(AreabyExperiment, Day == "Day 35"))
D35ANOVA <- summary(D35ANOVA)
test.area <- filter(AreabyExperiment, Day == "Day 35")
t.test(Area ~ Site, data = filter(test.area, Site != "UNC"))
t.test(Area ~ Site, data = filter(test.area, Site != "CN"))
t.test(Area ~ Site, data = filter(test.area, Site != "CHOP"))
CHOPvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "UNC"))
CHOPvUNC <- t.test(Area ~ Site, data = filter(test.area, Site != "CN"))
UNCvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "CHOP"))
area.test.table.tmp <- as.data.frame(rbind(CHOPvCN,CHOPvUNC,UNCvCN))
area.test.table.tmp <- rownames_to_column(area.test.table.tmp, var = "TestName") %>%
  mutate(TestName = paste("D35", TestName))
area.test.table <-rbind(area.test.table.tmp, area.test.table)

D56ANOVA <- aov(Area~Site, data = filter(AreabyExperiment, Day == "Day 56"))
D56ANOVA <- summary(D56ANOVA)
test.area <- filter(AreabyExperiment, Day == "Day 56")
t.test(Area ~ Site, data = filter(test.area, Site != "UNC"))
t.test(Area ~ Site, data = filter(test.area, Site != "CN"))
t.test(Area ~ Site, data = filter(test.area, Site != "CHOP"))
CHOPvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "UNC"))
CHOPvUNC <- t.test(Area ~ Site, data = filter(test.area, Site != "CN"))
UNCvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "CHOP"))
area.test.table.tmp <- as.data.frame(rbind(CHOPvCN,CHOPvUNC,UNCvCN))
area.test.table.tmp <- rownames_to_column(area.test.table.tmp, var = "TestName") %>%
  mutate(TestName = paste("D56", TestName))
area.test.table <-rbind(area.test.table.tmp, area.test.table)

AreaGrowthMeasures <- dplyr::select(AreabyExperiment, c(Day,Area,Site,Replicate)) %>%
  pivot_wider(names_from = Day, values_from = Area) %>%
  mutate(Growth714 = `Day 14`-`Day 07`) %>%
  mutate(Growth1435 = `Day 35`-`Day 14`) %>%
  mutate(Growth3556 = `Day 56`-`Day 35`) %>%
  mutate(Experiment = paste0(Site,Replicate))
#write_csv(AreaGrowthMeasures, "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/AverageAreaMeasures.csv")

D714ANOVA <- aov(Growth714~Site, data = AreaGrowthMeasures)
D714ANOVA <- summary(D714ANOVA)
t.test(Growth714 ~ Site, data = filter(AreaGrowthMeasures, Site != "UNC"))
t.test(Growth714 ~ Site, data = filter(AreaGrowthMeasures, Site != "CN"))
t.test(Growth714 ~ Site, data = filter(AreaGrowthMeasures, Site != "CHOP"))
CHOPvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "UNC"))
CHOPvUNC <- t.test(Area ~ Site, data = filter(test.area, Site != "CN"))
UNCvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "CHOP"))
area.test.table.tmp <- as.data.frame(rbind(CHOPvCN,CHOPvUNC,UNCvCN))
area.test.table.tmp <- rownames_to_column(area.test.table.tmp, var = "TestName") %>%
  mutate(TestName = paste("D714 growth", TestName))
area.test.table <-rbind(area.test.table.tmp, area.test.table)

D1435ANOVA <- aov(Growth1435~Site, data = AreaGrowthMeasures)
D1435ANOVA <- summary(D1435ANOVA)
t.test(Growth1435 ~ Site, data = filter(AreaGrowthMeasures, Site != "UNC"))
t.test(Growth1435 ~ Site, data = filter(AreaGrowthMeasures, Site != "CN"))
t.test(Growth1435 ~ Site, data = filter(AreaGrowthMeasures, Site != "CHOP"))
CHOPvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "UNC"))
CHOPvUNC <- t.test(Area ~ Site, data = filter(test.area, Site != "CN"))
UNCvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "CHOP"))
area.test.table.tmp <- as.data.frame(rbind(CHOPvCN,CHOPvUNC,UNCvCN))
area.test.table.tmp <- rownames_to_column(area.test.table.tmp, var = "TestName") %>%
  mutate(TestName = paste("D1435 growth", TestName))
area.test.table <-rbind(area.test.table.tmp, area.test.table)

D3556ANOVA <- aov(Growth3556~Site, data = AreaGrowthMeasures)
D3556ANOVA <- summary(D3556ANOVA)
t.test(Growth3556 ~ Site, data = filter(AreaGrowthMeasures, Site != "UNC"))
t.test(Growth3556 ~ Site, data = filter(AreaGrowthMeasures, Site != "CN"))
t.test(Growth3556 ~ Site, data = filter(AreaGrowthMeasures, Site != "CHOP"))
CHOPvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "UNC"))
CHOPvUNC <- t.test(Area ~ Site, data = filter(test.area, Site != "CN"))
UNCvCN <- t.test(Area ~ Site, data = filter(test.area, Site != "CHOP"))
area.test.table.tmp <- as.data.frame(rbind(CHOPvCN,CHOPvUNC,UNCvCN))
area.test.table.tmp <- rownames_to_column(area.test.table.tmp, var = "TestName") %>%
  mutate(TestName = paste("D3556 growth", TestName))
area.test.table <-rbind(area.test.table.tmp, area.test.table)

area.test.table <- mutate(area.test.table, p.adjust = p.adjust(p.value, method="BH")) %>%
  dplyr::select(TestName, p.value,p.adjust)

anova.area.test.table <- as_tibble(rbind(unlist(D7ANOVA)[9],unlist(D14ANOVA)[9],
      unlist(D35ANOVA)[9],unlist(D56ANOVA)[9],
      unlist(D714ANOVA)[9],unlist(D1435ANOVA)[9],unlist(D3556ANOVA)[9]))
anova.area.test.table$Test <- c("D7","D14","D35","D56","D714","D1435","D3556")
anova.area.test.table <- mutate(anova.area.test.table, p.adjust = p.adjust(`Pr(>F)1`, method="BH"))

#write_csv(area.test.table,"/work/users/r/o/roseg/IDDRC/IDDRCPlots/Area/MultipleComparisonCorrectPvalAreaMeasurements.csv")

#Cross-site differences in size
# PLOTS ###########################################################################
#Statistical Testing for Site Differences

ggplot(D7, aes(x=Site, y=AreaMM2,fill=Site)) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = TRUE) +
  geom_point(aes(shape=Replicate),size=0.5,
             position=position_jitterdodge(dodge.width=0.01), alpha=0.7)+
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=c(1,2,3,4,15,16,17))+
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw(base_size = 18) +
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(0.7,0.9,1.1)) + ggtitle("Day 7 Area")
ggsave("/work/users/r/o/roseg/IDDRC/IDDRCPlots/Area/CN_Registered_D7_Area.pdf",
       width = 6, height = 6, units = "in")

ggplot(D14, aes(x=Site, y=AreaMM2,fill=Site)) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = TRUE) +
  geom_point(aes(shape=Replicate),size=0.5,
             position=position_jitterdodge(dodge.width=0.01), alpha=0.7)+
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=c(1,2,3,4,15,16,17))+
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw(base_size = 18) +
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(2,3,4)) + ggtitle("Day 14 Area")
ggsave("/work/users/r/o/roseg/IDDRC/IDDRCPlots/Area/CN_Registered_D14_Area.pdf",
       width = 6, height = 6, units = "in")

ggplot(D35, aes(x=Site, y=AreaMM2, fill=Site)) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = TRUE) +
  geom_point(aes(shape=Replicate),size=0.5,
             position=position_jitterdodge(dodge.width=0.01), alpha=0.7)+
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=c(1,2,3,4,15,16,17))+
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw(base_size = 18) +
  ylab("Cross-Sectional Area (mm2)")+  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(5,5.5,6))+ggtitle("Day 35 Area")
ggsave("/work/users/r/o/roseg/IDDRC/IDDRCPlots/Area/CN_Registered_D35_Area.pdf",
       width = 6, height = 6, units = "in")

ggplot(D56, aes(x=Site, y=AreaMM2, fill=Site)) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = TRUE) +
  geom_point(aes(shape=Replicate),size=0.5,
             position=position_jitterdodge(dodge.width=0.01), alpha=0.7)+
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=c(1,2,3,4,15,16,17))+
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw(base_size = 18) +
  ylab("Cross-Sectional Area (mm2)")+
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(10,11,12))+ggtitle("Day 56 Area")
ggsave("/work/users/r/o/roseg/IDDRC/IDDRCPlots/Area/CN_Registered_D56_Area.pdf",
       width = 6, height = 6, units = "in")

ggplot(AreaGrowthMeasures, aes(x=Site, y=Growth714, fill=Site)) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = TRUE) +
  geom_point(aes(shape=Replicate),size=0.5,
             position=position_jitterdodge(dodge.width=0.01), alpha=0.7)+
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=c(1,2,3,4,15,16,17))+
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw(base_size = 18) +
  ggtitle("Day 7 to 14 Growth")+
  ylab("Cross-Sectional Area Growth (mm2)")+
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(0.6,0.7,0.8))
ggsave("/work/users/r/o/roseg/IDDRC/IDDRCPlots/Area/CN_Registered_D714Growth_Area.pdf",
       width = 6, height = 6, units = "in")

ggplot(AreaGrowthMeasures, aes(x=Site, y=Growth1435, fill=Site)) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = TRUE) +
  geom_point(aes(shape=Replicate),size=0.5,
             position=position_jitterdodge(dodge.width=0.01), alpha=0.7)+
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=c(1,2,3,4,15,16,17))+
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw(base_size = 18) +
  ggtitle("Day 14 to 35 Growth")+
  ylab("Cross-Sectional Area Growth (mm2)")+
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(1.6,1.7,1.8))
ggsave("/work/users/r/o/roseg/IDDRC/IDDRCPlots/Area/CN_Registered_D1435Growth_Area.pdf",
       width = 6, height = 6, units = "in")

ggplot(AreaGrowthMeasures, aes(x=Site, y=Growth3556, fill=Site)) +
  geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = TRUE) +
  geom_point(aes(shape=Replicate),size=0.5,
             position=position_jitterdodge(dodge.width=0.01), alpha=0.7)+
  scale_fill_manual(values = Site.color) +
  scale_shape_manual(values=c(1,2,3,4,15,16,17))+
  theme(axis.text.x = element_text(angle = 45)) +
  theme_bw(base_size = 18) +
  ggtitle("Day 35 to 56 Growth")+
  ylab("Cross-Sectional Area Growth (mm2)")+
  geom_signif(test = "t.test",comparisons = list(c("CN", "CHOP"),
                                                 c("CN", "UNC"),
                                                 c("UNC","CHOP")),
              y_position = c(5,5.2,5.4))
ggsave("/work/users/r/o/roseg/IDDRC/IDDRCPlots/Area/CN_Registered_D3556Growth_Area.pdf",
       width = 6, height = 6, units = "in")
