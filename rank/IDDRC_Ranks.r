library(tidyverse)
library(readxl)
library(plotrix)
#Set you working directory to be where your files are
setwd("/work/users/r/o/roseg/EVOSAnalysis/IDDRC_EVOS/Ranking")

#load in the ranking data
Rankings <- read_csv("IDDRCTidyRanking_March30.csv")
#Add a unique Experiment ID 
Rankings <- Rankings %>% mutate(Rankings,
                                Experiment = paste(Site,Replicate))
#making sure you can splits the ranking data by day
D14Rank <- filter(Rankings, Day =="D14")
D35Rank <- filter(Rankings, Day =="D35")
D56Rank <- filter(Rankings, Day=="D56")

DaystoTest <- c('D14','D35','D56')

#check for count of hCO in each experiment, if normally distributed, and if biased by Ranker
for (T in 1:length(DaystoTest)){
  PD <- filter(Rankings, Day == DaystoTest[T])
 
byExperiment <- ggplot(PD, aes(x=Experiment, fill=Ranker))+
  theme(axis.text.x = element_text(angle = 45))+
  geom_bar()

histo <- ggplot(PD, aes(x=Rank))+
    geom_histogram(stat = "count")

byRanker <- (ggplot(PD, aes(x=Ranker, y=Rank, fill=Ranker))) +
  #geom_boxplot(outlier.shape = 1, lwd=0.2) + 
  geom_violin()+
  geom_point(size = 1, position=position_jitterdodge(), alpha=0.9) +
  ylim(0,8) +  theme(axis.text.x = element_text(angle = 45)) 

byExperiments <- ggplot(PD, aes(x=Experiment, y=Rank, color=Ranker))+
  geom_violin()+
  #geom_boxplot(outlier.shape = 1, lwd=0.2) + 
  geom_point(size = 1, position=position_jitterdodge(), alpha=0.9) +
  ylim(0,8) +  theme(axis.text.x = element_text(angle = 45)) 

byAssignment <- (ggplot(PD, aes(x=Assignment, y=Rank, color=Ranker))) +
  #geom_boxplot(outlier.shape = 1, lwd=0.2) + 
  geom_violin()+
  geom_point(size = 1, position=position_jitterdodge(), alpha=0.9) +
  ylim(0,8) +  theme(axis.text.x = element_text(angle = 45)) 

bySite <- (ggplot(PD, aes(x=Site, y=Rank, color=Ranker))) +
  #geom_boxplot(outlier.shape = 1, lwd=0.2) + 
  geom_violin()+
  geom_point(size = 1, position=position_jitterdodge(), alpha=0.9) +
  ylim(0,8) +  theme(axis.text.x = element_text(angle = 45)) 

##save graphs output to a PDF, changing the working directory to a place where I want the graphs
setwd("/work/users/r/o/roseg/EVOSAnalysis/IDDRC_EVOS/Ranking/ForEmma")
pdf(paste(DaystoTest[T],".pdf"), onefile = TRUE, width = 12, height = 8)
  print(byExperiment)
  print(histo)
  print(byExperiments)
  print(byRanker)
  print(byAssignment)
  print(bySite)
dev.off()

}

###plots averaging rank by replicate
IDConvert <-  distinct(select(Rankings, c(Site,Replicate,Experiment)))

#you'll need to change this value for each timepoint (D14, D35, D56)
MeanRank <- D14Rank
T <- "D14Rank"

#summarize the data by Experiment
MRD <- MeanRank %>%
  dplyr::group_by(Experiment) %>%
  dplyr::summarise(RankMean = mean(Rank, na.rm = TRUE),
                   Rank.sem = std.error(Rank, na.rm = TRUE))
#add the Site and Replicate information back
MRD <- inner_join(MRD, IDConvert)

#Graph average rank of each replicate by Site 
bySite <- (ggplot(MRD, aes(x=Site, y=RankMean, fill=Site))) +
  geom_boxplot(lwd=0.2) +  ylim(0,8) +
  scale_shape_manual(values=c(1,2,3,4,15,16,17))+
  geom_point(aes(shape = Replicate),size = 1, position=position_jitterdodge(),alpha=0.9) +
  theme(axis.text.x = element_text(angle = 45)) 
bySite

#save graph as .pdf
pdf(paste(T,"AverageReplicatebySite",".pdf"), onefile = TRUE, width = 12, height = 8)
bySite
dev.off()

####Test for Site differences in Rank, you can view in R but we're also saving it as a .csv
SiteANOVA <- aov(RankMean~Site, data = MRD)
SiteANOVA <- summary(SiteANOVA)
capture.output(SiteANOVA,file =  paste(T,"SiteAnova.csv"))

#do some post-hoc testing to tell which site is different from which site
CHOP <- filter(MRD, Site =="CHOP")
UNC <- filter(MRD, Site =="UNC")
CN <- filter(MRD, Site=="CN")
t.test(CHOP$RankMean,UNC$RankMean)
capture.output(SiteANOVA, file = paste(T,"CHOPvUNCTTest.csv"))
t.test(CN$RankMean,UNC$RankMean)
capture.output(SiteANOVA,file =  paste(T,"CNvUNCTTest.csv"))
t.test(CN$RankMean,CHOP$RankMean)
capture.output(SiteANOVA,file =  paste(T,"CNvCHOPTTest.csv"))


###########################################
### Calculating intraclass-correlation coefficient for D14
library(psych)
library(qpcR)
#Get a list of all possible Site, Replicates
IDconvert <-  distinct(dplyr::select(Rankings, c(Site,Replicate,Experiment)))
#Pull just the information to test for intra class correlation coefficient, Where EBID is a "rater" of rank for each DCCID
D14ICC <- unique(dplyr:::select(D14Rank, c(Rank, Site, Replicate)))
#make raters vector
raters <- c("rater1", "rater2", "rater3", "rater4", "rater5", "rater6, rater7")
#make a vector of all DCCIDs with rank information
Site <- as.vector(unique(D14Rank$Site))
#create an object to slot for loop information into. 
D14ICC.output <- tibble()

#for every DCCID create a new column with rater label for every EBID
for (i in seq_along(Site)){
  nm <- Site[i]
  ReplicatetoRater <- filter(IDconvert, Site == nm)
  IDDRC.DCCIDfilter <- as_tibble(qpcR:::cbind.na(ReplicatetoRater, raters))
  D14ICC.output <- bind_rows(D14ICC.output, IDDRC.DCCIDfilter)
}

#Add Average Rank associated with each EBID
IVIV.iccready <- distinct(inner_join(D14ICC.output, D14ICC))
#filter R object to just DCCID, raters and RankD14
IVIV.iccready <- dplyr:::select(IVIV.iccready, c(Site, raters, Rank))
#make Robject have "raters" as variables and Donors as rows, then remove DCCID label 
IVIV.icc <- IVIV.iccready %>%
  pivot_wider(names_from = raters, values_from = Rank, values_fn = {mean})
IVIV.icc <-dplyr:::select(IVIV.icc, -c(Site))

#perform intraclass correlation coefficent calculation
#average_random_raters because were are interested in the mean of all EBID ranks and the ranks come from unique experiments
result.ICC <- ICC(IVIV.icc)
result.ICC
#for D14 this is ICC = 0.33

detach(package:qpcR,unload=TRUE)
detach(package:psych,unload=TRUE)
detach(package:lme4,unload=TRUE)
detach(package:MASS,unload=TRUE)
