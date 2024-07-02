library(Seurat)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(here)
library(reshape2)
library(readxl)
library(speckle)
library(limma)
library(plotrix)

setwd("/work/users/r/o/roseg/")
HERE <- "/work/users/r/o/roseg/single-cell_reproducibility/"
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"

# OUTPUT FILES #########################################################################################################

dir.pdf <- paste0(HERE,"doc/seurat/pdf/")

# INPUT FILES ##########################################################################################################

# integrated and clustered dataset, rescaled seurat object
seurat.named.rds <- paste0(HERE,"results/May12023_IDDRC_PGP1_UNCR3_reference_integrated_clustered_seurat_object.rds")

# add names and order to propeller
name.convert <- read_excel(paste0(db.here,"CellTypeProportions/Propeller.Cluster.Name.Conversion.xlsx"))

#primary scRNAseq from human fetal tissue Polioudakis
gw1718 <- read_csv(paste0(db.here,"CellTypeProportions/PrimaryTissueforProportionComparisons.csv"))

# GLOBALS ##############################################################################################################
mixing_pallete <-c('#5F95B2', '#BCC6E5','#B09AB1','#F1C4DC', '#EA9F8B','#F8C893',
                   '#89A48C','#369F48','#DB7B87','#E12228','#B177B3', '#2179B4',
                   '#F47B20','#F89B40','#F15A29')

# Load PGP1 Integrated Data ############
seur.pgp1 <- readRDS(seurat.named.rds)

# use "Site" as group for propeller
seur.pgp1@meta.data$group <- seur.pgp1@meta.data$Site

df.cellData <- as_tibble(seur.pgp1@meta.data)

# Day 14 ##########################################################################################
df.cellData %>%
    filter(Site != "Zelda") %>%
    filter(Day == "D14") -> df.cellData.d14

IDDRCD14Prop <- propeller(clusters = df.cellData.d14$integrated_snn_res.0.6, sample = df.cellData.d14$Sample, group = df.cellData.d14$Site)
IDDRCD14Prop$clusters <- row.names(IDDRCD14Prop)
# write_csv(IDDRCD14Prop, "/work/users/r/o/roseg/single-cell_reproducibility/results/IDDRCD14Prop_BroadClusters.csv")

plotCellTypeProps(clusters = df.cellData.d14$integrated_snn_res.0.6, sample = df.cellData.d14$Sample) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values= mixing_pallete)

# build dataframe of sample proportions for plotting
props.d14 <- getTransformedProps(clusters = df.cellData.d14$integrated_snn_res.0.6, sample = df.cellData.d14$Sample)

#write.csv(props.d14$Proportions, paste0(HERE,"results/IDDRC_D14_Proportions_all_samples_BroadClusters.csv"))

df.props.d14 <- tibble(melt(props.d14$Proportions, value.name = "proportion"))
df.props.d14 %>%
    left_join(tibble(melt(props.d14$Counts, value.name = "counts")), by = c("clusters", "sample"))
df.props.d14 %>%
    left_join(tibble(melt(props.d14$TransformedProps, value.name = "transformed_props")), by = c("clusters", "sample"))

#make your own stacked bar plots
#df.props.d14 <- inner_join(df.props.d84, name.convert, by ="clusters")
#df.props.d14 <- mutate(df.props.d84, CTO = paste(CellTypeOrder, CellType))

df.props.d14 %>%
    ggplot(aes(fill=CTO, y=proportion, x=sample)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values= mixing_pallete) +
    geom_bar(position="stack", stat="identity")

#calculate Coefficient of Variation
p.cluster.CV <-  df.props.d14 %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(cluster.cv = sd(proportion) / mean(proportion) )

#write_csv(p.cluster.CV, paste0(db.here,"CellTypeProportions/IDDRC_D14_AllClusterCoefficentofVariation.csv"))


# Day 84 ################################################################################################
df.cellData %>%
    filter(!Site == "Zelda") %>%
    filter(Day == "D84") -> df.cellData.d84

IDDRCD84Prop <- propeller(clusters = df.cellData.d84$integrated_snn_res.0.6, sample = df.cellData.d84$Sample, group = df.cellData.d84$Site)
IDDRCD84Prop$clusters <- row.names(IDDRCD84Prop)
#write_csv(IDDRCD84Prop, "/work/users/r/o/roseg/single-cell_reproducibility/results/IDDRCD84Prop_BroadClusters.csv")

plotCellTypeProps(clusters = df.cellData.d84$integrated_snn_res.0.6, sample = df.cellData.d84$Sample) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values= mixing_pallete)

# build dataframe of sample proportions for plotting
props.d84 <- getTransformedProps(clusters = df.cellData.d84$integrated_snn_res.0.6, sample = df.cellData.d84$Sample)

#write.csv(props.d84$Proportions, paste0(HERE,"results/IDDRC_D84_Proportions_all_samples_broadsamples.csv"))

df.props.d84 <- tibble(melt(props.d84$Proportions, value.name = "proportion"))
df.props.d84 %>%
    left_join(tibble(melt(props.d84$Counts, value.name = "counts")), by = c("clusters", "sample"))
df.props.d84 %>%
    left_join(tibble(melt(props.d84$TransformedProps, value.name = "transformed_props")), by = c("clusters", "sample"))

#df.props.d84 <- inner_join(df.props.d84, name.convert, by ="clusters")
#df.props.d84 <- mutate(df.props.d84, CTO = paste(CellTypeOrder, CellType))

#make your own stacked bar plots
df.props.d84 %>%
    ggplot(aes(fill=CTO, y=proportion, x=sample)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values= mixing_pallete) +
    geom_bar(position="stack", stat="identity")

#calculate Coefficient of Variation
p.cluster.CV <-  df.props.d84 %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(cluster.cv = sd(proportion) / mean(proportion))

#write_csv(p.cluster.CV, paste0(db.here,"CellTypeProportions/IDDRC_D84_AllClusterCoefficentofVariation.csv"))

# for ZELDA/D56 ###########################################################################################
df.cellData %>%
    filter(Site == "Zelda") -> df.cellData.d14

#only 3 samples (N of 1 per group) so can't perform any ANOVA
# build dataframe of sample proportions for plotting
props.d14 <- getTransformedProps(clusters = df.cellData.d14$integrated_snn_res.0.6, sample = df.cellData.d14$Sample)

#write.csv(props.d14$Proportions, paste0(HERE,"results/ZELDA_D56_Proportions_all_samples.csv"))

df.props.d14 <- tibble(melt(props.d14$Proportions, value.name = "proportion"))
df.props.d14 %<>%
    left_join(tibble(melt(props.d14$Counts, value.name = "counts")), by = c("clusters", "sample"))
df.props.d14 %<>%
    left_join(tibble(melt(props.d14$TransformedProps, value.name = "transformed_props")), by = c("clusters", "sample"))

df.props.d14 <- inner_join(df.props.d14, name.convert, by ="clusters")
df.props.d14 <- mutate(df.props.d14, CTO = paste(CellTypeOrder, CellType))

#make your own stacked bar plots
df.props.d14 %>%
    ggplot(aes(fill=CTO, y=proportion, x=sample)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values= mixing_pallete) +
    geom_bar(position="stack", stat="identity")

#calculate Coefficient of Variation
p.cluster.CV <-  df.props.d14 %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(cluster.cv = sd(proportion) / mean(proportion))

#write_csv(p.cluster.CV, paste0(db.here,"CellTypeProportions/ZeldaAllClusterCoefficentofVariation.csv"))

# For comparing Day 56 to Day 84 ################################################################################################
df.cellData %>%
  filter(Day != "D14") %>%
  filter(Site != "CHOP") %>%
  filter(Site != "CN")-> df.cellData.d5684

IDDRCD84Prop <- propeller(clusters = df.cellData.d5684$integrated_snn_res.0.6, sample = df.cellData.d5684$Sample, group = df.cellData.d5684$Site)
IDDRCD84Prop$clusters <- row.names(IDDRCD84Prop)
#write_csv(IDDRCD84Prop,"/work/users/r/o/roseg/IDDRC/IDDRCDatabase/CellTypeProportions/ProportionChangesfromD56toD84.csv")

#for Primary Tissue Comparison ###########################################################
mixing_pallete <-c("#82463D","#DB6767","#C47C6E","#B177B3","#5C2D7F",
                   "#336699","#383838","#7D7D7D","#E87862","#C05285",
                   "#4D55A5","#5BBD7C","#A55528","#EFB6CD","#8C9AED","#2B8540")
#make label with Donor, GW, library
#add counts for each cluster & total cell count for each donor
gw1718 <- mutate(gw1718, label = paste("Donor",Donor, "GW",Gestation_week, Library, sep = " "))
gw1718counts <- gw1718 %>% count(label, Cluster, sort = TRUE)
gw1718counts <- dplyr::rename(gw1718counts, CountCluster = n)
gw1718countsTotal <- gw1718counts %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(TotalCell = sum(CountCluster))
gw1718counts <- inner_join(gw1718counts, gw1718countsTotal, by = "label")
gw1718counts <- mutate(gw1718counts, PercentCluster = CountCluster/TotalCell)


GW1718Plot <- gw1718counts %>%
    ggplot(aes(fill=Cluster, y=PercentCluster, x=label)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    scale_fill_manual(values= mixing_pallete) +
    geom_bar(position="stack", stat="identity")+
    xlab("") + ylab("Percent Cluster")
GW1718Plot

gw1718prop <- inner_join(gw1718counts,gw1718)
gw1718prop <- distinct(dplyr::select(gw1718prop, c("label","Cluster","PercentCluster","Donor")))

gw1718.cluster.SEM <-  gw1718prop %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(cluster.sem = std.error(PercentCluster))

#calculate Coefficient of Variation
primary.cluster.CV <-  gw1718prop %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(cluster.cv = sd(PercentCluster) / mean(PercentCluster) * 100)

#write_csv(gw1718.cluster.SEM, paste0(db.here,"PrimaryAllClusterClusterCoefficentofVariation.csv"))


#Cell Type Proportions for Bhaduri #################################################################
seurat.kreig.rds <- "/proj/steinlab/projects/IVIV_scRNA/youngsook_pine/extData_annotation/kreigstein_primary_seuratObject_SCTv2/Kreigstein_primary_integrated.rds"
seur.kreig <- readRDS(seurat.kreig.rds)

df.primary <- seur.kreig@meta.data

#each age looks pretty different
plotCellTypeProps(clusters = df.primary$clustInfo, sample = df.primary$sampleInfo) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

props.d14 <- getTransformedProps(clusters = df.primary$clustInfo, sample = df.primary$sampleInfo)

#write.csv(props.d14$Proportions, paste0(HERE,"results/Primary_Sample_Proportions_Bhaduri.csv"))

df.props.d14 <- tibble(melt(props.d14$Proportions, value.name = "proportion"))
df.props.d14 %<>%
    left_join(tibble(melt(props.d14$Counts, value.name = "counts")), by = c("clusters", "sample"))
df.props.d14 %<>%
    left_join(tibble(melt(props.d14$TransformedProps, value.name = "transformed_props")), by = c("clusters", "sample"))

# Kreigstein (Bhaduri) cell types

reduce.cluster <- read_csv(paste0(db.here, "PrimaryTissueReducedCluster.csv"))
reduce.cluster <- dplyr::rename(reduce.cluster, clusters = Original)
df.props.d14 <- full_join(reduce.cluster, df.props.d14)
df.props.d14 <- mutate(df.props.d14, OriginalName = paste(Class,Type,Subtype, sep="."))

#make your own stacked bar plots
mixing_pallete.p <-c("#51261A","#2791CC","#E72071","#000000","#D76327","#AAAAAA","#53AF59")
df.props.d14 %>%
    ggplot(aes(fill=BroadCluster, y=proportion, x=sample)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values= mixing_pallete.p)

#calculate S.E.M for all clusters and broad clusters

p.cluster.SEM <-  df.props.d14 %>%
    dplyr::group_by(OriginalName) %>%
    dplyr::summarise(cluster.sem = std.error(proportion))

p.broadcluster.SEM <-  df.props.d14 %>%
    dplyr::group_by(BroadCluster) %>%
    dplyr::summarise(cluster.sem = std.error(proportion))

#calculate Coefficient of Variation
primary.cluster.CV <-  gw1718prop %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(cluster.cv = sd(PercentCluster) / mean(PercentCluster) * 100)

write_csv(primary.cluster.CV, paste0(db.here,"BhaduriClusterClusterCoefficentofVariation.csv"))

# GRUFFI ########################################################################################
df.cellData <- read.csv("/work/users/r/o/roseg/IDDRC/IDDRCDatabase/CellTypeProportions/GRUFFI_iddrc_metadata_june21.csv")
df.cellData %>%
    dplyr::filter(Site != "Zelda") %>%
    dplyr::filter(Day == "D84")%>%
    dplyr::filter(integrated_snn_res.0.6 == "9") -> df.cellData.d84

#mutate is.stressed to be a character
df.cellData.d84$is.Stressed <- as.character(df.cellData.d84$is.Stressed)

#Error in tmixture.matrix(out$t, stdev.unscaled, df.total, proportion,  :
#Dims of tstat and stdev.unscaled don't match

IDDRCD84Prop <- propeller(clusters = df.cellData.d84$is.Stressed, sample = df.cellData.d84$Sample, group = df.cellData.d84$Site)
IDDRCD84Prop$clusters <- row.names(IDDRCD84Prop)

#fist plot of cells
plotCellTypeProps(clusters = df.cellData.d84$is.Stressed, sample = df.cellData.d84$Sample) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values= mixing_pallete)


# ICC ############################################################################################
###########################################
### Calculating intraclass-correlation coefficient for D14/D84
library(psych)
library(qpcR)

#Pull just the information to test for intra class correlation coefficient,
# where experiment is a "rater" of celltype proportion for each Site
#must do one proportion at a time
df.props.d84 <- mutate(df.props.d84, Site = str_extract(df.props.d84$sample, ".*_"))
props.d84 <- filter(df.props.d84, CellType == "Upper layer neuron I Layer 2/3/4")
D84ICC <- unique(dplyr:::select(props.d84, c(proportion, Site,sample)))
#make raters vector
raters <- c("rater1", "rater2", "rater3", "rater4", "rater5")
#make a vector of all Site with rank information
Site <- as.vector(unique(D84ICC$Site))
#create an object to slot for loop information into.
D84ICC.output <- tibble()

#for every Site create a new column with rater label for every sample
for (i in seq_along(Site)){
    nm <- Site[i]
    SampletoRater <- dplyr::select(props.d84, c(sample,Site))
    SampletoRater <- unique(SampletoRater)
    SampletoRater <- dplyr::filter(SampletoRater, Site == nm)
    Sitefilter <- as_tibble(qpcR:::cbind.na(SampletoRater, raters))
    D84ICC.output <- bind_rows(D84ICC.output, Sitefilter)
}

IDDRC.iccready <- dplyr::distinct(inner_join(D84ICC.output, D84ICC))
#filter R object to just DCCID, raters and RankD14
IDDRC.iccready <- dplyr:::select(IDDRC.iccready, c(Site, raters, proportion))
#make Robject have "raters" as variables and Donors as rows, then remove DCCID label
IDDRC.icc <- IDDRC.iccready %>%
    pivot_wider(names_from = raters, values_from = proportion, values_fn = {mean})
IDDRC.icc <-dplyr:::select(IDDRC.icc, -c(Site))

#perform intraclass correlation coefficent calculation
#average_random_raters because were are interested in the mean of all EBID ranks and the ranks come from unique experiments
result.ICC <- ICC(IDDRC.icc)
result.ICC

detach(package:qpcR,unload=TRUE)
detach(package:psych,unload=TRUE)
detach(package:lme4,unload=TRUE)
detach(package:MASS,unload=TRUE)
