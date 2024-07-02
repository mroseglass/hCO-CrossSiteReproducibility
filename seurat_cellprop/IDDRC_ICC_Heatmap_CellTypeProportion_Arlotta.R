library(Seurat)
library(tidyverse)
library(readr)
library(magrittr)
library(ggplot2)
library(here)
library(reshape2)
library(readxl)
library(speckle)
library(plotrix)
library(irr)
library(circlize)
library(ComplexHeatmap)

setwd("/work/users/r/o/roseg/")
HERE <- "/work/users/r/o/roseg/single-cell_reproducibility/"
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
rds.here <- "/work/users/r/o/roseg/IDDRC/ArlottaGrouphCOData/"

# OUTPUT FILES #########################################################################################################
#Or try Velasco data
df.cellData.all <- readRDS(paste0(rds.here,"3mo_harm_111120.rds"))
df.cellData.wt <- as_tibble(df.cellData.all@meta.data)

ArlottaProp <- propeller(clusters = df.cellData.wt$FinalName, sample = df.cellData.wt$org, group = df.cellData.wt$orig.ident)
ArlottaProp$clusters <- row.names(ArlottaProp)
#check if similar to paper
plotCellTypeProps(clusters = df.cellData.wt$FinalName, sample = df.cellData.wt$org) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values= mixing_pallete)

props.d14 <- getTransformedProps(clusters = df.cellData.wt$FinalName, sample = df.cellData.wt$org)
df.props.d14 <- tibble(melt(props.d14$Proportions, value.name = "proportion"))
df.props.d14 %>%
  left_join(tibble(melt(props.d14$Counts, value.name = "counts")), by = c("clusters", "sample"))
df.props.d14 %>%
  left_join(tibble(melt(props.d14$TransformedProps, value.name = "transformed_props")), by = c("clusters", "sample"))

dt <- as_tibble(props.d14$Proportions)
dt <- dt %>% pivot_wider(names_from = clusters, values_from = n)
dt <- dplyr::rename(dt, Experiment = sample)
rownames(dt) <- dt$Experiment
#dt <-  dt %>% dplyr::select(all_of(cluster.list))

output<-matrix(NA,nrow=length(dt$Experiment),ncol=length(dt$Experiment))
colnames(output) <- dt$Experiment
rownames(output) <- dt$Experiment
nc <- ncol(dt)-1

for(j in 1:nrow(output)){
  for(k in 1:nrow(output)){
    if(j==k){
      output[j,k] <-NA
    }else{
      tmp <- data.frame(t(dt[dt$Experiment %in% c(colnames(output)[j],rownames(output)[k]),c(1:nc)]))
      colnames(tmp) <- tmp[1,]
      tmp <- tmp[-1,]
      tmp[,1] <- as.numeric(tmp[,1])
      tmp[,2] <- as.numeric(tmp[,2])
      icc.res<-icc(tmp,model='twoway', type='agreement',unit='single')$value #twoway
      output[j,k] = icc.res
    }
  }
}
output.dt<-data.frame(output)

Dset <- sub('[0-9]','',rownames(output.dt))

#set colors for interpreting ICC
col_icc <- colorRamp2(c(0.4, 0.5, 1), c("#2179B4",  "white","#DDCB77"))

#pdf(paste0('Seurat/plots/ICC.heatmap.res',res,'.scale.data.pdf'),useDingbats=F)
p<- Heatmap(as.matrix(output.dt),name='ICC',col = col_icc,
            cluster_rows=F,cluster_columns=F,
            row_names_gp=gpar(fontsize=8),column_names_gp=gpar(fontsize=8),
            right_annotation = rowAnnotation(Dataset=Dset),
            top_annotation = HeatmapAnnotation(Dataset=Dset)) +
  ggtitle("Arlotta ICC")
print(p)

library(correlation)
velasco.icc <- output.dt
velasco.icc.z <- as.vector(as.matrix(z_fisher(r = velasco.icc)))
iddrc.d84.icc <- output.dt

#within site ICC for D84
a <- iddrc.d84.icc[1:5,1:5]
b <- iddrc.d84.icc[6:9,6:9]
c <- iddrc.d84.icc[10:13,10:13]
iddrc.d84.icc.z.within <- c(as.vector(as.matrix(z_fisher(r = a))),as.vector(as.matrix(z_fisher(r = b))),as.vector(as.matrix(z_fisher(r = c))))

#across site ICC for D84
a <- iddrc.d84.icc[1:5,6:13]
b <- iddrc.d84.icc[6:9,1:5]
c <- iddrc.d84.icc[6:9,10:13]
d <- iddrc.d84.icc[10:13,1:9]

iddrc.d84.icc.z.across <- c(as.vector(as.matrix(z_fisher(r = a))),as.vector(as.matrix(z_fisher(r = b))),
                            as.vector(as.matrix(z_fisher(r = c))),as.vector(as.matrix(z_fisher(r = d))))

iddrc.d14.icc <- output.dt
#within site ICC for D14
a <- iddrc.d14.icc[1:5,1:5]
b <- iddrc.d14.icc[6:9,6:9]
c <- iddrc.d14.icc[10:12,10:12]
iddrc.d14.icc.z.within <- c(as.vector(as.matrix(z_fisher(r = a))),as.vector(as.matrix(z_fisher(r = b))),as.vector(as.matrix(z_fisher(r = c))))

#across site ICC for D84
a <- iddrc.d14.icc[1:5,6:12]
b <- iddrc.d14.icc[6:9,1:5]
c <- iddrc.d14.icc[6:9,10:12]
d <- iddrc.d14.icc[10:12,1:9]

iddrc.d14.icc.z.across <- c(as.vector(as.matrix(z_fisher(r = a))),as.vector(as.matrix(z_fisher(r = b))),
                            as.vector(as.matrix(z_fisher(r = c))),as.vector(as.matrix(z_fisher(r = d))))

#remove duplicate values from velasco.icc
velasco.icc.z <- unique(velasco.icc.z)
max_length <- max(length(velasco.icc.z), length(iddrc.d14.icc.z.within), length(iddrc.d14.icc.z.across),
                  length(iddrc.d84.icc.z.within),length(iddrc.d84.icc.z.across))
length(velasco.icc.z) <- max_length                      
length(iddrc.d84.icc.z.within) <- max_length
length(iddrc.d84.icc.z.across) <- max_length
length(iddrc.d14.icc.z.within) <- max_length 
length(iddrc.d14.icc.z.across) <- max_length 

icc.val <- rbind(velasco.icc.z, iddrc.d84.icc.z.across,iddrc.d84.icc.z.within,iddrc.d14.icc.z.across,iddrc.d14.icc.z.within)
icc.val[,1]<-rownames(icc.val)
icc.val <- as.data.frame(icc.val)
icc.val <- dplyr::rename(icc.val, "dataset" = "V1")
icc.val <- pivot_longer(icc.val, cols = !dataset, values_to = "icc")
icc.val$icc <- as.numeric(icc.val$icc)

library(ggsignif)
#Violin plot for icc
ggplot(icc.val, aes(x=dataset, y=icc, fill=dataset)) +
  geom_violin() + 
  geom_point(size=0.2,position = position_jitterdodge()) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  ggtitle("ICC Z-scores") + 
  ylab("ICC Z Score") +
  geom_signif(comparisons = list(c("iddrc.d14.icc.z.within", "iddrc.d14.icc.z.across")), y_position= 2.4) +
  geom_signif(comparisons = list(c("iddrc.d84.icc.z.within", "iddrc.d84.icc.z.across")), y_position= 2.4) +
  geom_signif(comparisons = list(c("iddrc.d14.icc.z.within", "iddrc.d84.icc.z.within")), y_position= 3) +
  geom_signif(comparisons = list(c("iddrc.d14.icc.z.across", "iddrc.d84.icc.z.across")), y_position= 2.7) +
  geom_signif(comparisons = list(c("iddrc.d14.icc.z.within", "velasco.icc.z")), y_position= 4.2) + 
  geom_signif(comparisons = list(c("iddrc.d14.icc.z.across", "velasco.icc.z")), y_position= 3.8) + 
  geom_signif(comparisons = list(c("iddrc.d84.icc.z.within", "velasco.icc.z")), y_position= 4.8) + 
  geom_signif(comparisons = list(c("iddrc.d84.icc.z.across", "velasco.icc.z")), y_position= 4.5) + 
  theme_bw() + ylim(0,5.5)


#####################
#primary scRNAseq from human fetal tissue Polioudakis
gw1718 <- read_csv(paste0(db.here,"CellTypeProportions/PrimaryTissueforProportionComparisons.csv"))

gw1718 <- mutate(gw1718, label = paste("Donor",Donor, "GW",Gestation_week, Library, sep = " "))
gw1718counts <- gw1718 %>% count(label, Cluster, sort = TRUE)
gw1718counts <- dplyr::rename(gw1718counts, CountCluster = n)
gw1718countsTotal <- gw1718counts %>%
  dplyr::group_by(label) %>%
  dplyr::summarise(TotalCell = sum(CountCluster))
gw1718counts <- inner_join(gw1718counts, gw1718countsTotal, by = "label")
gw1718counts <- mutate(gw1718counts, PercentCluster = CountCluster/TotalCell)

dt <- dplyr::select(gw1718counts,c(label,Cluster,PercentCluster))
dt <- dt %>% pivot_wider(names_from = Cluster, values_from = PercentCluster)
dt <- dplyr::rename(dt, Experiment = label)
rownames(dt) <- dt$Experiment
#dt <-  dt %>% dplyr::select(all_of(cluster.list))

output<-matrix(NA,nrow=length(dt$Experiment),ncol=length(dt$Experiment))
colnames(output) <- dt$Experiment
rownames(output) <- dt$Experiment
nc <- ncol(dt)-1

for(j in 1:nrow(output)){
  for(k in 1:nrow(output)){
    if(j==k){
      output[j,k] <-NA
    }else{
      tmp <- data.frame(t(dt[dt$Experiment %in% c(colnames(output)[j],rownames(output)[k]),c(1:nc)]))
      colnames(tmp) <- tmp[1,]
      tmp <- tmp[-1,]
      tmp[,1] <- as.numeric(tmp[,1])
      tmp[,2] <- as.numeric(tmp[,2])
      icc.res<-icc(tmp,model='twoway', type='agreement',unit='single')$value
      output[j,k] = icc.res
    }
  }
}
output.dt<-data.frame(output)

Dset <- rownames(output.dt)

#set colors for interpreting ICC
col_icc <- colorRamp2(c(0.4, 0.7, 1), c("#2179B4",  "white","#DDCB77"))

#pdf(paste0('Seurat/plots/ICC.heatmap.res',res,'.scale.data.pdf'),useDingbats=F)
p<- Heatmap(as.matrix(output.dt),name='ICC',col = col_icc,
            cluster_rows=F,cluster_columns=F,
            row_names_gp=gpar(fontsize=8),column_names_gp=gpar(fontsize=8),
            right_annotation = rowAnnotation(Dataset=Dset),
            top_annotation = HeatmapAnnotation(Dataset=Dset)) +
  ggtitle("Primary ICC")
print(p)



