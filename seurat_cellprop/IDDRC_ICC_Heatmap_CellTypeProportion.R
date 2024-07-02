#try with IDDRC ICC matrix as heatmap
#using cell type proportion
# cell types as rows, Experiment as columns
library(tidyverse)
library(irr)
library(ComplexHeatmap)
library(circlize)

setwd("/work/users/r/o/roseg/")
HERE <- "/work/users/r/o/roseg/"
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
rds.here <- "/work/users/r/o/roseg/IDDRC/IDDRC_EVOS/"

df.props.d84 <- read_csv(paste0(db.here,"CellTypeProportions/D84ProportionsbyExperiment.csv"))
df.props.d84 <- df.props.d84[-1,]

df.props.d14.r <- read_csv(paste0(db.here,"CellTypeProportions/D14ProportionsbyExperiment.csv"))
df.props.d14.r <- df.props.d14.r[-1,]

df.prop.d56.r <- read_csv(paste0(db.here,"CellTypeProportions/ZELDA_D56_Proportions_all_samples.csv"))
colnames(df.prop.d56.r) <-  c("CellType","UNC_Day56_Rep1","UNC_Day56_Rep2","UNC_Day56_Rep3")
df.prop.d56.r <- t(df.prop.d56.r)  
colnames(df.prop.d56.r) <- df.prop.d56.r[1,]
df.prop.d56.r <- df.prop.d56.r[-1,]
Experiment <- rownames(df.prop.d56.r)
df.prop.d56.r <- cbind(df.prop.d56.r,Experiment) %>%
  as_data_frame()

df.props.d14 <- read_csv(paste0(HERE,"results/IDDRC_D14_Proportions_all_samples_BroadClusters.csv"))
df.props.d14 <- t(df.props.d14)
df.props.d14 <- df.props.d14[-1,]
Experiment <- rownames(df.props.d14)
df.props.d14 <- cbind(as.data.frame(df.props.d14),Experiment)

df.props.d84 <- read_csv(paste0(HERE,"results/IDDRC_D84_Proportions_all_samples_broadsamples.csv"))
df.props.d14 <- t(df.props.d14)
df.props.d14 <- df.props.d14[-1,]
Experiment <- rownames(df.props.d14)
df.props.d14 <- cbind(as.data.frame(df.props.d14),Experiment)


dt <- df.prop.d56.r

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
#iddrc.d84.icc <- output.dt
#iddrc.d14.icc <- output.dt

Site <- sub('_.*','',rownames(output.dt))
#TotalCell <- dt$TotalCell
#day_col <-  setNames(c(t(my.colors[1:3])),c('Day14','Day56','Day84'))
#change to site colors in paper
site_col <-setNames(c("#12783D", "#882155", "#332F85"),c("CHOP","CN","UNC"))

#set colors for interpreting ICC
col_icc <- colorRamp2(c(0.4, 0.5, 1), c("#2179B4",  "white","#DDCB77"))

#pdf(paste0('Seurat/plots/ICC.heatmap.res',res,'.scale.data.pdf'),useDingbats=F)
p<- Heatmap(as.matrix(output.dt),name='ICC',col = col_icc,
            cluster_rows=F,cluster_columns=F,
            row_names_gp=gpar(fontsize=8),column_names_gp=gpar(fontsize=8),right_annotation = rowAnnotation(Site=Site,col=list(Site=site_col)),
            top_annotation = HeatmapAnnotation(Site=Site,col=list(Site=site_col))) + ggtitle("IDDRC Day 84")
print(p)
#dev.off()

#< 0.5 poor
# Everything above moderate for 
#0.5 to 07.75 moderate
#0.75 to 0.9 is good
#0.9+ is excellent