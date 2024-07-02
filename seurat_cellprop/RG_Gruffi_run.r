library('pdist')
library('raster')
library('rgl')
library("Stringendo")
library("CodeAndRoll2")
library("ReadWriter")
library("MarkdownHelpers")
library("Markdownreports")
library("ggExpress")
library("Seurat.utils")
library("gruffi")

my.colors<-c('#5F95B2','#BCC6E5','#B09AB1','#F1C4DC','#EA9F8B','#F8C893','#89A48C','#369F48',
             '#DB7B87', '#E12228','#B177B3','#2179B4','#F47B20','#F89B40','#F15A29')

sessionInfo()

ww.assign_to_global <- function(name, value,  pos = 1, verbose = TRUE){
  if (verbose) iprint(name, "defined as:", value) # , "is a new global environment variable"
  assign(name, value, envir=as.environment(pos) )
}

OutputDir <- '/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/Gruffi/20230109.keep/pgp1_day14_n_day84/9_9/'
MarkdownReports::create_set_OutDir(setDir = F, getwd(),'/',OutputDir)

inputfn <- '/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/20230321_PGP1_UNCR3_reference_integrated_SCTRescale_seurat_object.rds'
if(file.exists(inputfn)){
  S.obj <- readRDS(inputfn)
}
S.obj <- SetIdent(S.obj, value = S.obj@meta.data$integrated_snn_res.0.6)

#check plotting
DimPlot(S.obj, reduction = "tsne", raster = FALSE, cols = mixing_pallete, label=FALSE)
FeaturePlot(S.obj, features = 'Score.GO.0042063')

### load GO-terms and gene-sets
########################
message('loding GO data')
load(file="/work/users/r/o/roseg/single-cell_reproducibility/data/Gruffi.stress.marker.rdata")

# Gruffi works best if you partition cells into groups of 100-200 cells. Find the corresponding clustering resolution by:
S.obj <- aut.res.clustering(obj = S.obj)
granule.res.4.gruffi <- S.obj@misc$gruffi$'optimal.granule.res'

message(granule.res.4.gruffi)

S.obj <- reassign.small.clusters(S.obj, ident = granule.res.4.gruffi) # will be stored in meta data column as "seurat_clusters.reassigned"
#Above, granules with <30 cells are cell-by-cell re-assigned to a neighboring granule (by default based on Euclidean distance between the mean of cell groups in 3dim UMAP space).
#The reassigned granules are suffixed as :

granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')
message(granule.res.4.gruffi)

########################
### 3. Pathway scoring
#######################

# Glycolytic process	GO:0006096
S.obj <- GO_score_evaluation(obj = S.obj, GO_term = go1, save.UMAP = TRUE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
n_go1 <- S.obj@misc$GO[[sub(':','\\.',go1)]] %>% length()

# ER stress 	GO:0034976
S.obj <- GO_score_evaluation(obj = S.obj, GO_term = go2, save.UMAP = FALSE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
n_go2 <- S.obj@misc$GO[[sub(':','\\.',go2)]] %>% length()

# Gliogenesis		GO:0042063
S.obj <- GO_score_evaluation(obj = S.obj, GO_term = go3, save.UMAP = FALSE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
n_go3 <- S.obj@misc$GO[[sub(':','\\.',go3)]] %>% length()

#These functions store the resulting scores in S.obj@meta.data.
GO.info <- S.obj@meta.data[grepl('GO',names(S.obj@meta.data))]
#saveRDS(GO.info,file="/work/users/r/o/roseg/single-cell_reproducibility/data/gruffiJun15.1.GO.metadata.rds")

### 4. Stress filtering
########################
# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go3))

#Threshold proposal: GO:0006096: 5.65, GO:0034976: 3.43, GO:0042063: 6.33
GO96.5per <- quantile(GO.info$Score.GO.0006096, probs = seq(0, 1, 1/20))
GO76.5per <- quantile(GO.info$Score.GO.0034976, probs = seq(0, 1, 1/20))
GO63.5per <- quantile(GO.info$Score.GO.0042063, probs = seq(0, 1, 1/20))

S.obj <- Shiny.GO.thresh(obj = S.obj,
                         stress.ident1 = i1,
                         stress.ident2 = i2,
                         notstress.ident3 = i3,
                         plot.cluster.shiny = "orig.ident")

#used default thresholds. Need to save in Shiny
stressed.anno <- S.obj@meta.data$is.Stressed
GruffiMetaData <- S.obj@meta.data
write.csv(GruffiMetaData, "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/CellTypeProportions/GRUFFI_iddrc_metadata_june21.csv")

cluster.anno <- read_csv('/work/users/r/o/roseg/IDDRC/IDDRCscRNA/ClustertoCellType.csv')
cluster.anno <- rename(cluster.anno, cluster = Cluster)
cluster.anno$cluster <- as.character(cluster.anno$cluster)
Stressedcell.stats <- S.obj@meta.data %>% group_by(integrated_snn_res.0.6) %>%
mutate(cluster=integrated_snn_res.0.6,totalCellC = n()) %>%
  group_by(cluster,totalCellC,is.Stressed) %>%
  summarize(UnStressed=n()) %>% filter(is.Stressed==F) %>%
  dplyr::select(cluster,totalCellC,UnStressed) %>%
  mutate(Stressed = totalCellC - UnStressed) %>%
  mutate(PcStressed = paste0(signif((Stressed/totalCellC)*100,digits=3),'%'))

mt <- Stressedcell.stats %>% ungroup() %>% dplyr::select(cluster,UnStressed,Stressed) %>% reshape2::melt()
#mt <- inner_join(mt, cluster.anno, by = "cluster")
ggplot(mt)+geom_bar(aes(x=cluster,y=value,fill=variable),stat="identity") +
  geom_text(data=Stressedcell.stats,aes(x=cluster,y=Stressed+200,label=PcStressed),hjust=0.5,vjust=0.5,size=2.5,color="white")+
  scale_x_discrete(labels=paste0(cluster.anno$cluster,' (',cluster.anno$CellType,')'))+
  scale_fill_manual(name="Cell Identification",values=c('#6F99AD','#BC3C29FF')) + labs(y="# of cells")+
  theme_bw()+theme(legend.pos="bottom",axis.text.x=element_text(angle=60,hjust=1))

#brokendown by sample
Stressedcell.stats.perSample <- S.obj@meta.data %>% group_by(Sample) %>%
  mutate(totalCellC = n()) %>%
  group_by(Sample,is.Stressed,totalCellC) %>%
  summarize(UnStressed=n()) %>% filter(is.Stressed==F) %>%
  ungroup() %>%
  dplyr::select(Sample,totalCellC,UnStressed) %>%
  mutate(Stressed = totalCellC - UnStressed) %>%
  mutate(PcStressed = paste0(signif((Stressed/totalCellC)*100,digits=3),'%'))

mt <- Stressedcell.stats.perSample %>% ungroup() %>% dplyr::select(Sample,UnStressed,Stressed) %>% reshape2::melt()
mt$Site <- gsub("_.*","",mt$Sample)
mt$Day <- sub('.*D','Day',mt$Sample)
mt$Day <- gsub("_.*","",mt$Day)
mt$Rep <- sub('.*R','Rep',mt$Sample)
#objects for plotting with readable names?
#Stressedcell.stats.perSample$donorID <- sub('DCCID_(\\d+)_.*','\\1',Stressedcell.stats.perSample$sampleNames)
Stressedcell.stats.perSample$Day <- sub('.*Day','Day',Stressedcell.stats.perSample$Sample)
#Stressedcell.stats.perSample$sampleNames <- sub('_Day.*','',sub('DCCID_','',Stressedcell.stats.perSample$sampleNames))

p1<- mt %>% filter(Day == "Day14") %>%ggplot() +
  geom_bar(aes(Sample,y=value,fill=variable),stat="identity") +
  geom_text(data=(Stressedcell.stats.perSample %>% filter(Day=="Day14")),
            aes(x=Sample,y=Stressed+200,label=PcStressed),hjust=0.5,vjust=0.5,size=2.5,color="black",angle=90)+
  scale_fill_manual(name="Cell Identification",values=c('#6F99AD','#BC3C29FF')) + labs(y="# of cells")+
  theme_bw()+theme(legend.pos="bottom",axis.text.x=element_text(angle=60,hjust=1)) +
  facet_grid(.~Site,scales="free",space="free")

p2<- mt %>% filter(Day == "Day84") %>%ggplot() +
  geom_bar(aes(Sample,y=value,fill=variable),stat="identity") +
  geom_text(data=(Stressedcell.stats.perSample %>% filter(Day=="Day84")) ,aes(x=Sample,y=Stressed+200,label=PcStressed),hjust=0.5,vjust=0.5,size=2.5,color="black",angle=90)+
  scale_fill_manual(name="Cell Identification",values=c('#6F99AD','#BC3C29FF')) + labs(y="# of cells")+
  theme_bw()+theme(legend.pos="bottom",axis.text.x=element_text(angle=60,hjust=1)) +
  facet_grid(.~Site,scales="free",space="free")

p3<- Stressedcell.stats.perSample %>% ggplot(aes(x=100*(Stressed/totalCellC),fill=Day))+
  geom_histogram(binwidth=0.5,color="black")+
  labs(subtitle="proportion of stressed cells")+
  scale_fill_manual(values=c('Day14'='navy','Day84'='orange'))+
  theme_bw()

#brokendown by sample and cluster #make new colummn in metadata with sample_day_cluster
S.obj@meta.data <- mutate(S.obj@meta.data, ExpDayCluster = paste(Sample,integrated_snn_res.0.6, sep= "!"))
Stressedcell.stats.perSample.cl <- S.obj@meta.data %>% group_by(ExpDayCluster) %>%
    mutate(totalCellC = n()) %>%
    group_by(ExpDayCluster,is.Stressed,totalCellC) %>%
    summarize(UnStressed=n()) %>% filter(is.Stressed==F) %>%
    ungroup() %>%
    dplyr::select(ExpDayCluster,totalCellC,UnStressed) %>%
    mutate(Stressed = totalCellC - UnStressed) %>%
    mutate(PcStressed = paste0(signif((Stressed/totalCellC)*100,digits=3),'%'))

write_csv(Stressedcell.stats.perSample.cl, file = "/work/users/r/o/roseg/single-cell_reproducibility/data/June16.Gruffi.Suggested.Threshold.IDDRC.Stress.Experiment.Day.Cluster.csv")

df.cluster.site <- read_csv("/work/users/r/o/roseg/single-cell_reproducibility/data/June16.Gruffi.Suggested.Threshold.IDDRC.Stress.Experiment.Day.Cluster.csv")
df.cluster.site <- mutate(df.cluster.site, Site = gsub("_.*","",df.cluster.site$ExpDayCluster))
df.cluster.site <- mutate(df.cluster.site, Cluster = gsub(".*!","",df.cluster.site$ExpDayCluster))
df.cluster.site <- mutate(df.cluster.site, Day = gsub(".*_D", "\\1", df.cluster.site$ExpDayCluster))
df.cluster.site <- mutate(df.cluster.site, Day = substr(df.cluster.site$Day, 1, 2) )
df.cluster.site <- mutate(df.cluster.site, Rep = gsub("!.*","",df.cluster.site$ExpDayCluster))
df.cluster.site <- mutate(df.cluster.site, Rep = gsub(".*_","",df.cluster.site$Rep))

# remove Zelda
df.cluster.site <- filter(df.cluster.site, Site != "Zelda")
df.cluster.site<- mutate(df.cluster.site, PerStress = Stressed/totalCellC)
check.stress.only <- filter(df.cluster.site, Stressed >0)
cl.9 <- filter(df.cluster.site, Cluster == 9)
cl.11 <- filter(df.cluster.site, Cluster == 11)

CL9D84 <- cl.9 %>% filter(Day == "84") %>%
    ggplot(aes(fill=Site, y=PerStress, x=ExpDayCluster)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values= Site.color)+
    ylab("Percent Stressed Cells") + ylim(0,0.2) +
    ggtitle("Day 84 Cluster 9")

CL9D14 <- cl.9 %>% filter(Day == "14") %>%
    ggplot(aes(fill=Site, y=PerStress, x=ExpDayCluster)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values= Site.color)+
    ylab("Percent Stressed Cells") + ylim(0,0.2) +
    ggtitle("Day 14 Cluster 9")

CL11D84 <- cl.11 %>% filter(Day == "84") %>%
    ggplot(aes(fill=Site, y=PerStress, x=ExpDayCluster)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values= Site.color)+
    ylab("Percent Stressed Cells") + ylim(0,0.2) +
    ggtitle("Day 84 Cluster 11")

CL11D14 <- cl.11 %>% filter(Day == "14") %>%
    ggplot(aes(fill=Site, y=PerStress, x=ExpDayCluster)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values= Site.color)+
    ylab("Percent Stressed Cells") + ylim(0,0.2) +
    ggtitle("Day 14 Cluster 11")

SiteANOVA <- aov(PerStress~Site, data = filter(cl.11, Day == "84"))
SiteANOVA <- summary(SiteANOVA)

SiteANOVA <- aov(PerStress~Site, data = filter(cl.11, Day == "14"))
SiteANOVA <- summary(SiteANOVA)

SiteANOVA <- aov(PerStress~Site, data = filter(cl.9, Day == "84"))
SiteANOVA <- summary(SiteANOVA)

SiteANOVA <- aov(PerStress~Site, data = filter(cl.9, Day == "14"))
SiteANOVA <- summary(SiteANOVA)



#Feature Plotting removed for now
p1<- FeaturePlot(S.obj,reduction="tsne",
                 features=paste0("Score.",sub(':','\\.',go1)),min.cutoff="q05",max.cutoff="q95")+
  ggplot2::labs(title=paste0("Score.",go1,
                             " Glycolytic process"),
                caption = paste('Score calc. from ',n_go1,
                                "expr. genes from BioMart.", paste0("https://www.ebi.ac.uk/QuickGO/search/",go1)))


p1<- FeaturePlot(S.obj,reduction="tsne",features=paste0("Score.",sub(':','\\.',go1)),min.cutoff="q05",max.cutoff="q95")+
  ggplot2::labs(title=paste0("Score.",go1, " Glycolytic process"),caption = paste('Score calc. from ',n_go1, "expr. genes from BioMart.", paste0("https://www.ebi.ac.uk/QuickGO/search/",go1)))

p2<- FeaturePlot(S.obj,reduction="tsne",features=paste0("Score.",sub(':','\\.',go2)),min.cutoff="q05",max.cutoff="q95")+
  ggplot2::labs(title=paste0("Score.",go2, " ER stress"),caption = paste('Score calc. from ',n_go2, "expr. genes from BioMart.", paste0("https://www.ebi.ac.uk/QuickGO/search/",go2)))

p3<- FeaturePlot(S.obj,reduction="tsne",features=paste0("Score.",sub(':','\\.',go3)),min.cutoff="q05",max.cutoff="q95")+
  ggplot2::labs(title=paste0("Score.",go3, " Gliogenesis"),caption = paste('Score calc. from ',n_go3, "expr. genes from BioMart.", paste0("https://www.ebi.ac.uk/QuickGO/search/",go3)))

FeaturePlot(S.obj, reduction = "tsne", features = "integrated_snn_res.70.reassigned_cl.av_GO:0034976")
