#% export LD_LIBRARY_PATH=/nas/longleaf/rhel8/apps/gcc/11.2.0/lib64:/nas/longleaf/rhel8/apps/gcc/6.3.0/lib:/nas/longleaf/rhel8/apps/gcc/6.3.0/lib64:/nas/longleaf/rhel8/apps/r/4.1.3/lib64/R/lib
options(rgl.useNULL = TRUE) # without this option, ngl fails under my environment
#options(stringsAsFactors=F)
# set library path
#.libPaths(c(.libPaths(),'/work/users/n/a/nanam/tools/R/4.1/','/work/users/n/a/nanam/tools/R/4.2/','~/R/x86_64-pc-linux-gnu-library/4.1/'))
# load all nessesary tools
#devtools::load_all("/work/users/n/a/nanam/tools/GRUFFI/gruffi/")
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

########################
### load SeuratObject
########################
inputfn <- '/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/20230321_PGP1_UNCR3_reference_integrated_SCTRescale_seurat_object.rds'
if(file.exists(inputfn)){
 S.obj <- readRDS(inputfn)
}
 ########################
 ### redo UMAP projection
 ########################
 message('making 3d UMAP projection')

 # 1. backup old umap
 S.obj@misc$reductions.backup$umap2d <- S.obj@reductions$umap

 # 2. calculate and backup the new 3D umap | "dimensions=3"
 S.obj <- Seurat.utils::SetupReductionsNtoKdimensions(obj = S.obj, nPCs = 9, dimensions=3, reduction="umap")

 # 3. Recall the old 2D umap from @misc
 S.obj@reductions$umap <- S.obj@misc$reductions.backup$umap2d

 #Not sure what this is from Nana but I think Mike's object has the tsne info alrady

 #message('plotting')
 #not sure why Nana is doing this? Seurat.utils is not plottin
# Seurat.utils::clUMAP("integrated_snn_res.0.6", obj = S.obj, save.plot = T,reduction="tsne",prefix=paste0(OutputDir,'00'),cols=my.colors)
 S.obj <- SetIdent(S.obj, value = S.obj@meta.data$integrated_snn_res.0.6)

 S.obj <- RenameIdents(S.obj, `0` = "Neuroepithelial stem cells",
                       `6` = "Dividing neural progenitor cells, S",
                       `10` = "Dividing neural progenitor cells, G2",
                       `12` = "Medial Pallium/Marginal Zone",
                       `14` = "Dividing intermediate progenitors, S",
                       `4` = "Intermediate progenitors",
                       `11` = "Radial glia",
                       `5` = "Outer radial glia",
                       `2` = "Lower layer neuron a",
                       `3` = "Lower layer neuron b",
                       `13` = "Upper layer neuron a, Layer 2/3/4",
                       `1` = "Upper layer neuron b, layer 2/3",
                       `7` = "Pan cortical neuron",
                       `8` = "Pan neuron/Cajal - Retzius",
                       `9` = "Unspecified neuron")

DimPlot(S.obj, reduction = "tsne", raster = FALSE, cols = mixing_pallete, label=FALSE)

 #message('saving data')
 #saveRDS(S.obj, file='Seurat/post-processed/Gruffi/20230109.keep/pgp1_day14_n_day84/9_9/0.Res0.8.3D.rds')


########################
### load GO-terms and gene-sets
########################
message('loding GO data')
load(file="/work/users/r/o/roseg/single-cell_reproducibility/data/Gruffi.stress.marker.rdata")

# Gruffi works best if you partition cells into groups of 100-200 cells. Find the corresponding clustering resolution by:
S.obj <- aut.res.clustering(obj = S.obj)
granule.res.4.gruffi <- S.obj@misc$gruffi$'optimal.granule.res'

message(granule.res.4.gruffi)

S.obj.o <- reassign.small.clusters(S.obj, ident = granule.res.4.gruffi) # will be stored in meta data column as "seurat_clusters.reassigned"
#Above, granules with <30 cells are cell-by-cell re-assigned to a neighboring granule (by default based on Euclidean distance between the mean of cell groups in 3dim UMAP space).
#The reassigned granules are suffixed as :

granule.res.4.gruffi <- paste0(granule.res.4.gruffi, '.reassigned')
message(granule.res.4.gruffi)

########################
### 3. Pathway scoring
#######################

# Glycolytic process	GO:0006096
S.obj <- GO_score_evaluation(obj = S.obj, GO_term = go1, save.UMAP = FALSE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
n_go1 <- S.obj@misc$GO[[sub(':','\\.',go1)]] %>% length()
p1<- FeaturePlot(S.obj,reduction="tsne",features=paste0("Score.",sub(':','\\.',go1)),min.cutoff="q05",max.cutoff="q95")+
ggplot2::labs(title=paste0("Score.",go1, " Glycolytic process"),caption = paste('Score calc. from ',n_go1, "expr. genes from BioMart.", paste0("https://www.ebi.ac.uk/QuickGO/search/",go1)))

# ER stress 	GO:0034976
S.obj <- GO_score_evaluation(obj = S.obj, GO_term = go2, save.UMAP = FALSE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
n_go2 <- S.obj@misc$GO[[sub(':','\\.',go2)]] %>% length()
p2<- FeaturePlot(S.obj,reduction="tsne",features=paste0("Score.",sub(':','\\.',go2)),min.cutoff="q05",max.cutoff="q95")+
ggplot2::labs(title=paste0("Score.",go2, " ER stress"),caption = paste('Score calc. from ',n_go2, "expr. genes from BioMart.", paste0("https://www.ebi.ac.uk/QuickGO/search/",go2)))

# Gliogenesis		GO:0042063
S.obj <- GO_score_evaluation(obj = S.obj, GO_term = go3, save.UMAP = FALSE, new_GO_term_computation = T, clustering = granule.res.4.gruffi, plot.each.gene = F)
n_go3 <- S.obj@misc$GO[[sub(':','\\.',go3)]] %>% length()
p3<- FeaturePlot(S.obj,reduction="tsne",features=paste0("Score.",sub(':','\\.',go3)),min.cutoff="q05",max.cutoff="q95")+
ggplot2::labs(title=paste0("Score.",go3, " Gliogenesis"),caption = paste('Score calc. from ',n_go3, "expr. genes from BioMart.", paste0("https://www.ebi.ac.uk/QuickGO/search/",go3)))

#These functions store the resulting scores in S.obj@meta.data.
GO.info <- S.obj@meta.data[grepl('GO',names(S.obj@meta.data))]
saveRDS(GO.info,file="/work/users/r/o/roseg/single-cell_reproducibility/data/gruffiJun15.1.GO.metadata.rds")

library(patchwork)
png('Seurat/post-processed/Gruffi/20230109.keep/pgp1_day14_n_day84/9_9/3.GO.Score.TSNE.png',height=1000,width=700)
p1+p2+p3
dev.off()

########################
### 4. Stress filtering
########################
# Create score names:
(i1 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go1))
(i2 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go2))
(i3 <- Stringendo::kppu(granule.res.4.gruffi, 'cl.av', go3))

save(S.obj,i1,i2,i3,file="/work/users/r/o/roseg/single-cell_reproducibility/data/gruffiJun15.2.S.obj.w.PathwayScored.rdata")

# need to run through desktop terminal to open Shiny.app
if(FALSE){
S.obj <- Shiny.GO.thresh(obj = S.obj,
                                stress.ident1 = i1,
                                stress.ident2 = i2,
                                notstress.ident3 = i3,
                                plot.cluster.shiny = "orig.ident")

 saveRDS(S.obj, file="/work/users/r/o/roseg/single-cell_reproducibility/data/gruffiJun15.5.S.obj.w.Stressed.anno.rds")
 stressed.anno <- S.obj@meta.data$is.Stressed
 saveRDS(stressed.anno,file="Seurat/post-processed/Gruffi/20230109.keep/pgp1_day14_n_day84/9_9/6.Stressed.anno.rds")
}

## stats
if(FALSE){
 library(dplyr)
 library(ggplot2)
 library(ggsci)
 #stressed.anno <- readRDS('Seurat/post-processed/Gruffi/20230109.keep/pgp1_day14_n_day84/9_9/06.Stressed.anno.rds')
 #S.obj$is.Stressed <- stressed.anno
 S.obj <- readRDS('Seurat/post-processed/Gruffi/20230109.keep/pgp1_day14_n_day84/9_9/05.S.obj.w.Stressed.anno.rds')
 Seurat.utils::clUMAP(obj=S.obj,'is.Stressed', label =F,reduction="tsne",splitby='integrated_snn_res.0.6',nr.cols=4,h=8,cols=c('#6F99AD','#BC3C29FF',"black","red","blue","green"),prefix=paste0(OutputDir,'Stressed.split'))
 Seurat.utils::clUMAP(obj=S.obj,'is.Stressed', label =F,reduction="tsne",cols=c('#6F99AD','#BC3C29FF'),prefix=paste0(OutputDir,'Stressed'))

 cluster.anno <- read.delim('conf/Cell.Annotation.by.Rose.20230510.tsv',sep="\t")
 Stressedcell.stats <- S.obj@meta.data %>% group_by(integrated_snn_res.0.6) %>%
 mutate(cluster=integrated_snn_res.0.6,totalCellC = n()) %>%
 group_by(cluster,totalCellC,is.Stressed) %>%
 summarize(UnStressed=n()) %>% filter(is.Stressed==F) %>%
 select(cluster,totalCellC,UnStressed) %>%
 mutate(Stressed = totalCellC - UnStressed) %>%
 mutate(PcStressed = paste0(signif((Stressed/totalCellC)*100,digits=3),'%'))

 mt <- Stressedcell.stats %>% ungroup() %>% select(cluster,UnStressed,Stressed) %>% reshape2::melt()
 pdf('Seurat/post-processed/Gruffi/20230109.keep/pgp1_day14_n_day84/9_9/07.Stressed.cell.counts.pdf',width=9)
 ggplot(mt)+geom_bar(aes(x=cluster,y=value,fill=variable),stat="identity") +
 geom_text(data=Stressedcell.stats,aes(x=cluster,y=Stressed+200,label=PcStressed),hjust=0.5,vjust=0.5,size=2.5,color="white")+
 scale_x_discrete(labels=paste0(cluster.anno$cluster,' (',cluster.anno$ClusterName,')'))+
 scale_fill_manual(name="Cell Identification",values=c('#6F99AD','#BC3C29FF')) + labs(y="# of cells")+
 theme_bw()+theme(legend.pos="bottom",axis.text.x=element_text(angle=60,hjust=1))
 dev.off()

 library(ggplot2)

 Stressedcell.stats.perSample <- S.obj@meta.data %>% group_by(sampleNames) %>%
 mutate(totalCellC = n()) %>%
 group_by(sampleNames,is.Stressed,totalCellC) %>%
 summarize(UnStressed=n()) %>% filter(is.Stressed==F) %>%
 ungroup() %>%
 select(sampleNames,totalCellC,UnStressed) %>%
 mutate(Stressed = totalCellC - UnStressed) %>%
 mutate(PcStressed = paste0(signif((Stressed/totalCellC)*100,digits=3),'%'))

 mt <- Stressedcell.stats.perSample %>% ungroup() %>% select(sampleNames,UnStressed,Stressed) %>% reshape2::melt()
 mt$donorID <- sub('DCCID_(\\d+)_.*','\\1',mt$sampleNames)
 mt$Day <- sub('.*Day','Day',mt$sampleNames)
 mt$sampleNames <- sub('_Day.*','',sub('DCCID_','',mt$sampleNames))
 Stressedcell.stats.perSample$donorID <- sub('DCCID_(\\d+)_.*','\\1',Stressedcell.stats.perSample$sampleNames)
 Stressedcell.stats.perSample$Day <- sub('.*Day','Day',Stressedcell.stats.perSample$sampleNames)
 Stressedcell.stats.perSample$sampleNames <- sub('_Day.*','',sub('DCCID_','',Stressedcell.stats.perSample$sampleNames))
 p1<- mt %>% filter(Day == "Day14") %>%ggplot() +
 geom_bar(aes(sampleNames,y=value,fill=variable),stat="identity") +
 geom_text(data=(Stressedcell.stats.perSample %>% filter(Day=="Day14")) ,aes(x=sampleNames,y=Stressed+200,label=PcStressed),hjust=0.5,vjust=0.5,size=2.5,color="black",angle=90)+
 scale_fill_manual(name="Cell Identification",values=c('#6F99AD','#BC3C29FF')) + labs(y="# of cells")+
 theme_bw()+theme(legend.pos="bottom",axis.text.x=element_text(angle=60,hjust=1)) +
 facet_grid(.~donorID,scales="free",space="free")
 p2<- mt %>% filter(Day == "Day84") %>%ggplot() +
 geom_bar(aes(sampleNames,y=value,fill=variable),stat="identity") +
 geom_text(data=(Stressedcell.stats.perSample %>% filter(Day=="Day84")) ,aes(x=sampleNames,y=Stressed+200,label=PcStressed),hjust=0.5,vjust=0.5,size=2.5,color="black",angle=90)+
 scale_fill_manual(name="Cell Identification",values=c('#6F99AD','#BC3C29FF')) + labs(y="# of cells")+
 theme_bw()+theme(legend.pos="bottom",axis.text.x=element_text(angle=60,hjust=1)) +
 facet_grid(.~donorID,scales="free",space="free")
 p3<- Stressedcell.stats.perSample %>% ggplot(aes(x=100*(Stressed/totalCellC),fill=Day))+
 geom_histogram(binwidth=0.5,color="black")+
 labs(subtitle="proportion of stressed cells")+
 scale_fill_manual(values=c('Day14'='navy','Day84'='orange'))+
 theme_bw()
 pdf('Seurat/post-processed/Gruffi/20230109.keep/pgp1_day14_n_day84/9_9/08.Stressed.cell.counts.perSample.pdf',width=13)
 p1
 p2
 p3+plot_spacer()
 dev.off()
}
