#module load r/4.2.2
#inspired by https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/pseudobulk_DESeq2_scrnaseq.md and Nana Matoba
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(muscat))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))

run_n_plot_PCA<- function(cluster_metadata,cluster_counts,plot.label,outliers=NULL,plotDay=F,return.counts=F){
 # Check matching of matrix columns and metadata rows
 print(all(colnames(cluster_counts) == rownames(cluster_metadata)))
 # Create DESeq2 object
 dds <- DESeqDataSetFromMatrix(cluster_counts,
                               colData = cluster_metadata,
                               design = ~ as.factor(Site) +(1|sample_id))

 #design = ~ as.factor(Site)

 # Perform size factor normalization
 dds <- estimateSizeFactors(dds)
 # remove genes with no (low) expression (the following three lines is just for plot purpose to see how many donors can be survived
 keeped.genes <- do.call(rbind,lapply(c(0:10),
  function(x){keep = rowSums(counts(dds, normalized=TRUE) >= x ) >= ceiling(0.1*ncol(dds));
  data.frame(thr=x,geneN=sum(keep))}))
 keep <- rowSums(counts(dds, normalized=TRUE) >= 1 ) >= ceiling(0.1*ncol(dds))
 dds <- dds[keep,]
 message('set minimun sampleN as ', ceiling(0.1*ncol(dds)))
 vsd <- vst(dds)
 log2counts = assay(vsd)
 pca_logcpm <- prcomp(t(log2counts))
 #pca_points <- as.data.frame(pca_logcpm$x[,1:20])
 #not always calculating 20 PCs. set to lower number of PCs?
 pca_points <- as.data.frame(pca_logcpm$x[,1:5])
 pca_points$cluster_sample_id <- rownames(pca_points)
 pca_points <- pca_points %>% left_join(cluster_metadata)
 return(list(expr=log2counts,pca=pca_points))
}

#count_ls is list of cluster names?
counts_ls <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
#for(res in c(0.2,0.4,0.6,0.8,1)){
 load(paste0("/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/aggr.sum.cluster_sample.res",res,".rdata"))
 ### separate by Day
 res=0.6
  i = 1
  output <- vector(mode="list", 2*length(counts_ls))
  for(idxlist in names(counts_ls)){
   for(Day in c('D14','D84')){
    message(res,' ', idxlist, ' ', Day,' ',i)
    idx <- which(names(counts_ls) == idxlist)
    cluster_counts <- counts_ls[[idx]]
    cluster_metadata <- metadata_ls[[idx]]
    cluster_counts.tmp <- cluster_counts[,grep(Day,colnames(cluster_counts))]
    cluster_metadata.tmp <- cluster_metadata[grep(Day,rownames(cluster_metadata)),]
    message(nrow(cluster_metadata.tmp))
    if(nrow(cluster_metadata.tmp) == 0) {
     message('no cells are detected. skip.')
    }
    #removed because already have minimum of 3 reps per site
    #else if (length(unique(cluster_metadata.tmp$group_id)) < 10){
    # message(length(unique(cluster_metadata.tmp$group_id)))
    # message('less thana 10 donors are detected. skip.')
    #}
    else{
     message(length(unique(cluster_metadata.tmp$Site)))
     pca.res = run_n_plot_PCA(cluster_metadata.tmp,cluster_counts.tmp,paste0('cluster:',idxlist,' (',Day,')'),NULL,plotDay=F,return.counts=T)
     output[[i]]$gexpr = pca.res$expr
     output[[i]]$pca = pca.res$pca
    }
    i = i+1
   }
  }
 saveRDS(output,file=paste0('PseudoBulk/postGruffi/VST.Expr.PerSample_Day_Cluster.rm.lowExpGenes.',res,'.rds'))
#}

 #lots of DEG id'd by DESeq2 per cluster. Seem real? possibly set higher thresholds for cells?
