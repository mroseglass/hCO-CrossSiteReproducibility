library(tidyverse)
library(DESeq2)
#module load r/4.2.2
#inspired by https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/pseudobulk_DESeq2_scrnaseq.md and Nana Matoba
library(data.table)
library(Seurat)
library(SingleCellExperiment)
library(muscat)
library(stringr)
library(Matrix)
library(dplyr)
library(DESeq2)

results.here <- "/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/DESeq2/LTRtest_reduced1_july1823/"
plot.here <- "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/DESeq2/LTR_reducedTest/"
counts_ls <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
res=0.6
load(paste0("/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/aggr.sum.cluster_sample.res",res,".rdata"))

#for loop for cell type and cluster
output = list()
i = 1
idxlist = "0"
Day = "D14"
for(idxlist in names(counts_ls)){
  for(Day in c('D14','D84')){
    idx <- which(names(counts_ls) == idxlist)
    cluster_counts <- counts_ls[[idx]]
    cluster_metadata <- metadata_ls[[idx]]
    cluster_counts.tmp <- cluster_counts[,grep(Day,colnames(cluster_counts))]
    cluster_metadata.tmp <- cluster_metadata[grep(Day,rownames(cluster_metadata)),]
    message(nrow(cluster_metadata.tmp))

     message(length(unique(cluster_metadata.tmp$Site)))

     cluster_metadata.tmp.shrink <- dplyr::select(cluster_metadata.tmp, c(Site))
     print(all(colnames(cluster_counts.tmp) == rownames(cluster_metadata.tmp)))

     # Create DESeq2 object
     dds <- DESeqDataSetFromMatrix(cluster_counts.tmp,
                                   colData = cluster_metadata.tmp,
                                   design = ~ Site)
     PCA <- dds %>% vst %>% plotPCA(intgroup=c("Site", "Rep"), returnData=TRUE)
     percentVar <- round(100 * attr(PCA, "percentVar"))
     PCA <- ggplot(PCA, aes(PC1, PC2, color=Site, shape=Rep)) +
         geom_point(size=3) +
         xlab(paste0("PC1: ",percentVar[1],"% variance")) +
         ylab(paste0("PC2: ",percentVar[2],"% variance")) +
         coord_fixed()+
         ggtitle(paste("No Minimum Filtering",unique(cluster_metadata.tmp$Day),unique(cluster_metadata.tmp$cluster_id)))
          dds <- estimateSizeFactors(dds)
          # 5 reads in 3 samples
          keep <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3
          dds <- dds[keep,]
    PCA.wmin <- dds %>% vst %>% plotPCA(intgroup=c("Site", "Rep"), returnData=TRUE)
    percentVar <- round(100 * attr(PCA.wmin, "percentVar"))
    PCA.wmin <- ggplot(PCA.wmin, aes(PC1, PC2, color=Site, shape=Rep)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) +
        coord_fixed()+ggtitle(paste("Minimum 5 Reads 3 Samples",unique(cluster_metadata.tmp$Day),unique(cluster_metadata.tmp$cluster_id)))
              dds.d = DESeq(dds, test = "LRT", reduced = ~1)
          #for anova
          res = results(dds.d)
          #res.ph.UNCCHOP = results(dds.d, contrast = list("Site_UNC_vs_CHOP"))
          #res.ph.CNCHOP = results(dds.d, contrast = list("Site_CN_vs_CHOP"))
          #res.ph.CNUNC = results(dds.d, contrast = list("Site_CN_vs_CHOP",
          #                                            "Site_UNC_vs_CHOP"))
          #results(dds, contrast=c("Site","UNC","CN"))
          nameforDeseq2 <- as.character(cluster_metadata.tmp[1,1])
          #write out table from DESeq2 no need to overwrite other files
          write.csv(res, file= paste0(results.here,nameforDeseq2,".csv"))
          #write.csv(res, file= paste0(results.here,nameforDeseq2,".csv"))
          #write.csv(res.ph.UNCCHOP, file= paste0(results.here,"res.ph.UNCCHOP",nameforDeseq2,".csv"))
          #write.csv(res.ph.CNCHOP, file= paste0(results.here,"res.ph.CNCHOP",nameforDeseq2,".csv"))
          #write.csv(res.ph.CNUNC, file= paste0(results.here,"res.ph.CNUNC",nameforDeseq2,".csv"))
          pdf(paste(plot.here,"July19_IDDRC_DESeq2_LRT_",unique(cluster_metadata.tmp$Day),unique(cluster_metadata.tmp$cluster_id),".pdf"), onefile = TRUE, width = 12, height = 8)
          print(PCA)
          print(PCA.wmin)
          dev.off()
  }
  i = i+1
}
