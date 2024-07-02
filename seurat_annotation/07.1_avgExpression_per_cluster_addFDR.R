# Get average expression for each cluster from Kreigsten

library(Seurat)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
#library(here)
library(ComplexHeatmap)

#library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# directory for average expression across clusters at different resolutions
here <- "/work/users/r/o/roseg/single-cell_reproducibility/"
dir.expression <- paste0(here,"results/seurat/cluster_expression/")
dir.create(dir.expression, recursive = TRUE, showWarnings = FALSE)

dir.pdf <- paste0(here,"doc/seurat/pdf/")



# Kreigstein average expression (also for IMPORT)
kreig.avg.expression.data.output.rds <- paste0(here,"results/seurat/cluster_expression/20230321_kreigstein_avgExpression_SCT_Data.rds")
kreig.avg.expression.scaleData.output.rds <- paste0(here,"results/seurat/cluster_expression/20230321_kreigstein_avgExpression_SCT_scaleData.rds")


# INPUT FILES ##########################################################################################################
# integrated and clustered dataset
seurat.integrated.rds <- paste0(here,"results/seurat/20230320_PGP1_UNCR3_reference_integrated_clustered_seurat_object.rds")

# integrated and clustered dataset, rescaled seurat object
seurat.rescale.rds <- paste0(here,"results/seurat/20230321_PGP1_UNCR3_reference_integrated_SCTRescale_seurat_object.rds")

# kreigstein seurat object
seurat.kreig.rds <- "/proj/steinlab/projects/IVIV_scRNA/youngsook_pine/extData_annotation/kreigstein_primary_seuratObject_SCTv2/Kreigstein_primary_integrated.rds"

# Kreigstein average expression: 42 clusters, 5000 genes
#kreig.avg.expression.output.rds <- here("results/seurat/20230306_kreigstein_avgExpression_SCT_scaleData.rds")

# Kreigstein (Bhaduri) cell types
kreig.cell.type.txt <- "/proj/steinlab/projects/IVIV_scRNA/youngsook_pine/extData_annotation/markerList/markerList/bhaduri_primary.txt"

# gencode gene names, types, and ensgid
df.gencode.csv <- here("data/refgenome/gencode/gencode.v40.genes.csv.gz")

# marker genes
df.markers.csv <- here("data/marker_genes/marker_genes.csv")

# GLOBALS ##############################################################################################################

# Load Bhaduri Data ############
df.avgExp.bhaduri.data <- readRDS(kreig.avg.expression.data.output.rds)
df.avgExp.bhaduri.scaleData <- readRDS(kreig.avg.expression.scaleData.output.rds)

cell.info <- data.table::fread(kreig.cell.type.txt, data.table=F, header=T)
cell.info$cluster <- paste0('bhaduri_C',cell.info$V1)

# Load PGP1 Integrated Data ############
printMessage("Loading Seurat Object...")
seur.pgp1 <- readRDS(seurat.integrated.rds)
printMessage("Finished Loading.")

# get cluster resolution column names
cluster.colnames <- colnames(seur.pgp1@meta.data)[grepl("integrated_snn_res", colnames(seur.pgp1@meta.data))]

pdf(paste0(dir.pdf, "20230321_PGP1_UNCR3_reference_bhaduri_comparison.pdf"), height = 8, width = 12)

# ONLY Res0.6 now instead of loop over cluster resolutions
#for (i in 1:length(cluster.colnames)) {
i=4
    printMessage(paste("Plotting PGP1 v Bhaduri for resolution:", cluster.colnames[i]))

    # set clusters by selecting a resolution and setting idents
    Idents(seur.pgp1) <- seur.pgp1@meta.data[,match(cluster.colnames[i], colnames(seur.pgp1@meta.data))]
    print("Cluster idents:")
    print(summary(Idents(seur.pgp1)))

    # get avg expression
    df.avgExp.pgp1.data <- as.data.frame(AverageExpression(seur.pgp1, assays = "SCT", slot = "data", group.by = "ident")[[1]])
    colnames(df.avgExp.pgp1.data) <- sub('', 'pgp1_C', colnames(df.avgExp.pgp1.data))

    df.avgExp.pgp1.scaleData <- as.data.frame(AverageExpression(seur.pgp1, assays = "SCT", slot = "scale.data", group.by = "ident")[[1]])
    colnames(df.avgExp.pgp1.scaleData) <- sub('', 'pgp1_C', colnames(df.avgExp.pgp1.scaleData))


    # Scale.Data Slot

    # list of intersecting genes
    genelist <- intersect(rownames(df.avgExp.pgp1.scaleData), rownames(df.avgExp.bhaduri.scaleData))


    #calculates the pearson correlation but not significance values
    avgExp.bhaduri.tmp <- df.avgExp.bhaduri.scaleData[match(genelist, rownames(df.avgExp.bhaduri.scaleData)),]
    avgExp.pgp1.tmp <- df.avgExp.pgp1.scaleData[match(genelist, rownames(df.avgExp.pgp1.scaleData)),]
    cor.res.tmp <- cor(avgExp.bhaduri.tmp, avgExp.pgp1.tmp, method = "pearson") %>% as.matrix()
    print(max(cor.res.tmp))

    cor.res.tmp <- cor.res.tmp[match(cell.info$cluster, rownames(cor.res.tmp)),]
    cell.class <- cell.info$Class
    cell.state <- cell.info$State
    cell.type <- cell.info$Type
    cell.subtype<-cell.info$Subtype

    cell.class.color<-setNames(c('#5cc4c3','#2a4f72','#d1d1d1'),unique(cell.class))
    cell.state.color<-setNames(c('#dd7373','#3b3561','#ead94c','#d1d1d1'),unique(cell.state))
    cell.type.color<-setNames(c('#f75449','#2509b5','#6e88f5','#8dfab3','#23aa4e','#3fe72d','#fdfe61','#b29415','#f75449','#d1d1d1'),unique(cell.type))

    #update Bhaduri's data
    rownames(cor.res.tmp) <- paste(cell.class,cell.state,cell.type,cell.subtype,sep=' ')
    rownames(cor.res.tmp) <- sub(' Outlier Outlier Outlier','',rownames(cor.res.tmp))

#add FDR corrected pval with row x column names
    #calculate the p-vales for person correlations, save as vector, then performFDR correction
    Pearson.pval <- c("Call","Pval","Cor")
    #run for each PGP1 cluster
    for(x in 1:length(avgExp.pgp1.tmp)){
    for (i in 1:length(df.avgExp.bhaduri.scaleData)){
        cor.sigtest.tmp <- cor.test(avgExp.bhaduri.tmp[ ,i], avgExp.pgp1.tmp[ ,x], method = "pearson")
        #create table with BhaduriName, PGP1Name,
        Pearson.pval.t <- c(paste(colnames(avgExp.bhaduri.tmp[i]),colnames(avgExp.pgp1.tmp[x])),
                            cor.sigtest.tmp$p.value, cor.sigtest.tmp$estimate)
        Pearson.pval <- rbind(Pearson.pval,Pearson.pval.t)
    }
    }
    Pearson.pval <- as_tibble(Pearson.pval)
    colnames(Pearson.pval) <- Pearson.pval[1,]
    Pearson.pval <- Pearson.pval[-1,]
    Pearson.pval <- Pearson.pval %>% mutate_at(c('Pval', 'Cor'), as.numeric)

    #adjust p-val on number of correlations run
    #you'll want all the comparisons runs as vector for this function
    Pearson.pval <- mutate(Pearson.pval, FDR = p.adjust(Pval, method="BH"))



#add Heat map annotations for sample cell counts
    #cell.pc.perDay<-readRDS(paste0('Seurat/post-processed/UMAP/20230109.keep/RefN_2_0fv8qgdnjs/10_10/',DFiltering,'.UMAP.cell.embeddings.clusterinfo_res',res,'.rds'))
    cell.N <- seur.pgp1@meta.data %>% group_by(.dots=cluster.colnames[i]) %>% summarize(n=n()) %>%pull(n)
    #donor.N <- seur.pgp1@meta.data %>% group_by(Site, Day, Rep, seurat_clusters) %>% summarize(n=n()) %>% group_by(seurat_clusters, Day)%>%summarize(n=n())
    #seur.pgp1@meta.data %>% group_by(Day) %>% summarize(n=n()) %>% left_join(seur.pgp1@meta.data %>% group_by(seurat_clusters,Day)%>% summarize(nc=n())) %>% data.frame() %>% mutate(pc = 100*nc/n) %>% arrange(seurat_clusters) -> cell.pc.perDay

    #seur.pgp1@meta.data$pc <- perc

    # ha = HeatmapAnnotation( 'Cell proportion' = anno_barplot(matrix(nc = 2, c(seur.pgp1@meta.data %>% filter(Day=='D14') %>% pull(pc), seur.pgp1@meta.data %>% filter(Day=='D84') %>% pull(pc)),dimnames=list(unique(seur.pgp1@meta.data$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
    #                         'Donor #' = anno_barplot(matrix(nc = 2, c(donor.N %>% filter(Day=='D14') %>% pull(n), donor.N %>% filter(Day=='D84') %>% pull(n)),dimnames=list(unique(donor.N$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
    #                         'Cell counts' = anno_barplot(cell.N,gp=gpar(fill='darkgreen')))
    #
    ha = HeatmapAnnotation('Cell counts' = anno_barplot(cell.N,gp=gpar(fill='darkgreen')))

    #Heatmap(cor.res.tmp)

    p <- Heatmap(cor.res.tmp,name="Pearson's correlation",
            top_annotation=ha,
            right_annotation=rowAnnotation(Class=cell.class,state=cell.state,type=cell.type,
                                           col=list(Class=cell.class.color,state=cell.state.color,type=cell.type.color)),
            column_title = paste0(cluster.colnames[i], " scale.data pearson"),row_names_gp=gpar(fontsize=6),column_names_gp=gpar(fontsize=6),
            cell_fun = function(j, i, x, y, width, height, fill) {
                if(cor.res.tmp[i, j] > 0.3)
                    grid.text(sprintf("%.2f", cor.res.tmp[i, j]), x, y, gp = gpar(fontsize = 6))
            }
    )

    print(p)

    cor.res.tmp <- cor(avgExp.bhaduri.tmp, avgExp.pgp1.tmp, method = "spearman") %>% as.matrix()
    print(max(cor.res.tmp))

    cor.res.tmp <- cor.res.tmp[match(cell.info$cluster, rownames(cor.res.tmp)),]
    cell.class <- cell.info$Class
    cell.state <- cell.info$State
    cell.type <- cell.info$Type
    cell.subtype<-cell.info$Subtype

    cell.class.color<-setNames(c('#5cc4c3','#2a4f72','#d1d1d1'),unique(cell.class))
    cell.state.color<-setNames(c('#dd7373','#3b3561','#ead94c','#d1d1d1'),unique(cell.state))
    cell.type.color<-setNames(c('#f75449','#2509b5','#6e88f5','#8dfab3','#23aa4e','#3fe72d','#fdfe61','#b29415','#f75449','#d1d1d1'),unique(cell.type))

    #update Bhaduri's data
    rownames(cor.res.tmp) <- paste(cell.class,cell.state,cell.type,cell.subtype,sep=' ')
    rownames(cor.res.tmp) <- sub(' Outlier Outlier Outlier','',rownames(cor.res.tmp))

    #cell.pc.perDay<-readRDS(paste0('Seurat/post-processed/UMAP/20230109.keep/RefN_2_0fv8qgdnjs/10_10/',DFiltering,'.UMAP.cell.embeddings.clusterinfo_res',res,'.rds'))
    cell.N <- seur.pgp1@meta.data %>% group_by(.dots=cluster.colnames[i]) %>% summarize(n=n()) %>%pull(n)
    #donor.N <- seur.pgp1@meta.data %>% group_by(Site, Day, Rep, seurat_clusters) %>% summarize(n=n()) %>% group_by(seurat_clusters, Day)%>%summarize(n=n())
    #seur.pgp1@meta.data %>% group_by(Day) %>% summarize(n=n()) %>% left_join(seur.pgp1@meta.data %>% group_by(seurat_clusters,Day)%>% summarize(nc=n())) %>% data.frame() %>% mutate(pc = 100*nc/n) %>% arrange(seurat_clusters) -> cell.pc.perDay

    #seur.pgp1@meta.data$pc <- perc

    # ha = HeatmapAnnotation( 'Cell proportion' = anno_barplot(matrix(nc = 2, c(seur.pgp1@meta.data %>% filter(Day=='D14') %>% pull(pc), seur.pgp1@meta.data %>% filter(Day=='D84') %>% pull(pc)),dimnames=list(unique(seur.pgp1@meta.data$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
    #                         'Donor #' = anno_barplot(matrix(nc = 2, c(donor.N %>% filter(Day=='D14') %>% pull(n), donor.N %>% filter(Day=='D84') %>% pull(n)),dimnames=list(unique(donor.N$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
    #                         'Cell counts' = anno_barplot(cell.N,gp=gpar(fill='darkgreen')))
    #
    ha = HeatmapAnnotation('Cell counts' = anno_barplot(cell.N,gp=gpar(fill='darkgreen')))

    #Heatmap(cor.res.tmp)

    p <- Heatmap(cor.res.tmp,name="Spearman's correlation",
            top_annotation=ha,
            right_annotation=rowAnnotation(Class=cell.class,state=cell.state,type=cell.type,
                                           col=list(Class=cell.class.color,state=cell.state.color,type=cell.type.color)),
            column_title = paste0(cluster.colnames[i], " scale.data spearman"),row_names_gp=gpar(fontsize=6),column_names_gp=gpar(fontsize=6),
            cell_fun = function(j, i, x, y, width, height, fill) {
                if(cor.res.tmp[i, j] > 0.3)
                    grid.text(sprintf("%.2f", cor.res.tmp[i, j]), x, y, gp = gpar(fontsize = 6))
            }
    )

    print(p)

    # Data Slot

    # list of intersecting genes
    genelist <- intersect(rownames(df.avgExp.pgp1.data), rownames(df.avgExp.bhaduri.data))


    avgExp.bhaduri.tmp <- df.avgExp.bhaduri.data[match(genelist, rownames(df.avgExp.bhaduri.data)),]
    avgExp.pgp1.tmp <- df.avgExp.pgp1.data[match(genelist, rownames(df.avgExp.pgp1.data)),]
    cor.res.tmp <- cor(avgExp.bhaduri.tmp, avgExp.pgp1.tmp, method = "pearson") %>% as.matrix()
    print(max(cor.res.tmp))

    cor.res.tmp <- cor.res.tmp[match(cell.info$cluster, rownames(cor.res.tmp)),]
    cell.class <- cell.info$Class
    cell.state <- cell.info$State
    cell.type <- cell.info$Type
    cell.subtype<-cell.info$Subtype

    cell.class.color<-setNames(c('#5cc4c3','#2a4f72','#d1d1d1'),unique(cell.class))
    cell.state.color<-setNames(c('#dd7373','#3b3561','#ead94c','#d1d1d1'),unique(cell.state))
    cell.type.color<-setNames(c('#f75449','#2509b5','#6e88f5','#8dfab3','#23aa4e','#3fe72d','#fdfe61','#b29415','#f75449','#d1d1d1'),unique(cell.type))

    #update Bhaduri's data
    rownames(cor.res.tmp) <- paste(cell.class,cell.state,cell.type,cell.subtype,sep=' ')
    rownames(cor.res.tmp) <- sub(' Outlier Outlier Outlier','',rownames(cor.res.tmp))

    #cell.pc.perDay<-readRDS(paste0('Seurat/post-processed/UMAP/20230109.keep/RefN_2_0fv8qgdnjs/10_10/',DFiltering,'.UMAP.cell.embeddings.clusterinfo_res',res,'.rds'))
    cell.N <- seur.pgp1@meta.data %>% group_by(.dots=cluster.colnames[i]) %>% summarize(n=n()) %>%pull(n)
    #donor.N <- seur.pgp1@meta.data %>% group_by(Site, Day, Rep, seurat_clusters) %>% summarize(n=n()) %>% group_by(seurat_clusters, Day)%>%summarize(n=n())
    #seur.pgp1@meta.data %>% group_by(Day) %>% summarize(n=n()) %>% left_join(seur.pgp1@meta.data %>% group_by(seurat_clusters,Day)%>% summarize(nc=n())) %>% data.frame() %>% mutate(pc = 100*nc/n) %>% arrange(seurat_clusters) -> cell.pc.perDay

    #seur.pgp1@meta.data$pc <- perc

    # ha = HeatmapAnnotation( 'Cell proportion' = anno_barplot(matrix(nc = 2, c(seur.pgp1@meta.data %>% filter(Day=='D14') %>% pull(pc), seur.pgp1@meta.data %>% filter(Day=='D84') %>% pull(pc)),dimnames=list(unique(seur.pgp1@meta.data$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
    #                         'Donor #' = anno_barplot(matrix(nc = 2, c(donor.N %>% filter(Day=='D14') %>% pull(n), donor.N %>% filter(Day=='D84') %>% pull(n)),dimnames=list(unique(donor.N$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
    #                         'Cell counts' = anno_barplot(cell.N,gp=gpar(fill='darkgreen')))
    #
    ha = HeatmapAnnotation('Cell counts' = anno_barplot(cell.N,gp=gpar(fill='darkgreen')))

    #Heatmap(cor.res.tmp)

    p <- Heatmap(cor.res.tmp,name="Pearson's correlation",
            top_annotation=ha,
            right_annotation=rowAnnotation(Class=cell.class,state=cell.state,type=cell.type,
                                           col=list(Class=cell.class.color,state=cell.state.color,type=cell.type.color)),
            column_title = paste0(cluster.colnames[i], " data pearson"),row_names_gp=gpar(fontsize=6),column_names_gp=gpar(fontsize=6),
            cell_fun = function(j, i, x, y, width, height, fill) {
                if(cor.res.tmp[i, j] > 0.3)
                    grid.text(sprintf("%.2f", cor.res.tmp[i, j]), x, y, gp = gpar(fontsize = 6))
            }
    )

    print(p)

    cor.res.tmp <- cor(avgExp.bhaduri.tmp, avgExp.pgp1.tmp, method = "spearman") %>% as.matrix()
    print(max(cor.res.tmp))

    cor.res.tmp <- cor.res.tmp[match(cell.info$cluster, rownames(cor.res.tmp)),]
    cell.class <- cell.info$Class
    cell.state <- cell.info$State
    cell.type <- cell.info$Type
    cell.subtype<-cell.info$Subtype

    cell.class.color<-setNames(c('#5cc4c3','#2a4f72','#d1d1d1'),unique(cell.class))
    cell.state.color<-setNames(c('#dd7373','#3b3561','#ead94c','#d1d1d1'),unique(cell.state))
    cell.type.color<-setNames(c('#f75449','#2509b5','#6e88f5','#8dfab3','#23aa4e','#3fe72d','#fdfe61','#b29415','#f75449','#d1d1d1'),unique(cell.type))

    #update Bhaduri's data
    rownames(cor.res.tmp) <- paste(cell.class,cell.state,cell.type,cell.subtype,sep=' ')
    rownames(cor.res.tmp) <- sub(' Outlier Outlier Outlier','',rownames(cor.res.tmp))

    #cell.pc.perDay<-readRDS(paste0('Seurat/post-processed/UMAP/20230109.keep/RefN_2_0fv8qgdnjs/10_10/',DFiltering,'.UMAP.cell.embeddings.clusterinfo_res',res,'.rds'))
    cell.N <- seur.pgp1@meta.data %>% group_by(.dots=cluster.colnames[i]) %>% summarize(n=n()) %>%pull(n)
    #donor.N <- seur.pgp1@meta.data %>% group_by(Site, Day, Rep, seurat_clusters) %>% summarize(n=n()) %>% group_by(seurat_clusters, Day)%>%summarize(n=n())
    #seur.pgp1@meta.data %>% group_by(Day) %>% summarize(n=n()) %>% left_join(seur.pgp1@meta.data %>% group_by(seurat_clusters,Day)%>% summarize(nc=n())) %>% data.frame() %>% mutate(pc = 100*nc/n) %>% arrange(seurat_clusters) -> cell.pc.perDay

    #seur.pgp1@meta.data$pc <- perc

    # ha = HeatmapAnnotation( 'Cell proportion' = anno_barplot(matrix(nc = 2, c(seur.pgp1@meta.data %>% filter(Day=='D14') %>% pull(pc), seur.pgp1@meta.data %>% filter(Day=='D84') %>% pull(pc)),dimnames=list(unique(seur.pgp1@meta.data$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
    #                         'Donor #' = anno_barplot(matrix(nc = 2, c(donor.N %>% filter(Day=='D14') %>% pull(n), donor.N %>% filter(Day=='D84') %>% pull(n)),dimnames=list(unique(donor.N$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
    #                         'Cell counts' = anno_barplot(cell.N,gp=gpar(fill='darkgreen')))
    #
    ha = HeatmapAnnotation('Cell counts' = anno_barplot(cell.N,gp=gpar(fill='darkgreen')))

    #Heatmap(cor.res.tmp)

    p <- Heatmap(cor.res.tmp,name="Spearman's correlation",
            top_annotation=ha,
            right_annotation=rowAnnotation(Class=cell.class,state=cell.state,type=cell.type,
                                           col=list(Class=cell.class.color,state=cell.state.color,type=cell.type.color)),
            column_title = paste0(cluster.colnames[i], " data spearman"),row_names_gp=gpar(fontsize=6),column_names_gp=gpar(fontsize=6),
            cell_fun = function(j, i, x, y, width, height, fill) {
                if(cor.res.tmp[i, j] > 0.3)
                    grid.text(sprintf("%.2f", cor.res.tmp[i, j]), x, y, gp = gpar(fontsize = 6))
            }
    )

    print(p)

#}

dev.off()


