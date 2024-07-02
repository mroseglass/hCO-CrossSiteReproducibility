library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(ggpubr)
library(plotrix)
library(gprofiler2)
library(ggrepel)
library(EnhancedVolcano)
library(gprofiler2)
library(ggplot2)

# load in .csv from DESeq2 table
#INPUTS ##############################################################################################
db.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/"
#deg.here <- "/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/DESeq2/Min5read3xpt/"
#deg.here <- "/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/DESeq2/CNCHOP.ph.min/"
#deg.here <- "/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/DESeq2/CNUNC.ph.min/"
#deg.here <- "/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/DESeq2/UNCCHOP.ph.min/"
#deg.here <- "/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/DESeq2/LTRtest_reduced1_july1823_days/"
deg.here <- "/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/DESeq2/LTRtest_reduced1_july1823/"
plot.here <- "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/DESeq2/LRT_Day/"

CellTypeID <- read_csv(paste0(db.here, "ClustertoCellType.csv"))
CellTypeID$Cluster <- as.character(CellTypeID$Cluster)
# GLOBALS ##########################################################################################
Site.color <-c("#12783D", "#882155", "#4040C6")
Replicate.shapes <- c(0,1,2,3,4,15,16,17)
mixing_pallete <-c('#5F95B2','#BCC6E5','#B09AB1','#F1C4DC','#EA9F8B',
                   '#F8C893','#89A48C','#369F48','#DB7B87',
                   '#E12228','#B177B3','#2179B4','#F47B20','#F89B40','#F15A29')

#####################################################################################################
#Example opening of one file
List.deg.expt <- c(list.files(deg.here))
i=1
deseq2.result <- read_csv(paste0(deg.here,List.deg.expt[i]))
#drop NA values
deseq2.result <- deseq2.result%>% drop_na(`padj`)
#add name of test to table
filename <- gsub ("_R.*","",List.deg.expt[i])
filename <- gsub ("CHOP","",filename)
Add.DEGtest <- rep(c(filename),each=length(deseq2.result$log2FoldChange))
IDDRC.deg <- cbind(deseq2.result,Add.DEGtest)

#Add rest of files to same object
for (i in 2:length(List.deg.expt)){
  deseq2.result <- read_csv(paste0(deg.here,List.deg.expt[i]))
  #drop NA values
  deseq2.result <- deseq2.result%>% drop_na(`padj`)
  #add name of test to table
  filename <- gsub ("_R.*","",List.deg.expt[i])
  filename <- gsub ("CHOP","",filename)
  Add.DEGtest <- rep(c(filename),each=length(deseq2.result$log2FoldChange))
  CorrectFile2 <- cbind(deseq2.result,Add.DEGtest)
  IDDRC.deg <- rbind(CorrectFile2,IDDRC.deg)
}

#for recycle the scripts
#IDDRC.deg$Add.DEGtest <- gsub("Minimum5readsin3samples","",IDDRC.deg$Add.DEGtest)

#only 2 cluster/day combos with DEG surviving multiple correction testing
IDDRC.deg.005 <- filter(IDDRC.deg, padj < 0.05)
#write_csv(IDDRC.deg.005, file = paste0(db.here,"/IDDRC_LTRReduced1_alpha005padj_Min5readsin3samples_July17.csv"))
#IDDRC.deg.005<- read_csv(paste0(db.here,"/DEG/IDDRC_UNCCHOP.ph.min_alpha005padj_Min5readsin3samples_June17.csv"))

#create barchart of all DEG by Day
IDDRC.deg.005 <- mutate(IDDRC.deg.005, Cluster = gsub("._.*","",IDDRC.deg.005$Add.DEGtest))
#IDDRC.deg.005 <- mutate(IDDRC.deg.005, Cluster = gsub("res.ph.UNC","",IDDRC.deg.005$Cluster))
IDDRC.deg.005 <- mutate(IDDRC.deg.005, Day = gsub(".*_","",IDDRC.deg.005$Add.DEGtest))
IDDRC.deg.005.D14 <- filter(IDDRC.deg.005, Day == "D14")
D14.count <- IDDRC.deg.005.D14 %>% dplyr::count(Add.DEGtest)
IDDRC.deg.005.D14 <- inner_join(IDDRC.deg.005.D14, D14.count)
IDDRC.deg.005.D14 <- full_join(CellTypeID,IDDRC.deg.005.D14)
IDDRC.deg.005.D14 <- IDDRC.deg.005.D14 %>% replace_na(list(n=0))

IDDRC.deg.005.D84 <- filter(IDDRC.deg.005, Day == "D84")
D84.count <- IDDRC.deg.005.D84 %>% dplyr::count(Add.DEGtest)
IDDRC.deg.005.D84 <- inner_join(IDDRC.deg.005.D84, D84.count)
IDDRC.deg.005.D84 <- inner_join(IDDRC.deg.005.D84, CellTypeID)

IDDRC.deg.005.D14$CellType <- factor(IDDRC.deg.005.D14$CellType,
                       levels=c("Neuroepithelial stem cells","Dividing neural progenitor cells, S","Dividing neural progenitor cells, G2",
                                "Medial Pallium/Marginal Zone","Dividing intermediate progenitors, S","Intermediate progenitors",
                                "Radial glia","Outer radial glia","Lower layer neuron a",
                                "Lower layer neuron b","Upper layer neuron a, Layer 2/3/4",
                                "Upper layer neuron b, layer 2/3","Pan cortical neuron","Pan neuron/Cajal - Retzius",
                                "Unspecified neuron"))


ggplot(data = IDDRC.deg.005.D14, aes(y=CellType,
                             fill = CellType)) +
    geom_bar(position = position_dodge(preserve = "single")) +
    scale_fill_manual(values= mixing_pallete) +
    scale_y_discrete(limits=rev)+
    theme_bw() +
    ylab("") +
    xlab("") + theme(legend.position = "none")+
    ggtitle("Day 14 Count Differentially Expressed Genes")

IDDRC.deg.005.D84$CellType <- factor(IDDRC.deg.005.D84$CellType,
                                     levels=c("Neuroepithelial stem cells","Dividing neural progenitor cells, S","Dividing neural progenitor cells, G2",
                                              "Medial Pallium/Marginal Zone","Dividing intermediate progenitors, S","Intermediate progenitors",
                                              "Radial glia","Outer radial glia","Lower layer neuron a",
                                              "Lower layer neuron b","Upper layer neuron a, Layer 2/3/4",
                                              "Upper layer neuron b, layer 2/3","Pan cortical neuron","Pan neuron/Cajal - Retzius",
                                              "Unspecified neuron"))


ggplot(data = IDDRC.deg.005.D84, aes(y=CellType,
                                     fill = CellType)) +
    geom_bar(position = position_dodge(preserve = "single")) +
    scale_fill_manual(values= mixing_pallete) +
    scale_y_discrete(limits=rev)+
    theme_bw() +
    ylab("") +
    xlab("") + theme(legend.position = "none")+
    ggtitle("Day 84 Count Differentially Expressed Genes")

# GO:BP annotation of clusters with DEG #########################################
IDDRC.deg.005 <- inner_join(IDDRC.deg.005, CellTypeID)
IDDRC.deg.005 <- dplyr::rename(IDDRC.deg.005, gene = `...1`)
IDDRC.deg.005.D14 <- dplyr::rename(IDDRC.deg.005.D14, gene = `...1`)
IDDRC.deg.005.D84 <- dplyr::rename(IDDRC.deg.005.D84, gene = `...1`)

#run in 2 separate instances for D14 and D84
#Go D84 #####################################################################################
list.cluster.deg <- c(unique(IDDRC.deg.005.D84$CellType))
output <- list()
i=1
DEG.to.plot <- filter(IDDRC.deg.005.D84, CellType == unlist(list.cluster.deg[i]))
# Annotation of DEG in each cluster
gostres <- gost(query = DEG.to.plot$gene,
                organism = "hsapiens", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = FALSE,
                user_threshold = 0.05, correction_method = "g_SCS",
                domain_scope = "annotated", custom_bg = NULL,
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
DEG.ann <- gostres$result
CellType <- rep(list.cluster.deg[i],length(DEG.ann$p_value))
DEG.ann$query <- CellType
DEG.ann.all <- DEG.ann
output[[i]] <- DEG.ann.all

for (i in 2:length(list.cluster.deg)){
DEG.to.plot <- filter(IDDRC.deg.005.D84, CellType == unlist(list.cluster.deg[i]))
# Annotation of DEG in each cluster
gostres <- gost(query = DEG.to.plot$gene,
                organism = "hsapiens", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = FALSE,
                user_threshold = 0.05, correction_method = "g_SCS",
                domain_scope = "annotated", custom_bg = NULL,
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
DEG.ann <- gostres$result
CellType <- rep(list.cluster.deg[i],length(DEG.ann$p_value))
DEG.ann$query <- CellType
output[[i]] <- DEG.ann.all
DEG.ann.all <- rbind(DEG.ann.all,DEG.ann)
}
#saveRDS(output, paste0(db.here, "July18_LTRreduced1_Go_Annotation_Output_Day84_DEG_.rds"))
#write_csv(DEG.ann.all, paste0(db.here,"July18_LTRreduced1_Go_Annotation_Output_Day84_Output_ANOVAallsite.csv"))
DEG.ann.allKeep <- DEG.ann.all
DEG.ann.all <- DEG.ann.allKeep
#barcharts color by pval showing number of genes in category ####################################################
DEG.ann.all <- filter(DEG.ann.all, source == "GO:BP")
DEG.ann.all.r <- mutate(DEG.ann.all, negativelog10 = -log10(p_value) )
DEG.ann.all.r <- DEG.ann.all.r %>%
    arrange(desc(negativelog10)) %>%
    group_by(query) %>%
    slice(1:3)

#try a dotplot instead
ggplot(data = DEG.ann.all.r, aes(x=reorder(query,-negativelog10), y = reorder(term_name,negativelog10))) +
    geom_point(aes(color = query,size=negativelog10)) +
    scale_size(limits = c(0,16), breaks=c(0,4,8,12,16)) +
    theme_minimal()+ ylab("GO Biological Process") + xlab("Cell Type")+
    ggtitle("Day 84 GO Annotation")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#ggsave("D84GO_Dotplot_June28.pdf",plot = last_plot(), width = 9, height = 8, units = "in",
#       dpi = 300, device = "pdf",
#       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/DESeq2/")

#Save GO D14 ########################################################################################
IDDRC.deg.005.D14 <- IDDRC.deg.005.D14 %>% drop_na(gene)
list.cluster.deg <- c(unique(IDDRC.deg.005.D14$CellType))
output <- list()
i=1
DEG.to.plot <- filter(IDDRC.deg.005.D14, CellType == unlist(list.cluster.deg[i]))
# Annotation of DEG in each cluster
gostres <- gost(query = DEG.to.plot$gene,
                organism = "hsapiens", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = FALSE,
                user_threshold = 0.05, correction_method = "g_SCS",
                domain_scope = "annotated", custom_bg = NULL,
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
DEG.ann <- gostres$result
CellType <- rep(list.cluster.deg[i],length(DEG.ann$p_value))
DEG.ann$query <- CellType
DEG.ann.all <- DEG.ann
output[[i]] <- DEG.ann.all

for (i in 2:length(list.cluster.deg)){
    DEG.to.plot <- filter(IDDRC.deg.005.D14, CellType == unlist(list.cluster.deg[i]))
    # Annotation of DEG in each cluster
    gostres <- gost(query = DEG.to.plot$gene,
                    organism = "hsapiens", ordered_query = FALSE,
                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                    measure_underrepresentation = FALSE, evcodes = FALSE,
                    user_threshold = 0.05, correction_method = "g_SCS",
                    domain_scope = "annotated", custom_bg = NULL,
                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
    DEG.ann <- gostres$result
    CellType <- rep(list.cluster.deg[i],length(DEG.ann$p_value))
    DEG.ann$query <- CellType
    output[[i]] <- DEG.ann.all
    DEG.ann.all <- rbind(DEG.ann.all,DEG.ann)
}

#saveRDS(output, paste0(db.here, "July18_LTRreduced1_Go_Annotation_Output_Day14_DEG.rds"))
#write_csv(DEG.ann.all, paste0(db.here,"July18_LTRreduced1_Go_Annotation_D14vD84.csv"))

DEG.ann.allKeep <- DEG.ann.all
DEG.ann.all <- DEG.ann.allKeep
#barcharts color by pval showing number of genes in category ####################################################
DEG.ann.all <- filter(DEG.ann.all, source == "GO:BP")
DEG.ann.all.r <- mutate(DEG.ann.all, negativelog10 = -log10(p_value) )
DEG.ann.all.r <- DEG.ann.all.r %>%
    arrange(desc(negativelog10)) %>%
    group_by(query) %>%
    slice(1:3)

#try a dotplot instead
ggplot(data = DEG.ann.all.r, aes(x=reorder(query,-negativelog10), y = reorder(term_name,negativelog10))) +
    geom_point(aes(color = query,size=negativelog10)) +
    scale_size(limits = c(0,16), breaks=c(0,4,8,12,16)) +
    theme_minimal()+ ylab("GO Biological Process") + xlab("Cell Type")+
    ggtitle("Day 14 GO Annotation")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#ggsave("D14GO_Dotplot_June28.pdf",plot = last_plot(), width = 8, height = 8, units = "in",
#       dpi = 300, device = "pdf",
#       path = "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/DESeq2/")

#check an example to see if the DEG make sense
gexpr <- readRDS("/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/VST.Expr.PerSample_Day_Cluster.rm.lowExpGenes.0.6.rds")

#####################################################################################
# Each gene for each day for each cluster in a pdf
dc.exp <- as_tibble(gexpr[[1]]$gexpr, rownames = NA)
dc.exp <- rownames_to_column(dc.exp, var = "gene") %>% as_tibble()
#IDDRC.deg.005 <- dplyr::rename(IDDRC.deg.005, gene = `...1`)
dc.exp <- inner_join(IDDRC.deg.005, dc.exp, by = "gene")
list.cluster.deg <- unique(IDDRC.deg.005$Add.DEGtest)
i=1
for(i in 2:length(list.cluster.deg)){
dc.exp.p <- filter(dc.exp, Add.DEGtest == unlist(list.cluster.deg[i]))
plot <- dplyr::select(dc.exp.p, !c(log2FoldChange,lfcSE,stat,pvalue,Add.DEGtest,padj,baseMean,Cluster,Day))
#plot all DESeq2 anova sig gene  by site
plot.adjp <- dplyr::select(dc.exp.p,c(gene,padj))
plot <- pivot_longer(plot, cols = !gene, names_to = "Experiment", values_to = "VSTexpression")
plot <- filter(plot, Experiment != "baseMean")
plot$Experiment <- gsub("0.","",plot$Experiment)
plot <- mutate(plot, Site = gsub("_.*","",Experiment))
plot <- mutate(plot, Rep = gsub(".*_","",Experiment))
plot.r <- full_join(plot.adjp,plot, by = "gene")
forlabels <- distinct(dplyr::select(plot.r, c(gene,padj)))
gene.list <- c(forlabels$gene)
padj.list <-c(forlabels$padj)

pdf(paste0(plot.here,"July18_DESeq2_LTRreduced1_OneWayANOVA_cluster",unlist(list.cluster.deg[i]),".pdf"), onefile = TRUE, width = 12, height = 8)
for(a in 1:length(gene.list)){
    plot.m <- filter(plot, gene == gene.list[a])

    gene.p <- ggplot(plot.m, aes(x=Site, y=VSTexpression, fill = Site)) +
        geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) +
        theme(axis.text.x = element_text(angle = 45)) +
        scale_fill_manual(values = Site.color) +
        scale_shape_manual(values=Replicate.shapes)+
        theme_classic()+
        geom_point(aes(shape = Rep),size = 4, position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
        ylab("VST normalized expression")+
        ggtitle(paste(gene.list[a],"DESeq2 adjusted p value",round(padj.list[a], digits = 6)))
    print(gene.p)

}
dev.off()
}



#SCRATCH #################
#stack barplots fo GO terms
#DEG.ann.all.r$query <- factor(DEG.ann.all.r$query,
#                        levels=c('Dividing intermediate progenitors, S',
#                                 'Dividing neural progenitor cells, G2',
#                                 'Intermediate progenitors',
#                                 'Lower layer neuron b',
#                                 'Outer radial glia',
#                                 'Unspecified neuron',
#                                 'Upper layer neuron a, Layer 2/3/4'))
ggplot(data = DEG.ann.all.r, aes(x=negativelog10, y = reorder(term_name,negativelog10),
                                 fill = query)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    ylab("-log10pval") +
    xlab("GO Biological Process") +
    ggtitle("GO enrichment analysis D84")

#previous ploting for each gene
d14.nes.exp <- as_tibble(gexpr[[1]]$gexpr, rownames = NA)
d14.nes.exp <- rownames_to_column(d14.nes.exp, var = "gene") %>% as_tibble()
D14.nes <- dplyr::rename(D14.nes, gene = `...1`)
d14.anova.gene <- inner_join(D14.nes, d14.nes.exp, by = "gene")
plot <- dplyr::select(d14.anova.gene, !c(log2FoldChange,lfcSE,stat,pvalue,Add.DEGtest,padj))
#plot all DESeq2 anova sig gene  by site
plot.adjp <- dplyr::select(d14.anova.gene,c(gene,padj))
plot <- pivot_longer(plot, cols = !gene, names_to = "Experiment", values_to = "VSTexpression")
plot <- filter(plot, Experiment != "baseMean")
plot$Experiment <- gsub("0.","",plot$Experiment)
plot <- mutate(plot, Site = gsub("_.*","",Experiment))
plot <- mutate(plot, Rep = gsub(".*_","",Experiment))
plot.r <- full_join(plot.adjp,plot, by = "gene")
forlabels <- distinct(select(plot.r, c(gene,padj)))
gene.list <- c(forlabels$gene)
padj.list <-c(forlabels$padj)

pdf(paste(plot.here,"June7_DESeq2_OneWayANOVA_D14",".pdf"), onefile = TRUE, width = 12, height = 8)
for(i in 1:length(gene.list)){
    plot.m <- filter(plot, gene == gene.list[i])

    gene.p <- (ggplot(plot.m, aes(x=Site, y=VSTexpression, fill = Site))) +
        geom_boxplot(outlier.shape = NA, lwd=0.2, show.legend = FALSE) +
        theme(axis.text.x = element_text(angle = 45)) +
        scale_fill_manual(values = Site.color) +
        scale_shape_manual(values=Replicate.shapes)+
        theme_classic()+
        #geom_jitter(width = 0.1)+
        geom_point(aes(shape = Rep),size = 2, position=position_jitterdodge(dodge.width=0.01),alpha=0.7) +
        ylab("VST normalized expression")+
      ggtitle(paste(gene.list[i],"DESeq2 adjusted p value",round(padj.list[i], digits = 3)))
    print(gene.p)

    }
dev.off()


#SCRATCH#
#Volcano plot for each cluster ##################################################################
#for D84 cluster 9
D84.panneuron.all <- dplyr::rename(D84.panneuron.all, gene = `...1`)
tolabel <- subset(D84.panneuron, padj %in% 0.007:0.002)
EnhancedVolcano(D84.panneuron.all,
                #lab = D84.panneuron.all$gene,
                x = 'basemean',
                y = 'padj',
                title = 'Day 84 "Unspecified neuron',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 4,
                xlim=c(-4,4),
                ylim=c(0,4),
                ylab = "-Log10 Adjusted P value",
                col=c('black', 'black', 'blue3', 'red3'),
                colAlpha = 1)

D14.nes.all <- dplyr::rename(D14.nes.all, gene = `...1`)
D14.labelset <- c(D14.nes.all$...1)
EnhancedVolcano(D14.nes.all,
                lab = D14.nes.all$gene,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Day 14 NSC',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 4,
                xlim=c(-4,4),
                ylim=c(0,4),
                ylab = "-Log10 Adjusted P value",
                col=c('black', 'black', 'blue3', 'red3'),
                colAlpha = 1)

D84.9 <- filter(IDDRC.deg.005, Add.DEGtest == "9._D84")
D14.9 <- filter(IDDRC.deg.005, Add.DEGtest == "9._D14")
D84.8 <- filter(IDDRC.deg.005, Add.DEGtest == "8._D84")
D14.8 <- filter(IDDRC.deg.005, Add.DEGtest == "8._D14")
D84.7 <- filter(IDDRC.deg.005, Add.DEGtest == "7._D84")
D84.6 <- filter(IDDRC.deg.005, Add.DEGtest == "6._D84")
D14.6 <- filter(IDDRC.deg.005, Add.DEGtest == "6._D14")
D84.5 <- filter(IDDRC.deg.005, Add.DEGtest == "5._D84")
D14.5 <- filter(IDDRC.deg.005, Add.DEGtest == "5._D14")
D84.4 <- filter(IDDRC.deg.005, Add.DEGtest == "4._D84")
D14.4 <- filter(IDDRC.deg.005, Add.DEGtest == "4._D14")
D84.3 <- filter(IDDRC.deg.005, Add.DEGtest == "3._D84")
D14.3 <- filter(IDDRC.deg.005, Add.DEGtest == "3._D14")
D84.2 <- filter(IDDRC.deg.005, Add.DEGtest == "2._D84")
D84.14 <- filter(IDDRC.deg.005, Add.DEGtest == "14._D84")
D14.14 <- filter(IDDRC.deg.005, Add.DEGtest == "14._D14")
D84.13 <- filter(IDDRC.deg.005, Add.DEGtest == "13._D84")
D84.12 <- filter(IDDRC.deg.005, Add.DEGtest == "12._D84")
D14.12 <- filter(IDDRC.deg.005, Add.DEGtest == "12._D14")
D84.11 <- filter(IDDRC.deg.005, Add.DEGtest == "11._D84")
D14.10 <- filter(IDDRC.deg.005, Add.DEGtest == "10._D14")
D84.10 <- filter(IDDRC.deg.005, Add.DEGtest == "10._D84")
D84.1 <- filter(IDDRC.deg.005, Add.DEGtest == "1._D84")
D84.0 <- filter(IDDRC.deg.005, Add.DEGtest == "0._D84")
D14.0 <- filter(IDDRC.deg.005, Add.DEGtest == "0._D14")
