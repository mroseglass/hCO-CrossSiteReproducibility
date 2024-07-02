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
deg.here <- "/work/users/r/o/roseg/IDDRC/IDDRCDatabase/CellTypeProportions/D14genetoD84Prop/"
plot.here <- "/work/users/r/o/roseg/single-cell_reproducibility/doc/seurat/pdf/Figures/DESeq2"

CellTypeID <- read_csv(paste0(db.here, "ClustertoCellType.csv"))
CellTypeID$Cluster <- as.character(CellTypeID$Cluster)
# GLOBALS ##########################################################################################
Site.color <-c("#12783D", "#882155", "#4040C6")
Replicate.shapes <- c(0,1,2,3,4,15,16,17)

#####################################################################################################
#Example opening of one file
List.deg.expt <- c(list.files(deg.here))
i=1
deseq2.result <- read_csv(paste0(deg.here,List.deg.expt[i]))
#drop NA values
deseq2.result <- deseq2.result%>% drop_na(`padj`)
#add name of test to table
d14clustername <- gsub ("_clusterD14.*","",List.deg.expt[i])
d84clusternumber <- gsub (".*D14","",List.deg.expt[i])
d84clusternumber <- gsub ("\\.csv","",d84clusternumber)
DEGtest.d14cluster <- rep(c(d14clustername),each=length(deseq2.result$log2FoldChange))
IDDRC.deg <- cbind(deseq2.result,DEGtest.d14cluster)
DEGtest.d84cluster <- rep(c(d84clusternumber),each=length(deseq2.result$log2FoldChange))
IDDRC.deg <- cbind(IDDRC.deg,DEGtest.d84cluster)
DEGtest <- rep(c(List.deg.expt[i]),each=length(deseq2.result$log2FoldChange))
IDDRC.deg <- cbind(IDDRC.deg,DEGtest)


#Add rest of files to same object
for (i in 2:length(List.deg.expt)){
  deseq2.result <- read_csv(paste0(deg.here,List.deg.expt[i]))
  #drop NA values
  deseq2.result <- deseq2.result%>% drop_na(`padj`)
  #add name of test to table
  d14clustername <- gsub ("_clusterD14.*","",List.deg.expt[i])
  d84clusternumber <- gsub (".*D14","",List.deg.expt[i])
  d84clusternumber <- gsub ("\\.csv","",d84clusternumber)
  DEGtest.d14cluster <- rep(c(d14clustername),each=length(deseq2.result$log2FoldChange))
  CorrectFile2 <- cbind(deseq2.result,DEGtest.d14cluster)
  DEGtest.d84cluster <- rep(c(d84clusternumber),each=length(deseq2.result$log2FoldChange))
  CorrectFile2 <- cbind(CorrectFile2,DEGtest.d84cluster)
  DEGtest <- rep(c(List.deg.expt[i]),each=length(deseq2.result$log2FoldChange))
  CorrectFile2 <- cbind(CorrectFile2,DEGtest)
  IDDRC.deg <- rbind(CorrectFile2,IDDRC.deg)
}

# total test 2090315
# just padj significant
IDDRC.deg.005 <- filter(IDDRC.deg, padj < 0.05)
#just 298 genes
#write_csv(IDDRC.deg.005, paste0(db.here, "D14pseudbulked_genes_to_D84CellTypeProportion.csv"))
IDDRC.deg.005 <- read_csv(paste0(db.here, "D14pseudbulked_genes_to_D84CellTypeProportion.csv"))
IDDRC.deg.005 <- mutate(IDDRC.deg.005, DEGtest.d84cluster = gsub("c","",IDDRC.deg.005$DEGtest.d84cluster))
CellTypeID <- dplyr::rename(CellTypeID, DEGtest.d84cluster = Cluster)
IDDRC.deg.005 <- left_join(IDDRC.deg.005, CellTypeID, by = "DEGtest.d84cluster")

ggplot(data = IDDRC.deg.005, aes(y=DEGtest.d14cluster)) +ggplot(data = IDDRC.deg.005, aes(y=DEGtest.d14cluster)) +
  geom_bar() +
  theme_bw() +
  ylab("Day 14 Pseudobulked Cells") +
  xlab("Count Genes with FDR Significant Correlations to a D84 Cell Type") +
  ggtitle("Number D84 cell type gene correlations")

ggplot(data = IDDRC.deg.005, aes(y=CellType)) +
  geom_bar() +
  theme_bw() +
  ylab("Day 84 Cell Type") +
  xlab("Count Genes with FDR Significant Correlations to a D14 Cell Type") +
  ggtitle("Number D84 cell type gene correlations")

gostres <- gost(query = IDDRC.deg.005$gene,
                organism = "hsapiens", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = FALSE,
                user_threshold = 0.05, correction_method = "g_SCS",
                domain_scope = "annotated", custom_bg = NULL,
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
DEG.ann <- gostres$result
DEG.ann.all <- filter(DEG.ann, source == "GO:BP")

#just pull out progenitors cell types
#rg is just 5 gene 2 ENGS genes and other correlate to RG and NPCs
IDDRC.deg.005.rg <- filter(IDDRC.deg.005, DEGtest.d14cluster == "Radial.glia")
IDDRC.deg.005.nesc <- filter(IDDRC.deg.005, DEGtest.d14cluster == "Neuroepithelial.stem.cells")
IDDRC.deg.005.npc.s <- filter(IDDRC.deg.005, DEGtest.d14cluster == "Dividing.neural.progenitor.cells,.S")
IDDRC.deg.005.npc.g <- filter(IDDRC.deg.005, DEGtest.d14cluster == "Dividing.neural.progenitor.cells,.G2")
IDDRC.deg.005.prog <- rbind(IDDRC.deg.005.nesc,IDDRC.deg.005.npc.g,IDDRC.deg.005.npc.s,IDDRC.deg.005.rg)

IDDRC.deg.005.prog$DEGtest.d14cluster <- factor(IDDRC.deg.005.prog$DEGtest.d14cluster,
                       levels=c("Neuroepithelial.stem.cells","Dividing.neural.progenitor.cells,.S",
                                "Radial.glia","Dividing.neural.progenitor.cells,.G2"))

ggplot(data = IDDRC.deg.005.prog, aes(y=DEGtest.d14cluster)) +
    geom_bar() +
    theme_bw() +
    ylab("Day 14 Pseudobulked Cells") +
    xlab("Count Genes with FDR Significant Correlations to a D84 Cell Type") +
    ggtitle("Number D84 cell type gene correlations")

gostres <- gost(query = IDDRC.deg.005.prog$gene,
                organism = "hsapiens", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = FALSE,
                user_threshold = 0.05, correction_method = "g_SCS",
                domain_scope = "annotated", custom_bg = NULL,
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
DEG.ann <- gostres$result
DEG.ann.all <- filter(DEG.ann, source == "GO:BP")

#write_csv(IDDRC.deg.005.prog, paste0(db.here, "D14pseudbulked_genes_to_D84CellTypeProportion_JustProg.csv"))
IDDRC.deg.005.prog <- read_csv(paste0(db.here, "D14pseudbulked_genes_to_D84CellTypeProportion_JustProg.csv"))
#####################################################################################
#check an example to see if the DEG make sense
gexpr <- readRDS("/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/VST.Expr.PerSample_Day_Cluster.rm.lowExpGenes.0.6.rds")

#for NES
dc.exp <- as_tibble(gexpr[[1]]$gexpr, rownames = NA)
dc.exp <- rownames_to_column(dc.exp, var = "gene") %>% as_tibble()
IDDRC.deg.005.prog <- dplyr::rename(IDDRC.deg.005.prog, gene = `...1`)
dc.exp <- inner_join(IDDRC.deg.005.prog, dc.exp, by = "gene")

D84prop <- read_csv(paste0(db.here,"CellTypeProportions/Day84CellTypePropwithGRUFFI.csv"))

dc.exp.p <- filter(dc.exp, gene == "ATP11A")
plot <- dplyr::select(dc.exp.p, !c(CellType,log2FoldChange,lfcSE,stat,pvalue,DEGtest,padj,baseMean,DEGtest.d14cluster,DEGtest.d84cluster))
#plot all DESeq2 anova sig gene  by site
plot.adjp <- dplyr::select(dc.exp.p,c(gene,padj))
plot <- pivot_longer(plot, cols = !gene, names_to = "Experiment", values_to = "VSTexpression")
plot$Experiment <- gsub("0.","",plot$Experiment)
plot <- mutate(plot, Site = gsub("_.*","",Experiment))
plot <- mutate(plot, Rep = gsub(".*_","",Experiment))
plot.r <- full_join(plot.adjp,plot, by = "gene")
c8Prop <- dplyr::select(D84prop, c(Experiment.scRNAseq,"Pan neuron/Cajal - Retzius"))
c8Prop <- mutate(c8Prop, Experiment = gsub("D84","D14", c8Prop$Experiment.scRNAseq))
plot.rr <- inner_join(plot.r, c8Prop, by = "Experiment")

VSTcor <- cor.test(plot.rr$VSTexpression, plot.rr$`Pan neuron/Cajal - Retzius`, method="pearson")

cell.hm <- ggplot(plot.rr, aes(x = VSTexpression, y=`Pan neuron/Cajal - Retzius`, color=Site)) +
    scale_color_manual(values = Site.color) +
    scale_shape_manual(values=Replicate.shapes)+
    theme_classic()+
    geom_smooth(method = lm, color = "black") +
        #geom_point(aes(shape = Rep),size = 2) +
    geom_point(size = 2) +labs(title= "Day 14 ATP11A Expression to Pan CN/CR Proportion Day 84",
                               subtitle = "padj = 0.0457807440",
                               caption = paste0("r =",VSTcor$estimate))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
cell.hm

#for top hit
dc.exp <- as_tibble(gexpr[[19]]$gexpr, rownames = NA) #unspecified neurons
dc.exp <- rownames_to_column(dc.exp, var = "gene") %>% as_tibble()
IDDRC.deg.005 <- dplyr::rename(IDDRC.deg.005, gene = `...1`)
dc.exp <- inner_join(IDDRC.deg.005, dc.exp, by = "gene")

D84prop <- read_csv(paste0(db.here,"CellTypeProportions/Day84CellTypePropwithGRUFFI.csv"))

dc.exp.p <- filter(dc.exp, gene == "SH3GL2")
dc.exp.p <- filter(dc.exp.p, DEGtest.d14cluster == "Unspecified.neuron")
plot <- dplyr::select(dc.exp.p, !c(CellType,log2FoldChange,lfcSE,stat,pvalue,DEGtest,padj,baseMean,DEGtest.d14cluster,DEGtest.d84cluster))
#plot all DESeq2 anova sig gene  by site
plot.adjp <- dplyr::select(dc.exp.p,c(gene,padj))

plot <- pivot_longer(plot, cols = !gene, names_to = "Experiment", values_to = "VSTexpression")
plot$Experiment <- gsub("9.","",plot$Experiment)
plot <- mutate(plot, Site = gsub("_.*","",Experiment))
plot <- mutate(plot, Rep = gsub(".*_","",Experiment))
plot.r <- full_join(plot.adjp,plot, by = "gene")
c8Prop <- dplyr::select(D84prop, c(Experiment.scRNAseq,"Pan neuron/Cajal - Retzius"))
c8Prop <- mutate(c8Prop, Experiment = gsub("D84","D14", c8Prop$Experiment.scRNAseq))
plot.rr <- inner_join(plot.r, c8Prop, by = "Experiment")
VSTcor <- cor.test(plot.rr$VSTexpression, plot.rr$`Pan neuron/Cajal - Retzius`, method="pearson")

cell.hm <- ggplot(plot.rr, aes(x = VSTexpression, y=`Pan neuron/Cajal - Retzius`, color=Site)) +
    scale_color_manual(values = Site.color) +
    scale_shape_manual(values=Replicate.shapes)+
    theme_classic()+
    geom_smooth(method = lm, color = "black") +
    #geom_point(aes(shape = Rep),size = 2) +
    geom_point(size = 2) + labs(title="Day 14 SH3GL2 Expression to Pan CN/CR Proportion Day 84",
                                subtitle = "padj = 1.0163e-06",
                                caption = paste0("r =",VSTcor$estimate))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
cell.hm




forlabels <- distinct(select(plot.r, c(gene,padj)))
gene.list <- c(forlabels$gene)
padj.list <-c(forlabels$padj)


