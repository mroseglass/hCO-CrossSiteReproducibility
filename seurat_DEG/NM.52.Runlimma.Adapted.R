options(stringsAsFactors=F)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(limma))

parser<-ArgumentParser()
parser$add_argument('--month', default='24',help='MRI phenotype [%(default)s]')
parser$add_argument('--res', default='0.8',help='MRI phenotype [%(default)s]')
args<-parser$parse_args()
print(args$month)
print(args$res)

res<-args$res

# Load the phenotype file
phenotype <- readRDS('/proj/steinlab/projects/IVIV_scRNA/InfoVerifiedID/RObjects/IVIV_IBIS_Culture_Technical_VerifiedIDs.rds')
colnames(phenotype) <- gsub(' ','_',colnames(phenotype))
taqman.col <- colnames(phenotype)[grepl('Taq',colnames(phenotype))]
phenotype %>% select(!all_of(taqman.col)) %>%
mutate(EBID = tolower(EBID)) %>%
unique() ->  phenotype

common.col <- c('EBID','DCCID','Sex','Collection_Location')

setwd('PseudoBulk/postGruffi')

count.dt <- readRDS(paste0('VST.Expr.PerSample_Day_Cluster.rm.lowExpGenes.',res,'.rds'))
# 18 clusters; day14 and day84
# day 14 : odd 
# day 84 : even

pheno <- paste0('v3.13_TotSA_V',args$month)
age <- paste0('mri_parameter_form,Candidate_Age_V',args$month)

target.col <- c(common.col,pheno,age)
phenotype.tmp <- phenotype %>% 
select(all_of(target.col)) %>% 
rename_at(vars(5), ~"phenotype") %>% 
rename_at(vars(6), ~"age") %>% 
filter(!is.na(phenotype) & !is.na(age)) %>% 
mutate(phenotype = phenotype/1000) %>% unique() # scale

results = NULL
for(i in 1:length(count.dt)){
 clusterID <- ( i - 1) %/%2 
 Day.target <- ifelse(i%%2==1,'Day14','Day84')
 message('cluster [', clusterID, '] Day [',Day.target,'] Phenotype [',pheno,']')
 if(is.null(count.dt[[i]])){
  message('no count data is available skip..')
  next;
 }

 # Make phenotype matrix
 pc.info <- count.dt[[i]]$pca %>% mutate(idx=1:nrow(count.dt[[i]]$pca))
   
 # Extract data exists in both seqData and phenotype data
 covariate_data <- pc.info %>% 
 select(batchID,WTK,cell_count,idx,cluster_sample_id) %>% 
 mutate(batchID = ifelse(batchID=="mich","mich2",batchID)) %>%  ### update mich to mich2
 full_join(phenotype.tmp,by=c('batchID'='EBID')) %>% na.omit()
 message('target data size (brain measurement + seqinfo):', nrow(covariate_data))
 rownames(covariate_data) <- covariate_data$cluster_sample_id

 # Subset gene expression matrix
 gene_matrix <- count.dt[[i]]$gexpr[,covariate_data$idx]
 
 # Make sure that the order is correct
 message('order is matched?:', all(rownames(covariate_data)  == colnames(gene_matrix)))

 covariate_data$WTK <- as.factor(covariate_data$WTK)
 covariate_data$DCCID <- as.factor(covariate_data$DCCID)

 # Exclude samples with < 10 cells
 excl.samples <- which(covariate_data$cell_count < 10)
 if(length(excl.samples) !=0){
  message(length(excl.samples),' sample(s) will be excluded due to low cell counts')
  covariate_data <- covariate_data[-excl.samples,]
  gene_matrix <- gene_matrix[,-excl.samples]
 }
 # Check if we have enough donors included
 N.uniq.donors <- covariate_data$DCCID %>% unique() %>% length()
 message(N.uniq.donors)
 if(N.uniq.donors < 10){
  next;
  message('no enough donor detected')
 }

 # =================
 # limma 
 # =================
 design = model.matrix ( ~ phenotype + WTK + age, covariate_data)
 dupcor <- duplicateCorrelation(gene_matrix,design,block=covariate_data$DCCID)
 fit1 <- lmFit(gene_matrix,design,block=covariate_data$DCCID,correlation=dupcor$consensus)
 fit1 <- eBayes(fit1)
 results.tmp <- topTable(fit1,coef=2,adjust.method="BH",number=Inf) %>%
 tibble::rownames_to_column("gene") %>%
 mutate(cluster=clusterID,Day = Day.target)
 results <- rbind(results,results.tmp)
}

saveRDS(results,file=paste0("Corr.res/limma.res.",pheno,".wAge.combined.",res,".rds"))
