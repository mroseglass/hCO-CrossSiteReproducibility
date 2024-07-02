#module load r/4.2.2
#inspired by https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/pseudobulk_DESeq2_scrnaseq.md and Nana Matoba
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(muscat,lib="../../../tools/R/4.2/"))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(dplyr))

message('Loading s.obj')
S.obj <- readRDS('/work/users/r/o/roseg/single-cell_reproducibility/results/May12023_IDDRC_PGP1_UNCR3_reference_integrated_clustered_seurat_object.rds')

#for(res in c(0.2,0.4,0.6,0.8,1)){
#for(res in c(0.4,0.6,0.8,1)){
res =0.6
 message('Preparing metadata/reset')
 counts <- S.obj@assays$RNA@counts
 metadata <- S.obj@meta.data
 message(res)
 metadata$cluster_id <- factor(S.obj[[paste0('integrated_snn_res.',res)]][[1]])
 metadata %>% rename('sample_id'='Sample') -> metadata
 metadata <- mutate(metadata, SiteDay = paste0(Site,Day))

 #Site and day alrady columns in metadata
 #metadata %>% mutate(DCCID = sub('DCCID_','',DonorID)) -> metadata
 #metadata$group_id <- paste0(metadata$DCCID,'_',metadata$Day)

 message('Creating sce obj')
 # Create single cell experiment object
 sce <- SingleCellExperiment(assays = list(counts = counts),
                           colData = metadata)

 cluster_names <- levels(colData(sce)$cluster_id)
 message('# of clusters: ',length(cluster_names))
 sample_names <- levels(as.factor(colData(sce)$sample_id))
 message('# of samples: ',length(sample_names))

 # remove Zelda
 #pgp1_sph.cells <- sample_names[grepl('pgp1|sph',sample_names)]
 #sce.filtered <- subset(sce, ,!sample_id %in% pgp1_sph.cells)
 zelda.cells <- sample_names[grepl('Zelda',sample_names)]
 sce.filtered <- subset(sce, ,!sample_id %in% zelda.cells)

 message('remove Zelda')
 cluster_names <- levels(colData(sce.filtered)$cluster_id)
 message('# of clusters: ',length(cluster_names))
 sample_names <- levels(as.factor(colData(sce.filtered)$sample_id))
 message('# of samples: ',length(sample_names))

 pb <- aggregateData(sce.filtered,
     assay = "counts", fun = "sum",
     by = c("cluster_id", "SiteDay"))
#change above like to c("cluster_id","sample_id") to pseudobulk by expeirment

 # Get the assay data from object pb as a list
 filtered_assay <- cbind.data.frame(lapply(assays(pb), function(x) {
  # Remove columns with all 0s
  cols_to_keep <- colSums(x != 0) > 0
  x[, cols_to_keep]
 }))

 # Create a sparse matrix from the filtered data,
 num_mat <- as.matrix(filtered_assay)
 row_names <- rownames(filtered_assay)
 col_names <- colnames(filtered_assay)
 # Convert the filtered assay object to a sparse matrix while preserving row and column names
 sparse_mat <- sparseMatrix(i = row(num_mat), j = col(num_mat), x = as.numeric(num_mat),dimnames = list(row_names, col_names))
 aggr_counts <- drop0(sparse_mat)
 #print(colData(pb))
 rm(list=c('sparse_mat','filtered_assay','pb','num_mat','row_names','col_names'))

 #pseudobulk for each experiment
 print(tstrsplit(colnames(aggr_counts), "\\.") %>% str())
 # check cluster
 print(tstrsplit(colnames(aggr_counts), "\\.")[[1]] %>% unique())
 # check samples
 print(tstrsplit(colnames(aggr_counts), "\\.")[[2]] %>% unique())

 counts_ls <- list()

 for (i in 1:length(cluster_names)) {
  ## Extract indexes of columns in the global matrix that match a given cluster
  column_idx <- which(tstrsplit(colnames(aggr_counts), "\\.")[[1]] == cluster_names[i])
  ## Store corresponding sub-matrix as one element of a list
  counts_ls[[i]] <- aggr_counts[, column_idx]
  names(counts_ls)[i] <- cluster_names[i]
 }
 #str(counts_ls)

 message('check metadta')
 metadata <- metadata %>%
 select(sample_id, WTK_ID,cluster_id,Site,Day,Rep)  %>%
 filter(!grepl('Zelda',Site)) #%>% unique()

 cell.count.t <- table(metadata$sample_id,metadata$cluster_id)
 metadata <- metadata %>% unique()

 # Initialize an empty list to store metadata for each cluster
 metadata_ls <- list()

 # Loop through each cluster in the counts_ls list
 for (i in 1:length(counts_ls)) {
  df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))

  df$cluster_id <- tstrsplit(df$cluster_sample_id, "\\.")[[1]]
  df$sample_id  <- tstrsplit(df$cluster_sample_id, "\\.")[[2]]

  # Retrieve cell count information for this cluster from global cell count table
  cluster_name <- unique(df$cluster_id)
  cell_counts <- cell.count.t[, match(cluster_name, colnames(cell.count.t))]

  # Remove samples with zero cell contributing to the cluster
  cell_counts <- cell_counts[cell_counts > 0]

  # Match order of cell_counts and sample_ids
  sample_order <- match(df$sample_id, names(cell_counts))
  cell_counts <- cell_counts[sample_order]

  df$cell_count <- cell_counts

  # Join data frame (capturing metadata specific to cluster) to generic metadata
  df <- left_join(df, metadata, by = intersect(names(df), names(metadata)))
  rownames(df) <- df$cluster_sample_id

  # Store complete metadata for cluster i in list
  metadata_ls[[i]] <- df
  names(metadata_ls)[i] <- cluster_name
 }
 save(metadata_ls,counts_ls,file=paste0("/work/users/r/o/roseg/single-cell_reproducibility/results/seurat/aggr.sum.cluster_Site_Day.res",res,".rdata"))
#}

