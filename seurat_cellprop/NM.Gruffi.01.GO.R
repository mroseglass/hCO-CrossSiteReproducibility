options(stringsAsFactors=F)

########################
### Prepare GO-terms and gene-sets
########################

ensembl <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
go1 <- "GO:0006096" # Glycolysis
go2 <- "GO:0034976" # ER-stress
go3 <- "GO:0042063" # Gliogenesis, negative filtering

save(ensembl,go1,go2,go3,file="/work/users/r/o/roseg/single-cell_reproducibility/data/Gruffi.stress.marker.rdata")
