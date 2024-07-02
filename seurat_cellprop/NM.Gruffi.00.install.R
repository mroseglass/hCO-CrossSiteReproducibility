# removed specification to NANA's personal library
# depend on your enviroment you may need to install extra packages and do settings
#

install.packages('pdist')
install.packages('raster')
#install.packages('rgl',lib="/work/users/n/a/nanam/tools/R/4.2/",INSTALL_opts="--no-test-load") # needed to re-install

#install.packages('BiocManager')

# Install custom dependencies
#install.packages('devtools')
devtools::install_github(repo = "vertesy/Stringendo", upgrade = F,)
devtools::install_github(repo = "vertesy/CodeAndRoll2", upgrade = F)
devtools::install_github(repo = "vertesy/ReadWriter", upgrade = F)
devtools::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
devtools::install_github(repo = "vertesy/Markdownreports", upgrade = F)
## terminal ->  module load gcc/11.2.0
devtools::install_github(repo = "vertesy/ggExpress", upgrade = F)
####set LD path on terminal before installing Seurat.utils
##export LD_LIBRARY_PATH=/nas/longleaf/rhel8/apps/gcc/11.2.0/lib64:/nas/longleaf/rhel8/apps/gcc/6.3.0/lib:/nas/longleaf/rhel8/apps/gcc/6.3.0/lib64:/nas/longleaf/rhel8/apps/r/4.1.3/lib64/R/lib

devtools::install_github(repo = "vertesy/Seurat.utils", upgrade = F)

# Install gruffi
devtools::install_github(repo = "jn-goe/gruffi", upgrade = F)

#Installing package into ‘/nas/longleaf/home/nanam/R/x86_64-pc-linux-gnu-library/4.1’
#(as ‘lib’ is unspecified)
#ERROR: dependencies ‘Stringendo’, ‘CodeAndRoll2’, ‘ggExpress’, ‘MarkdownHelpers’, ‘MarkdownReports’, ‘Seurat.utils’ are not available for package ‘gruffi’
