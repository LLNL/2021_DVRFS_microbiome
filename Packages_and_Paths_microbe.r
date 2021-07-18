#  Packages and Paths

# Load Packages
packages <- c("affyPLM", "dada2", "ShortRead", "ggplot2", "phyloseq", "vegan", "DESeq2", "dendextend", "tidyr", "viridis", "reshape2", "stringr", "ggthemes", "pander", "plyr", "ranacapa", "ade4", "FactoMineR", "factoextra", "ggrepel", "Heatplus", "dbstats", "Rcpp", "ape", "dplyr", "forcats", "colorspace", "ggsci", "microbiome", "data.table", "metagMisc", "ampvis", "iNEXT", "ggpubr", "RCM", "SpiecEasi", "biomformat", "qiime2R", "hablar", "MuMIn", "mltools", "decontam", "dlookr", "bestNormalize", "PerformanceAnalytics", "compositions", "robCompositions", "MuMIn", "MASS", "Hmisc", "knitr", "ALDEx2", "CoDaSeq", "zCompositions", "igraph", "grDevices", "car", "propr", "ecodist", "pheatmap", "sva", "bapred", "glmulti", "gridExtra", "readxl", "corrplot", "pairwiseAdonis", "Matrix", "igraph", "doMC")

### Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

### Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}

### Install packages not yet installed
devtools::install_github("gauravsk/ranacapa")
devtools::install_github("vmikk/metagMisc")
devtools::install_github("MadsAlbertsen/ampvis")
devtools::install_github("jbisanz/qiime2R")
devtools::install_github('ggloor/CoDaSeq/CoDaSeq')
devtools::install_github("zdk123/SpiecEasi", force=TRUE)

### Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Set Paths
path <- "<path to main dir>"
quality <- file.path(paste(path, "/01.quality", sep=""))
filter <- file.path(paste(path, "/02.filter_trim", sep=""))
error <- file.path(paste(path, "/03.error", sep=""))
dereplication <- file.path(paste(path, "/04.dereplication", sep=""))
merging <- file.path(paste(path, "/05.merging", sep=""))
input <- file.path(paste(path, "/input", sep=""))
path_phy <- file.path(paste(path, "/06.phyloseq", sep=""))
tax_file <- "<path to silva_nr_v132_train_set.fa>"
species_file <- "<path to silva_species_assignment_v132.fa>"

dir.create(quality, showWarnings = FALSE)
dir.create(filter, showWarnings = FALSE)
dir.create(error, showWarnings = FALSE)
dir.create(dereplication, showWarnings = FALSE)
dir.create(merging, showWarnings = FALSE)
dir.create(path_phy, showWarnings = FALSE)
