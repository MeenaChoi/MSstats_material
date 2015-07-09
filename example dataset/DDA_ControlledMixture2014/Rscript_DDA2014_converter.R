#########################################
## Install MSstats (only need to install once)
#########################################
# 1. Please install all dependancy packages first.
install.packages(c("gplots","lme4","ggplot2","reshape","reshape2","data.table","Rcpp"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("limma","marray","preprocessCore","MSnbase"))

# 2. Install MSstats
install.packages(pkgs = "MSstats_3.0.6.tar.gz", repos = NULL, type ="source") 

# MSstats installation is complete ~~~

# Please exit R and re-enter R again. 
# Once re-enter R, please start with the below commands to exercute MSstats.



#########################################
## Run MSstats (start here if you already install MSstats)
#########################################
# load the library
library(MSstats)

# Help file
?MSstats

# 1. First, get protein ID information
proteinGroups<-read.table("proteinGroups.txt", sep="\t", header=TRUE)

# 2. Read in annotation including condition and biological replicates: annotation.csv
annot <- read.csv("annotation.csv", header=TRUE)


# 3. Read in MaxQuant file: evidence.txt
infile <- read.table("evidence.txt", sep="\t", header=TRUE)

# however, evidence.RData is available in MSstats material github due to limit of file size
infile<-load("evidence.RData")


# reformat for MSstats required input
# check options for converting format
?MaxQtoMSstatsFormat

msstats.raw<-MaxQtoMSstatsFormat(evidence=infile, annotation=annot, proteinGroups=proteinGroups)

head(msstats.raw)

save(msstats.raw, file="msstats.raw.RData")


