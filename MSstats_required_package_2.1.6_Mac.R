# step 1: install dependency packages
install.packages(c("gplots","lme4","ggplot2","reshape","data.table","Rcpp"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("limma","marray","preprocessCore","MSnbase"))

# step 2: install MSstats.daily - please chage the 'LocalPath' with your location.
install.packages(pkgs="LocalPath/MSstats.daily_2.1.6.tar.gz",repos=NULL,type="source")

# step 3: load the library of MSstats.daily
library(MSstats.daily)

# step 4: getting started
?MSstats.daily


install.packages(pkgs="/Users/meenachoi_iMac/Dropbox/MSstats_GitHub_private/MSstats.daily_2.3.7.tar.gz",repos=NULL,type="source")

install.packages(pkgs="/Users/Meena/Dropbox/MSstats_GitHub_private/MSstats.daily_3.0.3.tar.gz",repos=NULL,type="source")

install.packages("devtools")
library(devtools)

install_github("Cardinal", username="meena")
library(MSstats.daily)