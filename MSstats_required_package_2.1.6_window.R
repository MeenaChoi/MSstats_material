# step 1: install dependency packages
install.packages(c("gplots","lme4","ggplot2","reshape","data.table","Rcpp"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("limma","marray","preprocessCore","MSnbase"))

# step 2: select 'Install package(s) from local zip files...' under 'Packages' in Menu bar. Then select 'MSstats.daily_2.1.6.zip' from your local directory


# step 3: load the library of MSstats.daily
library(MSstats.daily)

# step 4: getting started
?MSstats.daily