#########################################
## Part 1: Install MSstats (only need to install once)
#########################################
# 1. Please install all dependancy packages first.
install.packages(c("gplots","lme4","ggplot2","reshape","data.table","Rcpp"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("limma","marray","preprocessCore","MSnbase"))

# 2. Install MSstats
install.packages(pkgs = "MSstats.daily_2.1.6.tar.gz", repos = NULL, type ="source") 

# MSstats installation is complete ~~~

# Please exit R and re-enter R again. 
# Once re-enter R, please start with the below commands to exercute MSstats.


#########################################
## Part 2: prepare Input data for MSstats
#########################################
# Three example datasets in the manuscript are generated from MultiQuant
# The RawData is ready for MSstats (which already add column Condition, BioReplicate, Run, Intensity)


#########################################
## Part 3: Run MSstats (start here if you already install MSstats)
#########################################
# load the library
library(MSstats.daily)

# Help file
?MSstats.daily

# Input data
raw<-read.csv("RawData.ovarian.csv")
head(raw)

#=====================
# Function: dataProcess
# pre-processing data: quality control of MS runs
? dataProcess
quantData<-dataProcess(raw)
quantData[1:5,]

# output QuantData
write.csv(quantData,file="OvarianCancer_QuantData.csv")


#=====================
# Function: dataProcessPlots
# visualization 
? dataProcessPlots
dataProcessPlots(data=quantData,type="ProfilePlot",address="OvarianCancer_")
dataProcessPlots(data=quantData,type="QCPlot",address="OvarianCancer_")
dataProcessPlots(data=quantData,type="ConditionPlot",address="OvarianCancer_")


#=====================
# Function: groupComparison
# generate testing results of protein inferences across concentrations
?groupComparison
levels(quantData$GROUP_ORIGINAL)
comparison<-matrix(c(-1,1),nrow=1)
row.names(comparison)<-"tumor-control"

resultOneComparison<-groupComparison(contrast.matrix=comparison,data=quantData)
resultOneComparison$ComparisonResult

write.csv(resultOneComparison$ComparisonResult, file="OvarianCancer_SignificanceTestingResult.csv")


#=====================
# Function: modelBasedQCPlots
# visualization for model-based quality control in fitting model.
?modelBasedQCPlots

modelBasedQCPlots(data=resultOneComparison$ModelQC,type="ResidualPlots",address="OvarianCancer_")

modelBasedQCPlots(data=resultOneComparison$ModelQC,type="QQPlots",address="OvarianCancer_")


#=====================
# Function: groupComparisonPlots
# visualization for testing results
?groupComparisonPlots

# Visualization 1: Volcano plot
# (1) default setup: FDR cutoff = 0.05; fold change cutoff = NA
groupComparisonPlots(data=resultOneComparison$ComparisonResult,type="VolcanoPlot",address="OvarianCancer_Ex1_")

# (2) Both FDR cutoff = 0.05; fold change cutoff = 1.5
groupComparisonPlots(data=resultOneComparison$ComparisonResult,type="VolcanoPlot",FCcutoff=1.5,ProteinName=FALSE,address="OvarianCancer_Ex2_")

# Visualization 2: Heatmap (required more than one comparisons)
# No heatmap is generate in Ovarian cancer study

# Visualization 3: Comparison plot
groupComparisonPlots(data=resultOneComparison$ComparisonResult,type="ComparisonPlot",address="OvarianCancer_Ex1_")


#=====================
# Function: designSampleSize
# calulate number of biological replicates per group for your next experiment
?designSampleSize

result.sample<-designSampleSize(data=quantData,numSample=TRUE,numPep=3,numTran=2,desiredFC=c(1.05,1.3),FDR=0.05,power=0.8)

result.sample


#=====================
# Function: designSampleSizePlots
# visualization for sample size calculation
?designSampleSizePlots

pdf("OvarianCancer_SampleSizePlot.pdf")
designSampleSizePlots(data=result.sample)
dev.off()


#=====================
# Function: quantification
# model-based quantification
?quantification

# (1) group quantification
groupQuant<-quantification(data=quantData,type="Group")
write.csv(groupQuant, file="OvarianCancer_GroupQuantification.csv")

# (2) sample quantification
sampleQuant<-quantification(data=quantData,type="Sample")
write.csv(sampleQuant, file="OvarianCancer_SampleQuantification.csv")


