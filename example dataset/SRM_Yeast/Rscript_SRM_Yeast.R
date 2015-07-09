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
raw<-read.csv("RawData.yeast.csv")
head(raw)

#=====================
# Function: dataProcess
# pre-processing data: quality control of MS runs
? dataProcess
quantData<-dataProcess(raw)
quantData[1:5,]

# output QuantData
write.csv(quantData,file="Yeast_QuantData.csv")


#=====================
# Function: dataProcessPlots
# visualization 
? dataProcessPlots
dataProcessPlots(data=quantData,type="ProfilePlot",address="Yeast_")
dataProcessPlots(data=quantData,type="QCPlot",address="Yeast_")
dataProcessPlots(data=quantData,type="ConditionPlot",address="Yeast_")


#=====================
# Function: groupComparison
# generate testing results of protein inferences across concentrations
?groupComparison
levels(quantData$GROUP_ORIGINAL)
comparison<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
row.names(comparison)<-"T7-T1"

resultOneComparison<-groupComparison(contrast.matrix=comparison,data=quantData)
resultOneComparison$ComparisonResult


# testing with more than one comparisons
comparison1<-matrix(c(-1,0,1,0,0,0,0,0,0,0),nrow=1)
comparison2<-matrix(c(-1,0,0,0,1,0,0,0,0,0),nrow=1)
comparison3<-matrix(c(-1,0,0,0,0,0,1,0,0,0),nrow=1)
comparison<-rbind(comparison1,comparison2,comparison3)
row.names(comparison)<-c("T3-T1","T5-T1","T7-T1")

resultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=quantData)
resultMultiComparisons$ComparisonResult

write.csv(resultMultiComparisons$ComparisonResult, file="Yeast_SignificanceTestingResult.csv")



#=====================
# Function: modelBasedQCPlots
# visualization for model-based quality control in fitting model.
?modelBasedQCPlots

modelBasedQCPlots(data=resultMultiComparisons$ModelQC,type="ResidualPlots",address="OvarianCancer_")

modelBasedQCPlots(data=resultMultiComparisons$ModelQC,type="QQPlots",address="OvarianCancer_")



#=====================
# Function: groupComparisonPlots
# visualization for testing results
?groupComparisonPlots


# Visualization 1: Volcano plot
# (1) default setup: FDR cutoff = 0.05; fold change cutoff = NA
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="VolcanoPlot", ylimUp=70, address="Yeast_Ex1_")

# (2) Both FDR cutoff = 0.05; fold change cutoff = 1.5
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="VolcanoPlot", ylimUp=70,FCcutoff=1.5,ProteinName=FALSE, address="Yeast_Ex2_")

# Visualization 2: Heatmap (required more than one comparisons)
# (1) default setup: FDR cutoff = 0.05; fold change cutoff = NA
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="Heatmap",address="Yeast_Ex1_")

# (2) Both FDR cutoff = 0.05; fold change cutoff = 1.5
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="Heatmap",FCcutoff=1.5,address="Yeast_Ex2_")

# Visualization 3: Comparison plot
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="ComparisonPlot",address="Yeast_Ex1_")


#=====================
# Function: designSampleSize
# calulate number of biological replicates per group for your next experiment
?designSampleSize

result.sample<-designSampleSize(data=quantData,numSample=TRUE,numPep=2,numTran=2,desiredFC=c(1.05,1.3),FDR=0.01,power=0.9)
result.sample


#=====================
# Function: designSampleSizePlots
# visualization for sample size calculation
?designSampleSizePlots

pdf("Yeast_SampleSizePlot.pdf")
designSampleSizePlots(data=result.sample)
dev.off()


#=====================
# Function: quantification
# model-based quantification
?quantification

# (1) group quantification
groupQuant<-quantification(data=quantData,type="Group")
write.csv(groupQuant, file="Yeast_GroupQuantification.csv")

# (2) sample quantification
sampleQuant<-quantification(data=quantData,type="Sample")
write.csv(sampleQuant, file="Yeast_SampleQuantification.csv")


