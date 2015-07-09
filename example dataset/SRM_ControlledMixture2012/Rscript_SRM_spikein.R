#########################################
## Part 1: Install MSstats (only need to install once)
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
## Part 2: prepare Input data for MSstats
#########################################
# Three example datasets in the manuscript are generated from MultiQuant
# The RawData is ready for MSstats (which already add column Condition, BioReplicate, Run, Intensity)


#########################################
## Part 3: Run MSstats (start here if you already install MSstats)
#########################################
# load the library
library(MSstats)

# Help file
?MSstats

# Input data
raw<-read.csv(file="RawData.SpikeIn.csv")
head(raw)

#=====================
# Function: dataProcess
# pre-processing data: quality control of MS runs
? dataProcess
quantData<-dataProcess(raw)
quantData$ProcessedData[1:5,]

# output QuantData
write.csv(quantData$ProcessedData,file="SpikeIn_QuantData.csv")

#=====================
# Function: dataProcessPlots
# visualization 
? dataProcessPlots
dataProcessPlots(data=quantData,type="ProfilePlot",address="SpikeIn_")
dataProcessPlots(data=quantData,type="QCPlot",address="SpikeIn_")
dataProcessPlots(data=quantData,type="ConditionPlot",address="SpikeIn_")


#=====================
# Function: modelBasedQCPlots
# visualization for model-based quality control in fitting model.
?modelBasedQCPlots

modelBasedQCPlots(data=quantData$ModelQC,type="ResidualPlots",address="SpikeIn_")
modelBasedQCPlots(data=quantData$ModelQC,type="QQPlots",feature.QQPlot="byFeature",address="SpikeIn_")


#=====================
# Function: groupComparison
# generate testing results of protein inferences across concentrations
?groupComparison
levels(quantData$ProcessedData$GROUP_ORIGINAL)
comparison<-matrix(c(0,0,-1,0,0,1),nrow=1)
row.names(comparison)<-"M6-M3"

resultOneComparison<-groupComparison(contrast.matrix=comparison,data=quantData)
resultOneComparison$ComparisonResult


# testing with more than one comparisons
comparison1<-matrix(c(0,0,-1,0,0,1),nrow=1)
comparison2<-matrix(c(0,0,0,-1,0,1),nrow=1)
comparison3<-matrix(c(0,0,0,0,-1,1),nrow=1)
comparison<-rbind(comparison1,comparison2,comparison3)
row.names(comparison)<-c("M6-M3","M6-M4","M6-M5")


resultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=quantData)
resultMultiComparisons$ComparisonResult

write.csv(resultMultiComparisons$ComparisonResult, file="SpikeIn_SignificanceTestingResult.csv")


#=====================
# Function: groupComparisonPlots
# visualization for testing results
?groupComparisonPlots

# Visualization 1: Volcano plot
# (1) default setup: FDR cutoff = 0.05; fold change cutoff = NA
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="VolcanoPlot",ProteinName=FALSE,address="SpikeIn_Ex1_")

# (2) Both FDR cutoff = 0.05; fold change cutoff = 2
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="VolcanoPlot",FCcutoff=2,ProteinName=FALSE,address="SpikeIn_Ex2_")

# Visualization 2: Heatmap (required more than one comparisons)
# (1) default setup: FDR cutoff = 0.05; fold change cutoff = NA
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="Heatmap", address="SpikeIn_Ex1_")

# (2) Both FDR cutoff = 0.05; fold change cutoff = 2
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="Heatmap",FCcutoff=2, address="SpikeIn_Ex2_")

# Visualization 3: Comparison plot
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="ComparisonPlot",address="SpikeIn_")



#=====================
# Function: designSampleSize
# calulate number of biological replicates per group for your next experiment
?designSampleSize

result.sample<-designSampleSize(data=resultMultiComparisons$fittedmodel,numSample=TRUE,desiredFC=c(1.2,1.5),FDR=0.05,power=0.8)
result.sample


#=====================
# Function: designSampleSizePlots
# visualization for sample size calculation
?designSampleSizePlots

pdf("SpikeIn_SampleSizePlot.pdf")
designSampleSizePlots(data=result.sample)
dev.off()


#=====================
# Function: quantification
# model-based quantification
?quantification

# (1) group quantification
groupQuant<-quantification(data=quantData$ProcessedData,type="Group")
write.csv(groupQuant, file="SpikeIn_GroupQuantification.csv")

# (2) sample quantification
sampleQuant<-quantification(data=quantData$ProcessedData,type="Sample")
write.csv(sampleQuant, file="SpikeIn_SampleQuantification.csv")


