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
raw<-read.csv("RawData.DDA.csv")
head(raw)

#=====================
# Function: dataProcess
# pre-processing data: quality control of MS runs
? dataProcess

quantData<-dataProcess(raw)
quantData$ProcessedData[1:5,]

# output QuantData
write.csv(quantData$ProcessedData,file="DDA_QuantData.csv")


#=====================
# Function: dataProcessPlots
# visualization 
? dataProcessPlots

dataProcessPlots(data=quantData,type="ProfilePlot",ylimUp=35,address="DDA_")
dataProcessPlots(data=quantData,type="QCPlot",ylimUp=35,address="DDA_")
dataProcessPlots(data=quantData,type="ConditionPlot",address="DDA_")


#=====================
# Function: modelBasedQCPlots
# visualization for model-based quality control in fitting model.
?modelBasedQCPlots

modelBasedQCPlots(data=quantData$ModelQC,type="ResidualPlots",featureName=FALSE,address="DDA_")
modelBasedQCPlots(data=quantData$ModelQC,type="QQPlots",feature.QQPlot="byFeature",address="DDA_")


#=====================
# Function: groupComparison
# generate testing results of protein inferences across concentrations
?groupComparison

levels(quantData$ProcessedData$GROUP_ORIGINAL)
comparison1<-matrix(c(1,-1,0,0,0,0),nrow=1)
comparison2<-matrix(c(1,0,-1,0,0,0),nrow=1)
comparison3<-matrix(c(1,0,0,-1,0,0),nrow=1)
comparison4<-matrix(c(1,0,0,0,-1,0),nrow=1)
comparison5<-matrix(c(1,0,0,0,0,-1),nrow=1)

comparison<-rbind(comparison1,comparison2,comparison3,comparison4,comparison5)
row.names(comparison)<-c("C1-C2","C1-C3","C1-C4","C1-C5","C1-C6")


resultMultiComparisons<-groupComparison(contrast.matrix=comparison,data=quantData)
resultMultiComparisons$ComparisonResult

write.csv(resultMultiComparisons$ComparisonResult, file="DDA_SignificanceTestingResult.csv")



#=====================
# Function: groupComparisonPlots
# visualization for testing results
?groupComparisonPlots


# Visualization 1: Volcano plot
# (1) default setup: FDR cutoff = 0.05; fold change cutoff = NA
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="VolcanoPlot", ylimUp=70, address="DDA_Ex1_")

# (2) Both FDR cutoff = 0.05; fold change cutoff = 1.5
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="VolcanoPlot", ylimUp=70,FCcutoff=1.5,ProteinName=FALSE, address="DDA_Ex2_")

# Visualization 2: Heatmap (required more than one comparisons)
# (1) default setup: FDR cutoff = 0.05; fold change cutoff = NA
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="Heatmap",address="DDA_Ex1_")

# (2) Both FDR cutoff = 0.05; fold change cutoff = 1.5
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="Heatmap",FCcutoff=1.5,address="DDA_Ex2_")

# Visualization 3: Comparison plot
groupComparisonPlots(data=resultMultiComparisons$ComparisonResult,type="ComparisonPlot",address="DDA_")


#=====================
# Function: designSampleSize
# calulate number of biological replicates per group for your next experiment
?designSampleSize

result.sample<-designSampleSize(data=resultMultiComparisons$fittedmodel,numSample=TRUE,desiredFC=c(1.5,2),FDR=0.05,power=0.9)
result.sample


#=====================
# Function: designSampleSizePlots
# visualization for sample size calculation
?designSampleSizePlots

pdf("DDA_SampleSizePlot.pdf")
designSampleSizePlots(data=result.sample)
dev.off()


#=====================
# Function: quantification
# model-based quantification
?quantification

# (1) group quantification
groupQuant<-quantification(data=quantData$ProcessedData,type="Group")
write.csv(groupQuant, file="DDA_GroupQuantification.csv")

# (2) sample quantification
sampleQuant<-quantification(data=quantData$ProcessedData,type="Sample")
write.csv(sampleQuant, file="DDA_SampleQuantification.csv")


