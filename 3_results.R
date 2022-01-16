library(flowClust)
library(flowCore)
library(flowStats)
library(flowViz)
library(ggcyto) # only for autoplot which outputs histograms for each parameter
library(rgl) # package for plot3d function
library(ggplot2)

# clust.large: contains clustering parameters used for splitting posFF -flowClust-
# clustDF.18: data frame containing all cells with matching 'Class' number
# clustDF.19: negative events also
# cols: twenty contrasting colors -vector-
# dfList: list of data frames (conversion from flowFrames)
# expID: unique identifier of image data -string-
# fcsFF: created from wholeDataSet.fcs -flowFrage-
# fcsPath: file path of wholeDataSet.fcs (cell intensity data) +string+
# ffList: contains one flowFrame for each cluster -list-
# gate.blue/red/green: negative gates for singles colors -rectangleGate-
# gate.all: gates for all colors combined -rectangleGate-
# intensityColNames: the three column names containing intensity data -vector-
# negposFS: contains all events split into neg. and pos. flowFrames *flowSet*
# posFF: contains all color-positive events. for convenience purposes *flow frame*
# trans: transformation matrix for transFF ~transformList~
# transFF: contains transformed data from wholeDataSet.fcs -flowFrame-
# wd: working directory *string* 
# ...

# define file path
setwd(wd <- "C:/Users/stadt/Desktop/Thesis Data & Scripts")
expID <- "210121_SEM6_DRAQ5"
# string created from paste() with two arguments. file.path() function replaces folder path string.
fcsPath <- paste(file.path(".","3_Merged Data", expID, "All Events"), "wholeDataSet.fcs", sep="/")
# read in FCS file
fcsFlowFrame <- read.FCS(fcsPath)
intensityColNames <- colnames(fcsFlowFrame)[startsWith(colnames(fcsFlowFrame), "Mean")] # later needed for flowClust

# transform FCS file
trans <- estimateLogicle(fcsFlowFrame, intensityColNames)
transFF <- transform(fcsFlowFrame, trans)

# use one-dimensional gating for pos/negative populations. rangeGate does this automatically
gate.blue <- rangeGate(transFF, intensityColNames[1], plot = T, sd = .4, filterId = "blue", positive = F)
gate.green <- rangeGate(transFF, intensityColNames[2],  plot = T, sd = .4, filterId = "green", positive = F)
gate.red <- rangeGate(transFF, intensityColNames[3], plot = T, sd = .4, filterId = "red", positive = F)
gate.all <- gate.blue * gate.green * gate.red

# split flowFrame "transFF" into two flowFrames contained in flowSet negposFS
negposFS <- split(transFF, filter(transFF, f = gate.all), flowSet = T, filterId = "negative")
posFF <- negposFS[[2]]

# create data frame for negative population and label all cells in the column "class" with 1
negDF <- as.data.frame(exprs(negposFS[[1]]))
negDF$Class <- c(1)
negDF$Color <- cols[1]

# cluster flowframe based on transformed mean intensities -> what is a good amount of clusters?
clust.large <- flowClust(posFF, expName = "test", varNames = intensityColNames, K = 18)
plot(18, criterion(clust.large, "BIC"))
ffList <- split(posFF, clust.large) # creates a list containing K flowFrames

# converting flowFrames into data frames. these are contained in dfList
dfList <- list(negDF)
clustDF.19 <- negDF

#vector with different colors. are assigned in loop
cols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

for(i in 2:(length(ffList) + 1)){
  dfList[[i]] <- as.data.frame(exprs(ffList[[i]]))
  dfList[[i]]$Class <- c(i)
  dfList[[i]]$Color <- c(cols[i])
  clustDF.19 <- rbind(clustDF.19, dfList[[i]])
}


### ------------------------------------------------------------------------------------------

# sorting data frame by EventNo 
clustDF.19 <- clustDF.19[order(clustDF.19$EventNo), ]

# interactive 3D plot withcolors from vector
with(clustDF.19, plot3d(clustDF.19[1:nrow(clustDF.19), intensityColNames], col = clustDF.19$Color))
legend3d("topright", legend = paste('Cluster', c(1:19)), pch= 20, col = cols, cex=1)

# create list for readouts


# how many events (outliers) are lost during clustering
missinglines <- nrow(posFF)-nrow(clustDF.19)


