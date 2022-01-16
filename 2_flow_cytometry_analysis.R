rm(list = ls()) # clear environment

library(stringr)
library(scatterplot3d)
library(rgl)
library('flowCore')
library('Biobase')
library('data.table')

setwd(wd <- "C:/Users/stadt/Desktop/Thesis Data & Scripts")

# expID: STRING with unique identifier of image data
# count: NUMBER of rows in wholeDataSet
# fcsSource: STRING with creation path of final .csv table (used for FCS conversion as well)
# fcsData: TABLE based on wholeDataSet data frame
# fcsData.ff: FLOW FRAME created from fcsData table
# metadata: metadata for created flow frame
# nameField: STRING with file name of merged .csv for 
# numFields: NUMBER of fields/imaged positions in dataset
# pathEnd: STRING with path of all event.csv
# pathOut: STRING with path where directory for results is created
# pathTables: STRING with path where raw .csv are located
# pathTbOutput: STRING with path for csv. file which stores only information about one field
# vecMrgNames: VECTOR with all file directories in pathOut
# vecTblNames: VECTOR with all file directories in pathTables
# wholeDataSet: DATA FRAME containing the final dataset

expID <- "210121_SEM6_DRAQ5"

pathOut <- file.path(".","3_Merged Data", expID)
dir.create(pathOut)

pathTables <- file.path(".","2_Data Extraction", expID) # path in which all unmerged tables are contained
vecTblNames <- list.files(path = pathTables, pattern=".csv", full.names = TRUE) # full.names=true gives whole path for file
numFields <- length(vecTblNames)/3 # fieldNum in old version

for(i in 0:(numFields-1)){
  
  df1 <- read.csv(vecTblNames[1+(i*3)], stringsAsFactors = FALSE)
  df2 <- read.csv(vecTblNames[2+(i*3)], stringsAsFactors = FALSE)
  df3 <- read.csv(vecTblNames[3+(i*3)], stringsAsFactors = FALSE)
  
  df1$Label <- NULL
  
  df2$Label <- NULL
  df2$Area <- NULL
  
  df3$Label <- NULL
  df3$Area <- NULL
  
  colnames(df1) <- c("ObjectNr", "Area", "MeanInten_b", "StdDev_b")
  colnames(df2) <- c("ObjectNr", "MeanInten_g", "StdDev_g")   
  colnames(df3) <- c("ObjectNr", "MeanInten_r", "StdDev_r") 
  
  df_interm <- merge(df1, df2, by="ObjectNr")
  df_final <- merge(df_interm, df3, by="ObjectNr")
  df_final$ObjectNr <- NULL
  
  nameField <- str_sub(vecTblNames[1+(i*3)], -30, -21) # assign unique identifier of the current image acquisition field as string
  nameField <- sub(pattern = c("00"), "", nameField) #get rid of some filler '0's net needed in a 96-well format
  nameField <- paste(nameField <- paste(nameField, "3Channels", sep="_"), "csv", sep=".") # paste() inside paste() to save lines
  
  pathTbOutput <- paste(pathOut, nameField, sep="/") # specifies whole output path of created data frame
  write.csv(df_final, pathTbOutput)
}

vecMrgNames <- list.files(path= pathOut, pattern=".csv", full.names = TRUE)
#create empty data frame to which data frames from each field are added
wholeDataSet <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("Area", "MeanInten_b", "MeanInten_g", "MeanInten_r"))

for(k in 1:numFields){
  
  dfWell <- read.csv(vecMrgNames[k], stringsAsFactors = FALSE)
  wholeDataSet <- rbind(wholeDataSet, dfWell)
}

wholeDataSet$X <- NULL

pathEnd <- file.path(".","3_Merged Data", expID, "All Events")
dir.create(pathEnd)

# get rid of StdDev in favor of second MeanInten column. this is a band-aid solution 
names(wholeDataSet)[names(wholeDataSet)=="StdDev_b"] <- "ut_MeanInten_b"
names(wholeDataSet)[names(wholeDataSet)=="StdDev_g"] <- "ut_MeanInten_g"
names(wholeDataSet)[names(wholeDataSet)=="StdDev_r"] <- "ut_MeanInten_r"
wholeDataSet$ut_MeanInten_b <- wholeDataSet$MeanInten_b
wholeDataSet$ut_MeanInten_g <- wholeDataSet$MeanInten_g
wholeDataSet$ut_MeanInten_r <- wholeDataSet$MeanInten_r

fcsSource <- paste(pathEnd, "wholeDataSet.csv", sep="/")
write.csv(wholeDataSet, fcsSource)

# create vector with colnames for plotting etc.
fcsColnames <- colnames(wholeDataSet)[startsWith(colnames(wholeDataSet), "Mean")] # later needed for flowClust

# different plots
count <- nrow(wholeDataSet) # variable for number of rows
plot3d(wholeDataSet[1:count, fcsColnames])
#plot3d(wholeDataSet[1:count, fcsColnames], xlim = c(0, 4000), ylim = c(0, 7000), zlim = c(0, 5000))
#plot3d(wholeDataSet[1:count, fcsColnames], xlim = c(0, 3600), ylim = c(0, 4500), zlim = c(0, 3000))
scatterplot3d(df_final[1:count, fcsColnames])


# read in csv file for fcs creation
fcsData <- fread(fcsSource, check.names = FALSE)
colnames(fcsData)[1] <- "EventNo" # label first column of the table
fcsSource <- gsub(".csv", "", fcsSource) # remove file extension from fcsSource

## Check data quality
dim(fcsData) # number of rows and columns

metadata <- data.frame(name=dimnames(fcsData)[[2]],desc=paste('column',dimnames(fcsData)[[2]],'from dataset'))

## Create FCS file metadata - ranges, min, and max settings
#metadata$range <- apply(apply(data_subset,2,range),2,diff)
metadata$minRange <- apply(fcsData,2,min)
metadata$maxRange <- apply(fcsData,2,max)

fcsData.ff <- new("flowFrame",exprs=as.matrix(fcsData), parameters=AnnotatedDataFrame(metadata)) # in order to create a flow frame, data needs to be read as matrix by exprs
head(fcsData.ff)
write.FCS(fcsData.ff, paste0(fcsSource, ".fcs"))
