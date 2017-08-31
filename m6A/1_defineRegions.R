# Define regions to query in multi alignment ###################################
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos ############

# Functions ####################################################################
# Saves tab separated table
mytable <- function(X, fileName){
  write.table(X,file = fileName, sep = "\t", col.names = F, row.names = F, quote = F)
}
# Function for generating a region around the first position of an input BED file
genRegion <- function(bedPeaks, d5, d3){
  start <- bedPeaks[,2] - d5
  end <- bedPeaks[,2] + d3 +1
  newbedPeaks <- bedPeaks
  newbedPeaks[,2] <- start
  newbedPeaks[,3] <- end
  return(newbedPeaks)
}

# Generate regions for BED files
BEDFILES <- c("m6Apeaks_human_3UTR.bed", "m6Apeaks_human_5UTR.bed")
d5 <- 5; d3 <- 5
for(i in BEDFILES){
  bedPeaks <- read.delim(i ,stringsAsFactors = F, header = F)
  newbedPeaks <- genRegion(bedPeaks, d5, d3)
  mytable(newbedPeaks, paste(gsub('.{4}$', '', i), "_regions.bed"))
}

## TEST SET #################
bedPeaks <- read.delim("testPeaks_human.bed" ,stringsAsFactors = F, header = F)
nBED <- genRegion(bedPeaks, 5, 5)
mytable(nBED, "testPeaks_human_regions.bed")
