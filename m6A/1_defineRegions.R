# Define regions to query in multi alignment ###################################
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos ############

# Function
# Saves tab separated table
mytable <- function(X, fileName){
  write.table(X,file = fileName, sep = "\t", col.names = F, row.names = F, quote = F)
}

# Generate regions for BED files
BEDFILES <- c("m6Apeaks_human_3UTR.bed", "m6Apeaks_human_5UTR.bed")

d5 <- 5; d3 <- 5
for(i in BEDFILES){
  bedPeaks <- read.delim(i ,stringsAsFactors = F, header = F)
  start <- bedPeaks[,2] - d5
  end <- bedPeaks[,3] + d3
  newbedPeaks <- bedPeaks
  newbedPeaks[,2] <- start
  newbedPeaks[,3] <- end
  write.table(x = newbedPeaks, 
              file = paste(gsub('.{4}$', '', i), "_regions.bed"), 
              quote = F, row.names = F, col.names = F, sep = "\t")
}

# Sample of human peaks BED file

humanBED <- read.delim("m6Apeaks_human_3UTR.bed" ,stringsAsFactors = F, header = F)
sampleBED <- humanBED[sample(1:nrow(humanBED), size = 100, replace = F), ]
mytable(sampleBED, "testingPeaks_human.bed")
