# Define regions to query in multi alignment ###################################
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos ############

# Functions ####################################################################
# Saves tab separated table
mytable <- function(X, fileName){
  write.table(X,file = fileName, sep = "\t", col.names = F, row.names = F, quote = F)
}

## TEST SET ####################################################################
# Sample of human peaks BED file
humanBED <- read.delim("m6Apeaks_human_3UTR.bed" ,stringsAsFactors = F, header = F)
sampleBED <- humanBED[sample(1:nrow(humanBED), size = 100, replace = F), ]
mytable(sampleBED, "testPeaks_human.bed")