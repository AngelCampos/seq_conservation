################################################################################
# Title: Define ?? regions to query in multi alignment
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos
################################################################################

bedPeaks <- read.delim("m6Apeaks_human.bed",stringsAsFactors = F, header = F)
left <- 5; right <- 5
start <- bedPeaks[,2] - left
end <- bedPeaks[,3] + right
newbedPeaks <- bedPeaks
newbedPeaks[,2] <- start
newbedPeaks[,3] <- end

write.table(x = newbedPeaks, file = "m6Apeaks_regions_human.bed", quote = F, row.names = F,
            col.names = F, sep = "\t")
