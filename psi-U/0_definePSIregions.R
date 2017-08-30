################################################################################
# Title: Define ?? regions to query in multi alignment
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos
################################################################################

bedPeaks <- read.delim("3UTR_peaks.bed",stringsAsFactors = F, header = F)
left <- 8; right <- 12
start <- bedPeaks[,2] - left
end <- bedPeaks[,3] + right
newbedPeaks <- bedPeaks
newbedPeaks[,2] <- start
newbedPeaks[,3] <- end
newbedPeaks <- newbedPeaks[,-c(5,6)]

write.table(x = newbedPeaks, file = "3UTR-PSI_regions.bed", quote = F, row.names = F,
            col.names = F, sep = "\t")
