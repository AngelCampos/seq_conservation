# Title: Generate random regions ###############################################
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos ############

# Get annotation and BED file with reference peaks
annot <- read.delim(file = "~/BIGDATA/UCSC/BED_annotations/3UTR_hg19_UCSCgenes.bed",
                    header = F, stringsAsFactors = F)
bed <- read.delim("testPeaks_human.bed", header = F, stringsAsFactors = F)

# Functions ####################################################################
# Generate a matrix with N number of random regions (rows) based on a BED format
# position within the specified annotation
genRanRegions <- function(x, d5, d3, nControls, annot){
  chr <- as.character(x[1])
  coor <- as.numeric(x[2])
  name <- as.character(x[4])
  sub <- annot[annot$V1 == chr,]
  inside <- NULL
  for(i in 1:nrow(sub)){
    inside[i] <- coor >= sub[i,2] & coor <= sub[i,3]
  }
  rangeRegion <- NULL
  lim5 <- min(sub[inside,2]) + d5
  lim3 <- max(sub[inside,3]) - d3
  if((coor - d5) - lim5 > 0){
    rangeRegion <- c(rangeRegion, lim5:(coor-d5*2))
  }
  if((coor + d3) - lim3 < 0){
    rangeRegion <- c(rangeRegion, (coor+d3*2):lim3)
  }
  nPeaks <- sample(rangeRegion, nControls, replace = F)
  nregions <- matrix(x, nrow = length(nPeaks), ncol = length(x), byrow = T)
  nregions[,2] <- nPeaks - d5
  nregions[,3] <- nPeaks + d3 + 1
  nregions[,5] <- 0
  nregions[,4] <- sapply(name, paste0,".R", 1:length(nPeaks))
  return(nregions)
}

# Saves tab separated table
mytable <- function(X, fileName){
  write.table(X,file = fileName, sep = "\t", col.names = F, row.names = F, quote = F)
}

# Executing ####################################################################
nCont <- 3
RegionsTable <- NULL
for (i in 1:nrow(bed)){
  RegionsTable <- rbind(RegionsTable, genRanRegions(bed[i,], 5, 5, nCont, annot))
}

mytable(RegionsTable, "randomRegions_testingPeaks_human.bed")
