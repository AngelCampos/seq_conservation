# Title: #######################################################################
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos ############

# Packages #####################################################################
library(stringdist)
library(magrittr)
library(RColorBrewer)

# Functions ####################################################################

# Reversing a string x
strReverse <- function(x){
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")  
}

# Remove last n characters from string x
rmTailString <- function(x, n){
  strReverse(substring(strReverse(x), n+1))
}

# Read MultiFasta
rFASTA <- function(filename){
  return(read.delim(file = filename, header = F, stringsAsFactors = F))
}

# Hamming distance of the input MULTIFASTA given in x
seqDistMulFasta <- function(x, refSpecies = "hg19"){
  c1 <- x[seq(1,nrow(x), 2),]; c2 <- x[seq(2,nrow(x), 2),]
  speciesList <- sapply(c1, FUN = strsplit, split = ".", fixed = T)
  c1 <- NULL
  for (a in 1:length(speciesList)) {
    c1[a] <- substring(speciesList[a][[1]][1], 2)
  }
  species <- unique(c1)
  flag <- c(grep(pattern = "hg19", c1), length(c1)+1); flag2 <- NULL
  for (i in 1:(length(flag)-1)) {flag2[i] <- flag[i+1] - flag[i]}
  colAssign <- NULL; for (p in 1:length(flag2)){
    colAssign <- c(colAssign, rep(p, flag2[p]))}
  seqMatrix <- matrix(NA, nrow = length(species), ncol = table(c1)["hg19"])
  rownames(seqMatrix) <- species
  for (j in 1:length(c1)){
    seqMatrix[c1[j],colAssign[j]] <- c2[j]
  }
  # Correcting empty cells
  for (l in 1:ncol(seqMatrix)) {
    if (sum(is.na(seqMatrix)) > 0) {
      x <- seqMatrix[,l]
      seqMatrix[,l][which(is.na(x))] <- paste(rep(x = "-", mean(na.omit(sapply(x, nchar)))), collapse = "")  
    }
  }
  # Collapse columns to 1 string
  fasSEQ <- apply(seqMatrix, MARGIN = 1, FUN = paste, collapse = "") %>% toupper()
  # Compare sequences with some metric ("Hamming") and build matrix for all results
  DIST <- stringdist(a = fasSEQ[refSpecies], b = fasSEQ, method = "hamming")
  names(DIST) <- names(fasSEQ)
  return(DIST)
}

# Distance Matrix of all files in fileList
distMATRIX <- function(fileList){
  allDist <- list()
  for (FILE in fileList) {
    if (file.size(FILE) > 0){
      DIST <- rFASTA(FILE) %>% seqDistMulFasta()
      allDist[[substring(FILE,1, nchar(FILE)-3)]] <- DIST  #Merge in list
    }
    else{
      allDist[[substring(FILE,1, nchar(FILE)-3)]] <- NA 
    }
  }
  allSpecies <- NULL # Species in allDist
  for (s in 1:length(allDist)){
    allSpecies <- unique(c(allSpecies, names(allDist[[s]])))
  }
  # Construct Distance matrix: Rows = species, Columns = Peaks
  distMATRIX <- matrix(data = NA, nrow = length(allSpecies), ncol = length(allDist))
  rownames(distMATRIX) <- allSpecies
  colnames(distMATRIX) <- names(allDist)
  for(z in names(allDist)){
    distMATRIX[,z] <- allDist[[z]][rownames(distMATRIX)]
  }
  return(distMATRIX)
}

# MAIN PROGRAM #################################################################
# Setting up directories and files 
iniDir <- getwd()
workDir <- paste(iniDir, "/fastaSeqs", sep ="")
setwd(workDir)
files <- list.files()
faFiles <- files[grep(pattern = ".fa", x = files)]

# File Names
ranControls <- sort(c(grep("R1", faFiles), grep("R2", faFiles), grep("R3", faFiles)))
peaks <- faFiles[-ranControls]
controls <- faFiles[ranControls] 
ranCon1 <- faFiles[grep("R1", faFiles)]
ranCon2 <- faFiles[grep("R2", faFiles)]
ranCon3 <- faFiles[grep("R3", faFiles)]

# Calculate empirically the mean distance between each species and human
spDistControl <- distMATRIX(controls) %>% apply(1, na.omit)
medDistControl <- sapply(spDistControl, FUN =  median)
meanDistControl <- sapply(spDistControl, FUN =  mean)
boxplot(spDistControl, las = 2, col= brewer.pal(11, "Spectral"), main = "Sequence distance - Controls")
barplot(medDistControl, las=2, col= brewer.pal(11, "Spectral"), main = "Median sequence distance - Controls")
barplot(meanDistControl, las=2, col= brewer.pal(11, "Spectral"), main = "Mean sequence distance - Controls")

spDistPeaks <- distMATRIX(peaks) %>% apply(1, na.omit) %>% sapply(FUN =  median)
barplot(spDistPeaks, las=2, col= brewer.pal(11, "Spectral"), main = "Median sequence distance - Peaks")

# Score Peaks for conservation of sequence

i <- 7
x <- rFASTA(peaks[i]) %>% seqDistMulFasta()
o <- rFASTA(ranCon1[i]) %>% seqDistMulFasta()
p <- rFASTA(ranCon2[i]) %>% seqDistMulFasta()
q <- rFASTA(ranCon3[i]) %>% seqDistMulFasta()

intSpec <- Reduce(f = intersect,x = list(names(x), names(o), names(p), names(q)))
RESMAT <- cbind(x[intSpec], o[intSpec], p[intSpec], q[intSpec])
colnames(RESMAT) <- c("Peak", "Rnd1", "Rnd2", "Rnd3")
gplots::heatmap.2(RESMAT)
colMeans(RESMAT)


row.distance = dist(RESMAT, method = "euclidean")
row.cluster = hclust(row.distance, method = "average")
plot(row.cluster)

col.distance = dist(t(RESMAT), method = "euclidean")
col.cluster = hclust(col.distance, method = "average")
plot(col.cluster)




