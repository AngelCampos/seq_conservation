################################################################################
# Title: Process FASTA blocks from FAM into score matrix
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos
################################################################################
#Load/install packages ####
library(stringdist)

# Custom functions ####
# collapseGaps()
# Function to collapse blocks of gaps in two alligned orthologous sequences
collapseGaps <- function(string1, string2){
  if(nchar(string1) != nchar(string2)){
    print("Error: Strings have not same length!")
  }
  else{
    gapblock <- NULL
    for(i in 1:nchar(string1)){
      if(substring(string1, i, i) == "-" & substring(string2, i, i) == "-"){
        gapblock <- c(gapblock, i)
      }
    }
    rString1 <- paste(unlist(strsplit(string1, ""))[-gapblock], collapse = "")
    rString2 <- paste(unlist(strsplit(string2, ""))[-gapblock], collapse = "")
    return (c(rString1, rString2))  
  }
}

# Setting Directories and files ####
iniDir <- getwd()
workDir <- "W:/schwartzlab/miguelg/BIGDATA/UCSC/MAF46_vertebrate/outDir/fastaOut"
setwd(workDir)
files <- list.files()
faFiles <- files[grep(pattern = ".fa", x = files)]

# Generates a list that contains for each peak the distance between the sequence
# of human and the other species
allDist <- list()
for (FILE in faFiles) {
  if (file.size(FILE) > 0){
    x <- read.delim(file = FILE, header = F, stringsAsFactors = F)
    c1 <- x[seq(1,nrow(x), 2),]; c2 <- x[seq(2,nrow(x), 2),]
    speciesList <- sapply(c1, FUN = strsplit, split = ".", fixed = T)
    c1 <- NULL
    for (a in 1:length(speciesList)) {
      c1[a] <- substring(speciesList[a][[1]][1], 2)
    }
    species <- unique(c1)
    flag <- c(grep(pattern = "hg19", c1), length(c1)+1); flag2 <- NULL
    for (i in 1:(length(flag)-1)) {flag2[i] <- flag[i+1] - flag[i]}
    colAssign <- NULL; for (p in 1:length(flag2)){colAssign <- c(colAssign, rep(p, flag2[p]))}
    # Sequence merging by species and order
    seqMatrix <- matrix(NA, nrow = length(species), ncol = table(c1)["hg19"])
    rownames(seqMatrix) <- species
    for (j in 1:length(c1)){
      seqMatrix[c1[j],colAssign[j]] <- c2[j]
    }
    # Correcting empty cells
    for (l in 1:ncol(seqMatrix)) {
      if (sum(is.na(seqMatrix))>0) {
        x <- seqMatrix[,l]
        seqMatrix[,l][which(is.na(x))] <- paste(rep(x = "-", mean(na.omit(sapply(x, nchar)))), collapse = "")  
      }
    }
    # Collapse columns to 1 string
    fasSEQ <- apply(seqMatrix, MARGIN = 1, FUN = paste, collapse = "")
    # Compare sequences with some metric ("Hamming") and build matrix for all results
    DIST <- NULL
    for (k in 1:length(fasSEQ)){
      if (nchar(fasSEQ[1] == nchar(fasSEQ[k]))){
        DIST[k] <- stringdist(a = fasSEQ[1], b = fasSEQ[k], method = "hamming")
      }
      else stop(print("ERROR: Impossible to compare different size sequences"))
    }
    names(DIST) <- names(fasSEQ)
    allDist[[substring(FILE,1, nchar(FILE)-3)]] <- DIST  #Merge in list
  }
  else{
    allDist[[substring(FILE,1, nchar(FILE)-3)]] <- NA 
  }
}

# Species in allDist
allSpecies <- NULL
for (s in 1:length(allDist)){
  allSpecies <- unique(c(allSpecies, names(allDist[[s]])))
}

# Construct Distance matrix: Rows = species, Columns = Peaks
distMATRIX <- matrix(data = NA, nrow = length(allSpecies), ncol = length(allDist))
rownames(distMATRIX) <- allSpecies; colnames(distMATRIX) <- names(allDist)
for(z in names(allDist)){
  distMATRIX[,z] <- allDist[[z]][rownames(distMATRIX)]
}

na_count <- apply(distMATRIX, function(y) sum(is.na(y)) , MARGIN = 2)
d_count <- length(allSpecies) - na_count

# Go back to initial working directory
setwd(iniDir)
# Checkpoint! ####
save.image("processFasta.Rdata")
# load("processFasta.Rdata")