# Generate Bed Files from supplementary tables #################################
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos

# Function
mytable <- function(X, fileName){
  write.table(humanBED,file = fileName, sep = "\t", col.names = F, row.names = F, quote = F)
}
# Human BED file
humantxt <- read.delim("m6Apeaks_human.txt", stringsAsFactors = F)
chr <- humantxt[,1]
pos <- humantxt[,2]
strand <- humantxt[,3]
cScore <- humantxt[,"Classifier.score"]
peakName <- as.vector(sapply("human_m6Apeak.", paste0, 1:length(chr)))
humanBED <- cbind(chr, pos, pos, peakName, cScore, strand)
mytable(humanBED, "m6Apeaks_human.bed")

# Mouse BED file
mousetxt <- read.delim("m6Apeaks_mouse.txt", stringsAsFactors = F)
chr <- mousetxt[,1]
pos <- mousetxt[,2]
strand <- mousetxt[,3]
cScore <- mousetxt[,9]
peakName <- as.vector(sapply("mouse_m6Apeak.", paste0, 1:length(chr)))
humanBED <- cbind(chr, pos, pos, peakName, cScore, strand)
mytable(mouseBED, "m6Apeaks_mouse.bed")