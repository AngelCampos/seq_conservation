# Methodology for evaluating the conservation of genomic intervals 
# through the use of MAF files and genomic coordinates from previous studies

# 0. Generate BED files from supplementariy table
Rscript 0_GenerateBEDfiles.R
bedtools sort -i m6Apeaks_human.bed > m6Apeaks_human.sorted.bed 
bedtools sort -i m6Apeaks_mouse.bed > m6Apeaks_mouse.sorted.bed 

# 2. Filter out by specific feature (UTR, exon, intron, etc...) OPTIONAL
annot1=$HOME/BIGDATA/UCSC/BED_annotations/3UTR_hg19_UCSCgenes.bed
annot2=$HOME/BIGDATA/UCSC/BED_annotations/5UTR_hg19_UCSCgenes.bed

bedtools intersect -a m6Apeaks_human.sorted.bed -b $annot1 -u > m6Apeaks_human_3UTR.bed
bedtools intersect -a m6Apeaks_human.sorted.bed -b $annot2 -u > m6Apeaks_human_5UTR.bed

# 1. Generate BED files that delimite the region we want to extract from the
# MAF file based on the human genome
Rscript 1.0_sampleTestPeaks.R
Rscript 1_defineRegions.R

# 2. Generate Random regions for controls
Rscript 2_generateRanRegions.R

## Testing
bedtools sort -i testPeaks_human_regions.bed > testPeaks_human_regions.sort.bed
bedtools sort -i randomRegions_testingPeaks_human.bed > test_RanRegions_human.sort.bed
# 3. Extract the blocks of multiple alignment
WORK_DIR=/$HOME/WORKSPACE_wexac/seq_conservation/m6A/
REGIONS=testPeaks_human_regions.sort.bed
REGIONS2=test_RanRegions_human.sort.bed
MAFS=mafs.txt
MAF_DIR=/$HOME/BIGDATA/UCSC/MAF46_vertebrate/
MAF_OUTDIR=MAFblocks

cat "${MAF_DIR}${MAFS}" | while read line
do
  mafsInRegion "${WORK_DIR}${REGIONS}" "${MAF_OUTDIR}" "${MAF_DIR}${line}" -outDir -keepInitialGaps 
done &

cat "${MAF_DIR}${MAFS}" | while read line
do
  mafsInRegion "${WORK_DIR}${REGIONS2}" "${MAF_OUTDIR}" "${MAF_DIR}${line}" -outDir -keepInitialGaps 
done &

# 4. Convert to FASTA sequence # Dependencies: bx-python
FASTA_OUTDIR=fastaSeqs 
mkdir "${WORK_DIR}${FASTA_OUTDIR}"
cd $MAF_OUTDIR	 
for file in *.maf
do
  maf_to_fasta.py < $file > "${WORK_DIR}${FASTA_OUTDIR}/${file/%maf/}fa"
done &

# 5. Perform analytic methodology
