# Methodology for evaluating the conservation of genomic intervals 
# through the use of MAF files and genomic coordinates from previous studies

# 1. Generate BED files that delimite the region we want to extract from the
# MAF file based on the human genome

# 2. Filter out by specific feature (UTR, exon, intron, etc...) OPTIONAL
# bedtools intersect -a peaks.sorted.bed -b 3UTR_hg19_UCSCgenes.bed -u > 3UTR_peaks.bed

# 3. Extract the blocks of multiple alignment
WORK_DIR=/home/labs/schwartzlab/miguelg/WORKSPACE_wexac/Exploration/seq_conservation/m6A/
REGIONS=m6Apeaks_regions_human.sorted.bed
MAFS=mafs.txt
MAF_DIR=/home/labs/schwartzlab/miguelg/BIGDATA/UCSC/MAF46_vertebrate/
MAF_OUTDIR=MAFblocks

cat "${MAF_DIR}${MAFS}" | while read line
do
  mafsInRegion "${WORK_DIR}${REGIONS}" "${MAF_OUTDIR}" "${MAF_DIR}${line}" -outDir -keepInitialGaps 
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
