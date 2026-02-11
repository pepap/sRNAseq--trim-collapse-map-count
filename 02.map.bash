#!/bin/bash

module add samtools

# @pepap : STAR-aligner executable
exeCMD="STAR"
# @pepap : REFERENCE sequence
REF="/path/to/GRCm38.fa"
# @pepap : STAR index repository
IND="/path/to/STAR_index"
# @pepap :  input reads' repository
STOR="/path/to/input/reads"
# @pepap : BASH array containing the sample files
FASTQS=(
  1  2  3
  4  5  6
  ...
  ...
  ...
)
# @pepap : BASH arrray containing the output PREFIX names
OUTBAM=(
 TraPR-dR-C5 TraPR-dR-A1 TraPR-dR-A4
 TraPR-dG-A1 TraPR-dG-B2 TraPR-dG-C4
 TraPR-dO-A2 TraPR-dO-C6 TraPR-dO-D1
 ...
 ...
 ...
)
# @pepap : SUFFIX of the input sample files
SUFF="-all.re-collapsed.fa.gz"

i=__INDEX__

${exeCMD} --runMode alignReads \
          --runThreadN  4 \
          --genomeDir ${IND} \
          --genomeLoad LoadAndRemove \
          --readFilesIn ${FASTQS[$i]}${SUFF} \
          --readFilesPrefix ${STOR}/ \
          --readFilesCommand zcat \
          --limitBAMsortRAM  20000000000 \
          --outFileNamePrefix ${OUTBAM[$i]}.se. \
          --outReadsUnmapped Fastx \
          --outSAMtype BAM SortedByCoordinate \
          --outFilterMultimapNmax 99999 \
          --outFilterMismatchNoverLmax 0.1 \
          --outFilterMatchNminOverLread 0.66 \
          --alignSJoverhangMin   999 \
          --alignSJDBoverhangMin 999 

samtools index ${OUTBAM[$i]}.se.Aligned.sortedByCoord.out.bam

gzip   ${OUTBAM[$i]}.se.Unmapped.out.mate1

