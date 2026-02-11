#!/bin/bash

# @pepap : cutAdapt version 2.5
metamodule add py-cutadapt
metamodule add fastx-0.0.14

CPU=8
STOR="/storage/brno12-cerit/home/pepap/brno1/Valeria.Buccheri/22.mESC--Dicer-loop-mutants--TraPR--20250516/FASTQ/2025_05_15_Buccheri_mESCs_loops_TraPR"
NRMV=`echo ${STOR} | wc -m`
SUFF="_L001_R1_001.fastq.gz"
FASTQ=( $(/bin/ls ${STOR}/*/*${SUFF} | colrm 1 ${NRMV} | sed "s/_ds[.].*$//g") )

# @pepap : NEXTflex Small RNA-Seq kit v3

#!#ADAP5P_X="GTTCAGAGTTCTACAGTCCGACGATCNNNN"
#!#ADAP3P_X="NNNNTGGAATTCTCGGGTGCCAAGG"
#!#ADAP5P_B="GTTCAGAGTTCTACAGTCCGACGATC"
#!#ADAP3P_B="AGATCGGAAGAGCACACGTCT"
# @pepap : 5'-adaptor sequence
ADAP5P_N="TCTTTCCCTACACGACGCTCTTCCGATCT"
# @pepap : 3'-adaptor sequence
ADAP3P_N="TGGAATTCTCGGGTGCCAAGG"
#!#TRUSEQ3P="GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
# @pepap : error rate
ERRLIM=0.075
# @pepap : minimal overlap
MINOVR=14

printf "LIBNAME,TOT,TRIMMED,UNTRIMMED\n" >  totN.dat

for ((i=1;i<=${#FASTQ[@]};i++))
do

f_i=`printf "%03d" ${i}`
j=$(($i-1))
inpSEQ=`/bin/ls ${STOR}/${FASTQ[${j}]}_ds.*/*${SUFF}`

echo ""
echo " ++ ${FASTQ[${j}]} : ${inpSEQ} ++"
echo ""

echo -n "${FASTQ[${j}]},"                >> totN.dat
TMPNUM=`zcat ${inpSEQ} | wc -l | awk '{ print $1/4 }'`
echo -n "${TMPNUM},"                     >> totN.dat

cutadapt --cores ${CPU} \
         --front=${ADAP5P_N} \
         --adapter=${ADAP3P_N} \
         --error-rate=${ERRLIM} \
         --times=2 \
         --overlap=${MINOVR} \
         --match-read-wildcards \
         --minimum-length=12 \
         --action trim \
         --trim-n \
         --max-n=0.1 \
         --output="${FASTQ[${j}]}-trim.fq" \
         --untrimmed-output="${FASTQ[${j}]}-untr.fq" \
         ${inpSEQ}

TMPNUM=`wc -l ${FASTQ[${j}]}-trim.fq | awk '{ print $1/4 }'`
echo -n "${TMPNUM},"                     >> totN.dat
TMPNUM=`wc -l ${FASTQ[${j}]}-untr.fq | awk '{ print $1/4 }'`
echo    "${TMPNUM}"                      >> totN.dat

cat ${FASTQ[${j}]}-trim.fq ${FASTQ[${j}]}-untr.fq | fastx_collapser -v -i - -o ${FASTQ[${j}]}-all.re-collapsed.fa

rm -f tmp.${f_i}.fastq
gzip -fv ${FASTQ[${j}]}-trim.fq ${FASTQ[${j}]}-untr.fq ${FASTQ[${j}]}-all.re-collapsed.fa

echo ""

done

