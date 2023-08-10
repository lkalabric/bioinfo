#!/bin/bash

# command name: assembly_by_reference.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 14 JUN 2023
# objetive: Give examples of assembly apps
# Syntax: ./assembly_byreference.sh

# Validating arguments
if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters"
    echo "Syntax: assembly_byreference.sh <-illumina | -minion> <SAMPLE_ID>"
    exit 0    
fi

# Declaring variables
SAMPLE_ID=$2
INPUT_DIR="${HOME}/bioinfo-results/${SAMPLE_ID}/qc-filter"
if [ ! -d ${INPUT_DIR} ]; then
    echo "Qc-filter results absent. Using qc-filter results from sample 0001.1 instead!"
    INPUT_DIR="${HOME}/bioinfo-results/0001.1/qc-filter" # If a bash variable is empty, let's use an example data
fi
#OUTPUT_DIR="${HOME}/qc-results/${SAMPLE_ID}"
OUTPUT_DIR="${HOME}/bioinfo-results/${SAMPLE_ID}/assembly_byreference"
[ -d ${OUTPUT_DIR} ] || mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

REFSEQ="${HOME}/data/REFSEQ/hbv/NC_003977.2.fasta"

case $1 in
  "-illumina")
  # 1) Use of bwa
  bwa index ${REFSEQ}
  bwa mem ${REFSEQ} ${INPUT_DIR}/output_forward_paired.fq ${INPUT_DIR}/output_reverse_paired.fq gzip -3 > aln-pe.sam.gz
    
;;
  "-minion")
# 1) Use of minimap
# Link: https://timkahlke.github.io/LongRead_tutorials/ASS_M.html

  
 ;;
  *)
    echo "Invalid parameter!"
    exit 1
esac
