#!/bin/bash

# command name: assembly_by_reference.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 14 JUN 2023
# objetive: Give examples of assembly apps
# Syntax: ./assembly_by_reference.sh

# Validating arguments
if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters"
    echo "Syntax: assembly_by_reference.sh <-illumina | -minion> <SAMPLE_ID>"
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
OUTPUT_DIR="${HOME}/bioinfo-results/${SAMPLE_ID}/assembly_by_reference"
[ -d ${OUTPUT_DIR} ] || mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

REFSEQ="${HOME}/data/REFSEQ/hbv/NC_003977.2.fasta"

case $1 in
  "-illumina")
  # 1) Use of bwa
  # Link: https://github.com/lh3/bwa
  # Link: https://www.htslib.org/workflow/fastq.html
  bwa index ${REFSEQ}
  # bwa mem ${REFSEQ} ${INPUT_DIR}/output_forward_paired.fq ${INPUT_DIR}/output_reverse_paired.fq | gzip -3 > aln-pe.sam.gz
  # Mapping
  bwa mem ${REFSEQ} ${INPUT_DIR}/output_forward_paired.fq ${INPUT_DIR}/output_reverse_paired.fq | \ 
  samtools view -bS -F4 - | \
  samtools sort - -o ${OUTPUT_DIR}/out.bam
;;
  "-minion")
# 1) Use of minimap
# Link: https://timkahlke.github.io/LongRead_tutorials/ASS_M.html
# Link: https://www.htslib.org/workflow/fastq.html
    minimap2 -t 8 -a -x sr ${REFSEQ} ${INPUT_DIR}/output_forward_paired.fq ${INPUT_DIR}/output_reverse_paired.fq  | \
    samtools fixmate -u -m - - | \
    samtools sort -u -@2 -T /tmp/example_prefix - | \
    samtools markdup -@8 --reference ${REFSEQ} - final.cram
 ;;
  *)
    echo "Invalid parameter!"
    exit 1
esac
