#!/bin/bash

# command name: assembly_de_novo.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 14 JUN 2023
# objetive: Give examples of assembly apps
# Syntax: ./assembly_de_novo.sh

# Validating arguments
if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters"
    echo "Syntax: assembly_denovo.sh <-illumina | -minion> <SAMPLE_ID>"
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
OUTPUT_DIR="${HOME}/bioinfo-results/${SAMPLE_ID}/assembly_denovo"
[ -d ${OUTPUT_DIR} ] || mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

REFSEQ="${HOME}/data/REFSEQ/hbv/NC_003977.2.fasta"

case $1 in
  "-illumina")
    # De Novo Assembly
    # Decompress input files from qc-filter dir
    gzip -d ${INPUT_DIR}/*.gz
    # 1) Use of spades
    # Link: https://github.com/ablab/spades/blob/spades_3.15.5/README.md
    # They recommend running SPAdes with BayesHammer/IonHammer to obtain high-quality assemblies.
    # Note: We decided to avoid unpaired reads!!!!
    # For single lib use -1 and -2
    # spades -1 ${INPUT_DIR}/output_forward_paired.fq -2 ${INPUT_DIR}/output_reverse_paired.fq -o ${OUTPUT_DIR}
    # For single lib all reads paired and unpaired use -s
    # spades -1 ${INPUT_DIR}/output_forward_paired.fq -2 ${INPUT_DIR}/output_reverse_paired.fq -s ${INPUT_DIR}/output_forward_unpaired.fq -s ${INPUT_DIR}/output_reverse_unpaired.fq -o ${OUTPUT_DIR}
    # Assembly by reference
    # spades -1 ${INPUT_DIR}/output_forward_paired.fq -2 ${INPUT_DIR}/output_reverse_paired.fq -o ${OUTPUT_DIR} --trusted-contigs ${REFSEQ}
    # 2) Use of bwa
    bwa index ${REFSEQ}
    bwa mem ${REFSEQ} ${INPUT_DIR}/output_forward_paired.fq ${INPUT_DIR}/output_reverse_paired.fq gzip -3 > aln-pe.sam.gz
    
;;
  "-minion")
  ##
  # MinIon data
  ##
  # De Novo Assembly
  # Use of spades
  # source activate spades
  # spades.py -o ${OUTPUT_DIR} --nanopore ${INPUT_DIR}/Filename
    # Assembly by reference
    # Use of bwa

# By reference
# 1) Use of minimap
# Link: https://timkahlke.github.io/LongRead_tutorials/ASS_M.html

  
 ;;
  *)
    echo "Invalid parameter!"
    exit 1
esac

