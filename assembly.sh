#!/bin/bash

# command name: assembly.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 14 JUN 2023
# objetive: Give examples of assembly apps
# Syntax: ./assembly.sh

# Validating arguments
if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters"
    echo "Syntax: assembly.sh <-illumina | -minion> <SAMPLE_ID>"
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
OUTPUT_DIR="${HOME}/bioinfo-results/${SAMPLE_ID}/assembly"
[ -d ${OUTPUT_DIR} ] || mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

case $1 in
  "-illumina")
    # De Novo Assembly
    # Decompress input files from qc-filter dir
    gzip -d ${INPUT_DIR}/*.gz
    # Use of spades
    # Link: https://github.com/ablab/spades/blob/spades_3.15.5/README.md
    # They recommend running SPAdes with BayesHammer/IonHammer to obtain high-quality assemblies.  
    source activate spades
    spades.py -s ${INPUT_DIR}/output_forward_paired.fq --pe1-2 ${INPUT_DIR}/output_reverse_paired.fq -o ${OUTPUT_DIR}
    
    # Assembly by reference
    # Use of bwa    
    
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

  
 ;;
  *)
    echo "Invalid parameter!"
    exit 1
esac

