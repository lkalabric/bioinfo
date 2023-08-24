#!/bin/bash

# command name: quality_filter.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 15 JUN 2023
# objetive: Give examples of quality filter apps
# Syntax: ./quality_filter.sh

# Validate arguments
if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters"
    echo "Syntax: quality_filter.sh <-illumina | -minion> <SAMPLE_ID> " 
    exit 0    
fi

# Declaring variables
SAMPLE_ID=$2
INPUT_DIR="${HOME}/data/hbv/${SAMPLE_ID}"
if [ ! -d ${INPUT_DIR} ]; then
    echo "Sample not identified. Using example data 0001.1 instead!"
    INPUT_DIR="${HOME}/data/hbv/0001.1" # If a bash variable is empty, let's use an example data
fi
#OUTPUT_DIR="${HOME}/qc-filter/${SAMPLE_ID}"
OUTPUT_DIR="${HOME}/bioinfo-results/${SAMPLE_ID}/quality-filter"
[ -d ${OUTPUT_DIR} ] || mkdir -p ${OUTPUT_DIR}
MINLENTGH=100

case $1 in
  "-illumina")
      ##
    # Illumina data
    ##
    # 1) Afterqc
    # Link: https://github.com/OpenGene/AfterQC
    # Requirements: Miniconda (Python)
    # Installation:
    # $ install_thirdparty.sh
    # source activate afterqc
    # Single-ended analysis
    # after.py -d ${INPUT_DIR} -1 *R1*
    # after.py -d ${INPUT_DIR} -2 *R2*
    # Pair-ended analysis
    # after.py -d ${INPUT_DIR} -1 *R1* -2 *R2*
    
    # 2) Trimmomatic
    # Use: Filtering and trimming Illumina data
    # Link: http://www.usadellab.org/cms/?page=trimmomatic
    # java -jar trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
    source activate trimmomatic
    trimmomatic PE ${INPUT_DIR}/*R1*.fastq.gz ${INPUT_DIR}/*R2*.fastq.gz ${OUTPUT_DIR}/output_forward_paired.fq.gz ${OUTPUT_DIR}/output_forward_unpaired.fq.gz ${OUTPUT_DIR}/output_reverse_paired.fq.gz ${OUTPUT_DIR}/output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:${MINLENTGH}

    # 3) Prinseq 
;;
  "-minion")
    ##
    # MinIon data
    ##
    # 1) Nanofilt
    # Use: Filtering and trimming MinIon data
    # Link: https://github.com/wdecoster/nanofilt
    # Requirements: Miniconda (Python)
    # Installation:
    # $ install_thirdparty.sh
    # Or: If you already have Conda env, run this
    # $ conda create -n nanofilt
    # $ source activate nanofilt
    # $ conda install -c bioconda nanofilt
    # source activate nanofilt
    # NanoFilt -l 100 -q 9 ...
    # The next step continues in other_filters.sh

    # 2) Porechop -barecode_threshold 85
    # Link: https://github.com/rrwick/Porechop
  ;;
  *)
  echo "Invalid parameter!"
  exit 1
esac
  







