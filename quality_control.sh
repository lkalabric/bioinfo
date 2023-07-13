#!/bin/bash

# command name: quality_control.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 14 JUN 2023
# objetive: Give examples of quality control apps
# Syntax: ./quality_control.sh

# This is the very first step after basecalling and demux steps which are undergone in apps from each NGS tech.
# There are different methods in bioinfo to do quality control. Let's see some of them: 


# Validating arguments
if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters"
    echo "Syntax: quality_control.sh <-illumina | -minion> <SAMPLE_ID>"
    exit 0    
fi

# Declaring variables
SAMPLE_DIR=$2
INPUT_DIR="${HOME}/data/hbv/${SAMPLE_ID}"
if [ ! -d ${INPUT_DIR} ]; then
    echo "Sample not identified. Using example data 0001.1 instead!"
    INPUT_DIR="${HOME}/data/hbv/0001.1" # If a bash variable is empty, let's use an example data
fi
#OUTPUT_DIR="${HOME}/qc-results/${SAMPLE_ID}"
OUTPUT_DIR="${HOME}/bioinfo-results/${SAMPLE_ID}/qc-control"
[ -d ${OUTPUT_DIR} ] || mkdir ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

case $1 in
  "-illumina")
    # Quality control only
    
    ##
    # Illumina data
    ##
    # 1) Fastqc
    # Link: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    # Requirements: Java
    # Installation:
    # $ sudo apt install default-jre
    # $ sudo apt install fastqc
    fastqc -o ${OUTPUT_DIR} -f fastq -c ${INPUT_DIR}/*.fastq
    
    # 2) Fastqcr 
    # Link: https://rpkgs.datanovia.com/fastqcr/index.html
    # Requirements: R
    # Installation:
    # Being developed
    # Rscript fastqcr-analysis.R # Need to be written!!!
    
    # 3) Afterqc
    # Link: https://github.com/OpenGene/AfterQC
    # Requirements: Miniconda (Python)
    # Installation:
    # $ install_thirdparty.sh
    source activate afterqc
    # Quality control only
    cd ${OUTPUT_DIR}
    # Single-ended analysis
    after.py --qc_only -d ${INPUT_DIR} -1 *R1*
    after.py --qc_only -d ${INPUT_DIR} -2 *R2*
    # Pair-ended analysis
    after.py --qc_only -d ${INPUT_DIR} -1 *R1* -2 *R2*

;;
  "-minion")
  ##
  # MinIon data
  ##
  # 1) pycoQC
  # Link: https://hpc.nih.gov/apps/pycoQC.html
  # Requirements: Miniconda (Python)
  # Installation:
  # $ install_thirdparty.sh
  # Or: If you already have Conda env, run this
  # $ conda create -n pycoqc
  # $ source activate pycoqc
  # $ conda install -c bioconda pycoqc
  # source activate pycoqc
  # The next step continues in quality_filter.sh
 ;;
  *)
    echo "Invalid parameter!"
    exit 1
esac
