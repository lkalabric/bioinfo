#!/bin/bash

# command name: quality_control.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 14 JUN 2023
# objetive: Give examples of quality control apps
# Syntax: ./quality_control.sh
# 

# There are different methods in bioinfo to do the quality control. Let's see some of them: 

# Example data
SAMPLE="0001.1"
INPUT_DIR="${HOME}/data/hbv/${SAMPLE}"
OUTPUT_DIR="${HOME}/qc-results/${SAMPLE}"
[ -d ${OUTPUT_DIR} ] || mkdir ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

# Preliminary data quality

# 1) Fastqc
# Link: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# Requirements: Java
# Installation:
# $ sudo apt install default-jre
# $ sudo apt install fastqc
fastqc -o ${OUTPUT_DIR} -f fastq -c ${INPUT_DIR}/*.fastq

# 2) Fastqrc 
# Link: https://rpkgs.datanovia.com/fastqcr/index.html
# Requirements: R
# Installation:


# 3) Afterqc
# Link: https://github.com/OpenGene/AfterQC
# Requirements: Miniconda (Python)
# Installation:
# $ install_thirdparty.sh
source activate afterqc
after.py --qc_only -1 ${INPUT_DIR}/*R1* -2 ${INPUT_DIR}/*R2*

# Quality filtering

# Use of timmomatic for Illumina data


# Use of Nanofilt for MinIon data
