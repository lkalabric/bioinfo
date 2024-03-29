#!/bin/bash

# command name: assembly_de_novo.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 14 JUN 2023
# objetive: Give examples of de novo assembly apps
# Syntax: ./assembly_de_novo.sh

# montadores mais utilizados para análise transcrptômica: trinity, spades, soapdenovo-trans
# montadores mais utilizados para genômica: spades, velvet

# Validating arguments
if [[ $# -ne 2 ]]; then
    echo "Illegal number of parameters"
    echo "Syntax: assembly_denovo.sh <-spades | -velvet | -minimap> <SAMPLE_ID>"
    exit 0    
fi

# Declaring variables
SAMPLE_ID=$2
INPUT_DIR="${HOME}/bioinfo-results/${SAMPLE_ID}/quality-filter"
if [ ! -d ${INPUT_DIR} ]; then
    echo "Qc-filter results absent. Using qc-filter results from sample 0001.1 instead!"
    INPUT_DIR="${HOME}/bioinfo-results/0001.1/quality-filter" # If a bash variable is empty, let's use an example data
fi
#OUTPUT_DIR="${HOME}/qc-results/${SAMPLE_ID}"
OUTPUT_DIR="${HOME}/bioinfo-results/${SAMPLE_ID}/assembly_de_novo"
[ -d ${OUTPUT_DIR} ] || mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

case $1 in
    ##
    # Illumina data
    ##
    "-spades")
        # Decompress input files from quality-filter dir
        # gzip -d ${INPUT_DIR}/*.gz
        # 1) Use of spades
        # Link: https://github.com/ablab/spades/blob/spades_3.15.5/README.md
        # They recommend running SPAdes with BayesHammer/IonHammer to obtain high-quality assemblies.
        # Note: We decided to avoid unpaired reads!!!!
        # For single lib use -1 and -2
        # spades -1 ${INPUT_DIR}/output_forward_paired.fq.gz -2 ${INPUT_DIR}/output_reverse_paired.fq.gz -o ${OUTPUT_DIR}
        # For single lib all reads paired and unpaired use -s
        spades -1 ${INPUT_DIR}/output_forward_paired.fq.gz -2 ${INPUT_DIR}/output_reverse_paired.fq.gz -s ${INPUT_DIR}/output_forward_unpaired.fq.gz -s ${INPUT_DIR}/output_reverse_unpaired.fq.gz -o ${OUTPUT_DIR}
    ;;
    "-velvet")
        # 2) Use of velvet
        # Link: https://github.com/dzerbino/velvet
        velveth ${OUTPUT_DIR} 21 -fastq -short ${INPUT_DIR}/output_forward_paired.fq.gz ${INPUT_DIR}/output_reverse_paired.fq.gz ${INPUT_DIR}/output_forward_unpaired.fq.gz ${INPUT_DIR}/output_reverse_unpaired.fq.gz
        velvetg ${OUTPUT_DIR} -cov_cutoff 4 -min_contig_lgth 100
    
    ;;
    ##
    # MinIon data
    ##
    "-minimap")
        # 1) Use of minimap
        # Link: https://timkahlke.github.io/LongRead_tutorials/ASS_M.html
        echo "Procedure under development..."
        exit 0;  
    ;;
    *)
        echo "Invalid parameter!"
        exit 1
esac

