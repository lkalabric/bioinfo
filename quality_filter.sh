#!/bin/bash

# command name: quality_filter.sh
# author: Luciano Kalabric Silva
# institution: Oswaldo Cruz Foundation, Goncalo Moniz Institute, Bahia, Brazil
# last update: 15 JUN 2023
# objetive: Give examples of quality filter apps
# Syntax: ./quality_filter.sh

# Variable definitions
INPUTDIR="${HOME}/data/hbv/$1"
OUTPUTDIR="${HOME}/qc-filter/hbv/$1"
[ ! -d ${OUTPUTDIR} ] && mkdir -p ${OUTPUTDIR}
MINLENTGH=100

# Validate arguments
if ["$#" -ne 2]; then
    echo "Illegal number of parameters"
    echo "Syntax: quality_filter.sh <SAMPLE_NAME> <-illumina | -minion>"
fi

case $2 in
  "-illumina")
      ##
    # Illumina data
    ##
    # 1) Afterqc
    # Link: https://github.com/OpenGene/AfterQC
    # Requirements: Miniconda (Python)
    # Installation:
    # $ install_thirdparty.sh
    source activate afterqc
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
    trimmomatic PE ${INPUTDIR}/*R1*.fastq.gz ${INPUTDIR}/*R2*.fastq.gz ${OUTPUTDIR}/output_forward_paired.fq.gz ${OUTPUTDIR}/output_forward_unpaired.fq.gz ${OUTPUTDIR}/output_reverse_paired.fq.gz ${OUTPUTDIR}/output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:${MINLENTGH}
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
  ;;
  *)
  echo "Invalid parameter!"
  exit 0
esac
  







