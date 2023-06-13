# There are different alternatives methods in bioinfo to do the quality control. Let's see some of them: 

# 1) Fastqc
# Link: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# Requirements: Java
# $ sudo apt install default-jre
# Installation:
# $ sudo apt install fastqc

# Using example data
fastqc -o qc-results -f fastq -c data/hbv/0001.1/*.fastq

# 2) Fastqrc 
# Link: https://rpkgs.datanovia.com/fastqcr/index.html
# Requirements: R


# 3) Afterqc
# Link: https://github.com/OpenGene/AfterQC
# Requirements: Miniconda
source activate afterqc
after.py -1 data/hbv/0001.1/A24_S6_L001_R1_001.fastq -2 data/hbv/0001.1/A24_S6_L001_R2_001.fastq
