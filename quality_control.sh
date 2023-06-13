# There are different alternatives methods in bioinfo to do the quality control. Let's see them both: 

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

