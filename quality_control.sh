# There are different alternatives methods in bioinfo to do the quality control. Let's see them both: 

# 1) Fastqc
# Link: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# Requirements: Java
# Installation:
# $ sudo apt install default-jre

# 2) Fastqrc 
# Link: https://rpkgs.datanovia.com/fastqcr/index.html
# Requirements: R

# Using example data
fastqc -o qc-results -f data/hbv/0001.1/*.fastq

# 3) Software favorito de Gilmar
# Link:
# Requirements:

