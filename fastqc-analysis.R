# Projeto: Bioinfo
# Autor: Luciano Kalabric
# Data: 13/06/2023

# Installation from CRAN
install.packages("fastqcr")

# Or, install the latest version from GitHub:
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/fastqcr")

# Main
library(fastqcr)

# Environment preparation
SAMPLE = '0001.1'


# Set working directory in Linux
setwd("~/") # You only need to do this once
INPUT_DIR = paste('data/hbv/',SAMPLE,sep="")
OUTPUT_DIR = paste('qc-results/',SAMPLE,sep="")

# Set working directory in Windows
#setwd("C:/Users/kalabric/Downloads") # You only need to do this once
#INPUT_DIR = paste(getwd(),SAMPLE,sep="/")
#OUTPUT_DIR = paste(getwd(),'qc-results',SAMPLE,sep="/")


# Run FastQC to generation FastQC Reports 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fastqc(fq.dir = INPUT_DIR, # FASTQ files directory
       qc.dir = OUTPUT_DIR, # Results direcory
       threads = 4                    # Number of threads
)

# Aggregating Multiple FastQC Reports into a Data Frame 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Demo QC directory containing zipped FASTQC reports
qc.dir <- system.file(INPUT_DIR, package = "fastqcr")
qc <- qc_aggregate(qc.dir)
qc

# Inspecting QC Problems
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# See which modules failed in the most samples
qc_fails(qc, "module")
# Or, see which samples failed the most
qc_fails(qc, "sample")

# Building Multi QC Reports
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qc_report(qc.dir, result.file = "multi-qc-report" )

# Building One-Sample QC Reports (+ Interpretation)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qc.file <- system.file(OUTPUT_DIR, "*R1*", package = "fastqcr")
qc_report(qc.file, result.file = "one-sample-report",
          interpret = TRUE)
