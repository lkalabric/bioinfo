library(fastqcr)

SAMPLE = '0001.1'
INPUT_DIR = paste('data/hbv/',SAMPLE)
OUTPUT_DIR = paste('qc-results/',SAMPLE)

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
