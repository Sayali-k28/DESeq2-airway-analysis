# Script 1: extract the count matrix and sample information from the airway package

# install packages if needed
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("airway", quietly = TRUE)) BiocManager::install("airway")
if (!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")

library(airway)

# load the dataset (a SummarizedExperiment object)
data(airway)
airway

# sample information: keep cell line and treatment columns
sample_info <- as.data.frame(colData(airway))
sample_info <- sample_info[, c(2, 3)]

# use clearer labels for the treatment groups
sample_info$dex <- gsub("trt", "treated", sample_info$dex)
sample_info$dex <- gsub("untrt", "untreated", sample_info$dex)
names(sample_info) <- c("cellLine", "dexamethasone")
write.table(sample_info, file = "sample_info.csv", sep = ",",
            col.names = TRUE, row.names = TRUE, quote = FALSE)

# raw count matrix (genes x samples)
countsData <- assay(airway)
write.table(countsData, file = "counts_data.csv", sep = ",",
            col.names = TRUE, row.names = TRUE, quote = FALSE)
