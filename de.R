#!/usr/bin/env Rscript
# de.R
library(knitr)
library(tximport)
library(readr)
library(DESeq2)

# TODO: update constants for your machine
# Define constants
TESTING <- FALSE# Change to FALSE if using entire Samples set
RESULTS_DIR <- "/home/dropathi.m/BINF6309/module04-meghana-dropathi"
AIPTASIA_DIR <- "/work/courses/BINF6309/AiptasiaMiSeq"

# for testing purposes - alternative samples table
testing_samples <- data.frame(Sample = c("Aip02", "Aip02", "Aip02", "Aip02"),
                              Menthol = c("Control", "Control", "Menthol", "Menthol"),
                              Vibrio = c("Control", "Vibrio", "Control", "Vibrio"))
head(testing_samples)

# True script begins
tx2gene <- read.csv(file.path(RESULTS_DIR, "tx2gene.csv"))
head(tx2gene)

if (TESTING) {
  print("***Running test with Aip02 only***")
  samples <- testing_samples
} else {
  samples <- read.csv(file.path(AIPTASIA_DIR, "Samples.csv"), header=TRUE)
}
head(samples)


files <- file.path(RESULTS_DIR, "quant", samples$Sample, "quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

dds <- DESeqDataSetFromTximport(txi, colData = samples, 
                                design = ~ Menthol + Vibrio)

dds$Vibrio <- relevel(dds$Vibrio, ref = "Control")
dds$Menthol <- relevel(dds$Menthol, ref = "Control")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)

padj <- .05
minLog2FoldChange <- .5
dfAll <- data.frame()
# Get all DE results except Intercept, and "flatten" into a single file.
for (result in resultsNames(dds)){
  if(result != 'Intercept'){
    res <- results(dds, alpha=.05, name=result)
    dfRes <- as.data.frame(res)
    dfRes <- subset(subset(dfRes, select=c(log2FoldChange, padj)))
    dfRes$Factor <- result
    dfRes$ko <- rownames(dfRes)
    dfAll <- rbind(dfAll, dfRes)
  }
}
rownames(dfAll) <- NULL
head(dfAll)

write.csv(dfAll, file=file.path(RESULTS_DIR, "dfAll.csv"))
# end of de.R script
# filtering according to padj values
padj_filter<- subset(dfAll, padj<.05)

#reading pathways as table and assigning column names
path_ways <- read.table("/work/courses/BINF6309/data_BINF6309/Module4/Annotation/path.txt", sep='\t', header=FALSE)
colnames(path_ways) <- c("ko", "pathway")

#reading ko as table and assigning colum names 
pathway_names <- read.table( "/work/courses/BINF6309/data_BINF6309/Module4/Annotation/ko", sep="\t", header=FALSE)
colnames(pathway_names) <- c("pathway", "description")

# Merging the required fields
merged_path_ko<- merge(path_ways, pathway_names)
# merge dfAll with merged_pko
deAnnotated <- merge(merged_path_ko, padj_filter)

# displaying the values as kable formatted table
kable(deAnnotated)
# writing as a file 
write.csv(deAnnotated, file=file.path(RESULTS_DIR,"deAnnotated.csv"),row.names=FALSE)