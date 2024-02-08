# Author 

**Meghana Reddy Dropathi**

# Differential Expression Analysis

# Contents
*This github repository contains five programs:*
1. de.R
2. deAnnotated.csv
3. methodsResults.html
5. methodsResults.Rmd

# Result files

1.deAnnotated.csv

2.methodsResults.html

**The deAnnotated.csv has the adjusted p-values where padj is less than 0.05 with pathway and pathway names merged with results.**

**methodsResults.html is created by rendering methodsResults.Rmd to html by knitting.**

# overview

In this project, I did differential expression analysis. We used Salmon(Patro et al. 2017) to find all the Aip samples in the “/work/courses/BINF6309/AiptasiaMiSeq/fastq/ directory” and aligned them to AipIndex.In order to use the Salmon input in tximport(Soneson, Love, and Robinson 2016) I created a table, mapping transcripts to genes. I used annotation files in “/work/courses/BINF6309/data_BINF6309/Module4/Annotation/” and I showed the table using “kable”(part of knitr library). This creates a tx2gene.csv file which can be used further in analysis.I create a script called de.R which uses DESeq2(Love, Huber, and Anders 2014) and tximport(Soneson, Love, and Robinson 2016) to get the output. The script uses the DESeqDataSetFromTximport function to create a DESeqDataSet object from the tximport output. It specifies the experimental design with “Menthol” and “Vibrio” as factors. Rows with counts less than 10 are removed to get better efficiency in the results as it removes low count genes. After this DESeq analysis is run.The dfAll.csv file is created. I filter the values of dfAll with padj less than 0.05 and we merge pathways,pathway names with results and write it into deAnnotated.csv. This csv file consists differentially expressed gene information with padj values less than 0.05.
