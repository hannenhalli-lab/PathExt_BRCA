### Link where I learned 	https://lashlock.github.io/compbio/R_presentation.html #######

library(DESeq2)

## read the raw count data as an input. make sure the values are in integer and no duplicate genes should be there. Fill the NAs with 0 if required ###
countData <- read.csv("DESeq2_input", header = TRUE, sep = ",")

###### prepare a metafile, where first column should be sample id, 2nd column shoud be condition (Normal, Tumor) and 3 column is about the tissue ####
metaData <- read.csv("DESeq2_metadata", header = TRUE, sep = ",")

#### Construct DESEQDataSet Object####### 
dds <- DESeqDataSetFromMatrix(countData=countData, colData=metaData, design=~dex, tidy = TRUE)

#### Now we’re ready to run DESEQ function #####
dds <- DESeq(dds)

#### Take a look at the results table ######
res <- results(dds)
head(res)

### Write the output ########
write.table(results(dds), file = "DESeq2_Out", sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)
