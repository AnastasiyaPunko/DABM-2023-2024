
rm(list=ls())
getwd()
## Gene-level differential expression analysis using DESeq2

#-----------------------Setup--------------------------------####
### Bioconductor and CRAN libraries used
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
install.packages('pheatmap')
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)


#-----------------------tximport--------------------------------####
## List all directories containing data  
samples <- list.files(path = "./data", full.names = T, pattern="kallisto$")
## Obtain a vector of all filenames including the path
files <- file.path(samples, "abundance.tsv")
## Since all quant files have the same name it is useful to have names for each element
names(files) <- str_replace(samples, "./data/", "") %>% 
  str_replace(".kallisto", "")

# Load the annotation table for GrCh38
tx2gene <- read.delim("tx2gene_grch38_ens94.txt")

# Take a look at it 
tx2gene %>% View()

# Run tximport
txi <- tximport(files, type="kallisto", tx2gene=tx2gene[,c("tx_id", "ensgene")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE)

head(txi$counts)
attributes(txi)

# Look at the counts
txi$counts %>% View()

# Write the counts to an object
data <- txi$counts %>% 
  round() %>% 
  data.frame()

## Create a sampletable/metadata
sampletype <- factor(c(rep("A",3), rep("B", 4), rep("A", 1)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

path_to_save <- "./meta"
write.table(meta, file = "meta.txt")    
#-----------------------DeSeq--------------------------------####
## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)
dds <- DESeq(dds)
plotDispEsts(dds)

path_to_save <- "./results"

View(counts(dds))
counts <- counts(dds)
write.csv(counts, "counts.csv", row.names=TRUE)

#Take a look at the results table
res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table

write.csv(res, "res.csv", row.names=TRUE)

result <- as.data.frame(res) %>% filter(res$padj < 0.05 & res$pvalue < 0.05 & res$log2FoldChange>1)
result %>% view()
#Summary of differential gene expression
summary(res) #summary of results
plotMA(res, ylim=c(-2,2)) #The genes that are significantly DE are colored 

#Sort summary list by p-value
res <- res[order(res$padj),]
head(res)

#plotCounts
#we can use plotCounts fxn to compare the normalized counts
#between treated and control groups for our top 6 genes
par(mfrow=c(2,3))

plotCounts(dds, gene="ENSG00000115850", intgroup="sampletype")
plotCounts(dds, gene="ENSG00000127831", intgroup="sampletype")
plotCounts(dds, gene="ENSG00000130234", intgroup="sampletype")
plotCounts(dds, gene="ENSG00000134028", intgroup="sampletype")
plotCounts(dds, gene="ENSG00000138792", intgroup="sampletype")
plotCounts(dds, gene="ENSG00000158578", intgroup="sampletype")


#Volcano Plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-15,15)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))

#PCA
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="sampletype") + coord_fixed(ylim =c(-30,30)) + geom_text(aes(label=name), hjust = 0, nudge_x = 3, size = 3) #using the DESEQ2 plotPCA fxn we can



#-----------------------Functional analysis--------------------------------####
# Create background dataset for hypergeometric testing using all genes tested for significance in the results
all_genes <- as.character(rownames(res))
# Extract significant results
signif_res <- res[res$padj < 0.05 & res$pvalue < 0.05 & res$log2FoldChange > 1 & !is.na(res$padj), ] 
signif_genes <- as.character(rownames(signif_res))
signif_genes %>% view()

write.csv(signif_genes, "signif_genes.csv", row.names=TRUE)
