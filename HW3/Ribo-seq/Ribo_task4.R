## Author: Anastasiya Punko
## First created: 27-11-2023
## Description: HW3 NGS
## Keywords: Ribo-seq, RNA-seq
rm(list = ls())

getwd()

#####--------------- Load libraries and functions ---------------#####
library(readr)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(tximport)
library(ggrepel)


theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())


#####--------------- Load Database ---------------#####
df<-readr::read_tsv("01. RiboSeq_RNASeq_HCC_counts.tsv")
df <- as.data.frame(df)

#####--------------- DE genes RNA-seq---------------#####
df_RNA <-  df %>% dplyr::select(geneID, geneSymbol,
                     `LC001-normal-RNA`, `LC033-normal-RNA`, `LC034-normal-RNA`, `LC501-normal-RNA`, `LC502-normal-RNA`, 
                     `LC505-normal-RNA`, `LC506-normal-RNA`, `LC507-normal-RNA`, `LC508-normal-RNA`, `LC509-normal-RNA`,
                     `LC001-tumor-RNA`, `LC033-tumor-RNA`, `LC034-tumor-RNA`, `LC501-tumor-RNA`, `LC502-tumor-RNA`,  
                     `LC505-tumor-RNA`, `LC506-tumor-RNA`, `LC507-tumor-RNA`, `LC508-tumor-RNA`, `LC509-tumor-RNA`)

df_RNA[is.na(df_RNA)] <- 0
rownames(df_RNA) = make.names(df_RNA$geneSymbol, unique=TRUE)
## Create a sampletable/metadata
sampletype <- factor(c(rep("normal",10), rep("tumor", 10)))
meta <- data.frame(sampletype, row.names = colnames(df_RNA[3:22]))

#-----------------------DeSeq--------------------------------####
## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = df_RNA[3:22], colData = meta, design = ~ sampletype)
dds <- DESeq(dds)

counts <- counts(dds)

#Take a look at the results table
res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table


#Summary of differential gene expression
summary(res) #summary of results
plotMA(res, ylim=c(-2,2)) #The genes that are significantly DE are colored 

#Sort summary list by p-value
res <- res[order(res$padj),]
head(res)


#Volcano Plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-10,10)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))

result <- as.data.frame(res) %>% filter(res$padj < 0.05 & res$pvalue < 0.05 & res$log2FoldChange>1)
result %>% view()
write.csv(result, "result_RNA.csv", row.names=TRUE)

#PCA
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="sampletype") + coord_fixed(ylim =c(-30,30)) + geom_text(aes(label=name), hjust = 0, nudge_x = 3, size = 3) #using the DESEQ2 plotPCA fxn we can


#####--------------- DE genes Ribo-seq---------------#####
df_RPF <-  df %>% dplyr::select(geneID, geneSymbol,
                         `LC001-normal-RPF`, `LC033-normal-RPF`, `LC034-normal-RPF`, `LC501-normal-RPF`, `LC502-normal-RPF`, 
                         `LC505-normal-RPF`, `LC506-normal-RPF`, `LC507-normal-RPF`, `LC508-normal-RPF`, `LC509-normal-RPF`,
                         `LC001-tumor-RPF`, `LC033-tumor-RPF`, `LC034-tumor-RPF`, `LC501-tumor-RPF`, `LC502-tumor-RPF`,  
                         `LC505-tumor-RPF`, `LC506-tumor-RPF`, `LC507-tumor-RPF`, `LC508-tumor-RPF`, `LC509-tumor-RPF`)

df_RPF[is.na(df_RPF)] <- 0
rownames(df_RPF) = make.names(df_RPF$geneSymbol, unique=TRUE)
## Create a sampletable/metadata
sampletype <- factor(c(rep("normal",10), rep("tumor", 10)))
meta <- data.frame(sampletype, row.names = colnames(df_RPF[3:22]))

#-----------------------DeSeq--------------------------------####
## Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = df_RPF[3:22], colData = meta, design = ~ sampletype)
dds <- DESeq(dds)

counts <- counts(dds)

#Take a look at the results table
res <- results(dds)
head(results(dds, tidy=TRUE)) #let's look at the results table

#Summary of differential gene expression
summary(res) #summary of results
plotMA(res, ylim=c(-2,2)) #The genes that are significantly DE are colored 

#Sort summary list by p-value
res <- res[order(res$padj),]
head(res)


#Volcano Plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-10,10)))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="red"))

result <- as.data.frame(res) %>% filter(res$padj < 0.05 & res$pvalue < 0.05 & res$log2FoldChange>1)
result %>% view()
write.csv(result, "result_Ribo.csv", row.names=TRUE)

#PCA
#First we need to transform the raw count data
#vst function will perform variance stabilizing transformation
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="sampletype") + coord_fixed(ylim =c(-30,30)) + geom_text(aes(label=name), hjust = 0, nudge_x = 3, size = 3) #using the DESEQ2 plotPCA fxn we can



