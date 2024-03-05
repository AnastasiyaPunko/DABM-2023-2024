if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("tximportData")


library(tximport)
library(tximportData)
library(DESeq2)

dir <- system.file("extdata", package = "tximportData")
vignette("tximportData")
library(readr)
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))


# Загрузка данных с использованием tximport
abundance <-  read_tsv("abundance.tsv")
abundance2 <-  read_tsv("abundance-2.tsv")
files <- c(abundance, abundance2 )  # замените путями к файлам abundance.tsv
tx2gene <- read.table("tx2gene.gencode.v27.csv", header=TRUE, stringsAsFactors=FALSE)
txi <- tximport(files, type="kallisto", tx2gene=tx2gene)

# Создание объекта DESeqDataSet
dds <- DESeqDataSetFromTximport(txi, colData=coldata, design=~group)


2. Провести дифференциальную экспрессию между образцами из групп A и B на уровне генов с использованием DESeq2:
  
  # Проведение дифференциальной экспрессии
  dds <- DESeq(dds)
res <- results(dds)





samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
names(files) <- paste0("sample",1:6)
# tx2gene links transcript IDs to gene IDs for summarization
tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v27.csv"))
txi <- tximport(files, type="salmon", tx2gene=tx2gene)





#-----------------####

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("DESeq2")

library(tximport)
library(DESeq2)


# Загрузите файлы abundance.tsv для групп A и B
abundance_A <- read.table("abundance.tsv", header = TRUE, row.names = 1, sep = "\t")
abundance_B <- read.table("abundance-2.tsv", header = TRUE, row.names = 1, sep = "\t")

# Создайте объект tx2gene, указывая соответствие транскриптов генам
tx2gene <- read.table(file.path(dir,"tx2gene.gencode.v27.csv", header = TRUE, sep = ","))

txi_B <- tximport(files, type = "kallisto", tx2gene = tx2gene)

BiocManager::install("tximportData")
library(tximportData)
dir <- system.file("extdata", package = "tximportData")
vignette("tximportData")

library(readr)
tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv", header = TRUE, sep = ","))


dir2 <- file.path('Users', 'anastasiyapunko', fsep="/")
files <-  
  
