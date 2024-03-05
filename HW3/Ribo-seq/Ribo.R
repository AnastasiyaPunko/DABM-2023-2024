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
library(corrplot)
library(RColorBrewer)

theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())


#####--------------- Load Database ---------------#####
df<-readr::read_tsv("01. RiboSeq_RNASeq_HCC_counts.tsv")
df <- as.data.frame(df)

df[is.na(df)] <- 0
#####---------------Correlation coefficient of the number of counts between RNA-Seq and Ribo-Seq experiments for each gene---------------#####
df <-  df %>% dplyr::select(geneSymbol, 
                      `LC001-normal-RPF`, `LC001-tumor-RPF`, `LC033-normal-RPF`, `LC033-tumor-RPF`, #Ribo
                      `LC034-normal-RPF`, `LC034-tumor-RPF`, `LC501-normal-RPF`, `LC501-tumor-RPF`,
                      `LC502-normal-RPF`, `LC502-tumor-RPF`, `LC505-normal-RPF`, `LC505-tumor-RPF`,
                      `LC506-normal-RPF`, `LC506-tumor-RPF`, `LC507-normal-RPF`, `LC507-tumor-RPF`, 
                      `LC508-normal-RPF`, `LC508-tumor-RPF`, `LC509-normal-RPF`, `LC509-tumor-RPF`,
                      `LC001-normal-RNA`, `LC001-tumor-RNA`, `LC033-normal-RNA`, `LC033-tumor-RNA`, #RNA
                      `LC034-normal-RNA`, `LC034-tumor-RNA`, `LC501-normal-RNA`, `LC501-tumor-RNA`,
                      `LC502-normal-RNA`, `LC502-tumor-RNA`, `LC505-normal-RNA`, `LC505-tumor-RNA`,
                      `LC506-normal-RNA`, `LC506-tumor-RNA`, `LC507-normal-RNA`, `LC507-tumor-RNA`, 
                      `LC508-normal-RNA`, `LC508-tumor-RNA`, `LC509-normal-RNA`, `LC509-tumor-RNA`)

test <-  df[3,2:40]
hist(as.numeric(test))
shapiro.test(as.numeric(test)) #If the value of p <= 0.05, then the hypothesis of normality will be rejected 

plot(df$`LC502-normal-RPF`, df$`LC502-normal-RNA`)
ggplot(df, aes(df$`LC502-normal-RPF`, df$`LC502-normal-RNA`)) + geom_text(aes(label=df$geneSymbol), hjust = 0, nudge_x = 3, size = 3) #using the DESEQ2 plotPCA fxn we can
df2 <- df[!(df$geneSymbol %in% c("ALB")),]

#Pearson correlation
corr_pearson <- data.frame(t(apply(df, 1, function(x) c(cor(as.numeric(x[2:21]),as.numeric(x[22:41])), method = 'pearson'))))

corr_pearson <- corr_pearson %>% select(1)
corr_pearson$Gene <- df$geneSymbol
colnames(corr_pearson)[1] <- 'Correlation'

corr_pearson <- na.omit(corr_pearson)
corr_pearson$Correlation <- as.numeric(corr_pearson$Correlation)
max(corr_pearson$Correlation) #1
min(corr_pearson$Correlation) #-0.4997039
mean(corr_pearson$Correlation) #0.5581791

write.csv(corr_pearson, "corr_pearson.csv", row.names=TRUE)

#Distribution histogram
hist(corr_pearson$Correlation)

#Top genes
corr_pearson <- corr_pearson[order(corr_pearson$Correlation),]
tail(corr_pearson, n=10L)
head(corr_pearson, n=10L)


#Spearman correlation
corr_spearman <- data.frame(t(apply(df, 1, function(x) c(cor(as.numeric(x[2:21]),as.numeric(x[22:41])), method = 'spearman'))))

corr_spearman <- corr_spearman %>% select(1)
corr_spearman$Gene <- df$geneSymbol
colnames(corr_spearman)[1] <- 'Correlation'

corr_spearman <- na.omit(corr_spearman)
corr_spearman$Correlation <- as.numeric(corr_spearman$Correlation)
max(corr_spearman$Correlation) #1
min(corr_spearman$Correlation) #-0.4997039
mean(corr_spearman$Correlation) #0.5581791

write.csv(corr_spearman, "corr_spearman.csv", row.names=TRUE)

#Distribution histogram
hist(corr_spearman$Correlation)

#Top genes
corr_spearman <- corr_spearman[order(corr_spearman$Correlation),]
tail(corr_spearman, n=10L)
head(corr_spearman, n=10L)

#####---------------Correlation coefficient in normal sample & tumor sample---------------#####
#norm
df_norm <-  df %>% dplyr::select(geneSymbol, 
                      `LC001-normal-RPF`, `LC033-normal-RPF`, 
                      `LC034-normal-RPF`, `LC501-normal-RPF`,
                      `LC502-normal-RPF`, `LC505-normal-RPF`,
                      `LC506-normal-RPF`, `LC507-normal-RPF`,
                      `LC508-normal-RPF`, `LC509-normal-RPF`,
                      `LC001-normal-RNA`, `LC033-normal-RNA`, 
                      `LC034-normal-RNA`, `LC501-normal-RNA`,
                      `LC502-normal-RNA`, `LC505-normal-RNA`,
                      `LC506-normal-RNA`, `LC507-normal-RNA`, 
                      `LC508-normal-RNA`, `LC509-normal-RNA`)

colnames(df_norm[2:11])
colnames(df_norm[12:21])

#Spearman correlation
corr_spearman <- data.frame(t(apply(df_norm, 1, function(x) c(cor(as.numeric(x[2:21]),as.numeric(x[22:41])), method = 'spearman'))))

corr_spearman <- corr_spearman %>% select(1)
corr_spearman$Gene <- df_norm$geneSymbol
colnames(corr_spearman)[1] <- 'Correlation'

corr_spearman <- na.omit(corr_spearman)
corr_spearman$Correlation <- as.numeric(corr_spearman$Correlation)
max(corr_spearman$Correlation) #1
min(corr_spearman$Correlation) #-0.4997039
mean(corr_spearman$Correlation) #0.5811202

#Distribution histogram
hist(corr_spearman$Correlation)

#Top genes
corr_spearman <- corr_spearman[order(corr_spearman$Correlation),]
tail(corr_spearman, n=10L)
head(corr_spearman, n=10L)

#tumor
df_t <-  df %>% dplyr::select(geneSymbol, 
                          `LC001-tumor-RPF`, `LC033-tumor-RPF`, 
                          `LC034-tumor-RPF`, `LC501-tumor-RPF`,
                          `LC502-tumor-RPF`, `LC505-tumor-RPF`,
                          `LC506-tumor-RPF`, `LC507-tumor-RPF`,
                          `LC508-tumor-RPF`, `LC509-tumor-RPF`,
                          `LC001-tumor-RNA`, `LC033-tumor-RNA`, 
                          `LC034-tumor-RNA`, `LC501-tumor-RNA`,
                          `LC502-tumor-RNA`, `LC505-tumor-RNA`,
                          `LC506-tumor-RNA`, `LC507-tumor-RNA`, 
                          `LC508-tumor-RNA`, `LC509-tumor-RNA`)

#Spearman correlation
corr_spearman <- data.frame(t(apply(df_t, 1, function(x) c(cor(as.numeric(x[2:21]),as.numeric(x[22:41])), method = 'spearman'))))

corr_spearman <- corr_spearman %>% select(1)
corr_spearman$Gene <- df_t$geneSymbol
colnames(corr_spearman)[1] <- 'Correlation'

corr_spearman <- na.omit(corr_spearman)
corr_spearman$Correlation <- as.numeric(corr_spearman$Correlation)
max(corr_spearman$Correlation) #1
min(corr_spearman$Correlation) #-0.4997039
mean(corr_spearman$Correlation) #0.5811202

#Distribution histogram
hist(corr_spearman$Correlation)

#Top genes
corr_spearman <- corr_spearman[order(corr_spearman$Correlation),]
tail(corr_spearman, n=10L)
head(corr_spearman, n=10L)


#####--------------- Common correlation ---------------#####
df[is.na(df)] <- 0

M <-cor(df[,unlist(lapply(df, is.numeric))])
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

RPF_norm_tumor <-c(df$`LC001-normal-RPF`, df$`LC001-tumor-RPF`, df$`LC033-normal-RPF`, df$`LC033-tumor-RPF`, 
                   df$`LC034-normal-RPF`, df$`LC034-tumor-RPF`, df$`LC501-normal-RPF`, df$`LC501-tumor-RPF`,
                   df$`LC502-normal-RPF`, df$`LC502-tumor-RPF`, df$`LC505-normal-RPF`, df$`LC505-tumor-RPF`,
                   df$`LC506-normal-RPF`, df$`LC506-tumor-RPF`, df$`LC507-normal-RPF`, df$`LC507-tumor-RPF`, 
                   df$`LC508-normal-RPF`, df$`LC508-tumor-RPF`, df$`LC509-normal-RPF`, df$`LC509-tumor-RPF`)

RNA_norm_tumor <-c(df$`LC001-normal-RNA`, df$`LC001-tumor-RNA`, df$`LC033-normal-RNA`, df$`LC033-tumor-RNA`, 
                   df$`LC034-normal-RNA`, df$`LC034-tumor-RNA`, df$`LC501-normal-RNA`, df$`LC501-tumor-RNA`,
                   df$`LC502-normal-RNA`, df$`LC502-tumor-RNA`, df$`LC505-normal-RNA`, df$`LC505-tumor-RNA`,
                   df$`LC506-normal-RNA`, df$`LC506-tumor-RNA`, df$`LC507-normal-RNA`, df$`LC507-tumor-RNA`, 
                   df$`LC508-normal-RNA`, df$`LC508-tumor-RNA`, df$`LC509-normal-RNA`, df$`LC509-tumor-RNA`)

cor(RPF_norm_tumor, RNA_norm_tumor, method = c("pearson")) # 0.7979086

cor(RPF_norm_tumor, RNA_norm_tumor, method = c("spearman")) # 0.859714
plot(x = RPF_norm_tumor, y = RNA_norm_tumor, xlab = 'Ribo-seq', ylab = 'RNA-seq', col='#0c4c8a', 
     main = 'Корреляция числа каунтов между экспериментами RNA-Seq и Ribo-Seq в нормальной ткани и в опухолевой')

RPF_norm <- c(df$`LC001-normal-RPF`, df$`LC033-normal-RPF`, 
              df$`LC034-normal-RPF`, df$`LC501-normal-RPF`,
              df$`LC502-normal-RPF`, df$`LC505-normal-RPF`,
              df$`LC506-normal-RPF`, df$`LC507-normal-RPF`, 
              df$`LC508-normal-RPF`, df$`LC509-normal-RPF`)

RPF_tumor <- c(df$`LC001-tumor-RPF`, df$`LC033-tumor-RPF`, 
               df$`LC034-tumor-RPF`, df$`LC501-tumor-RPF`,
               df$`LC502-tumor-RPF`, df$`LC505-tumor-RPF`,
               df$`LC506-tumor-RPF`, df$`LC507-tumor-RPF`, 
               df$`LC508-tumor-RPF`, df$`LC509-tumor-RPF`)

RNA_norm <- c(df$`LC001-normal-RNA`, df$`LC033-normal-RNA`, 
              df$`LC034-normal-RNA`, df$`LC501-normal-RNA`,
              df$`LC502-normal-RNA`, df$`LC505-normal-RNA`,
              df$`LC506-normal-RNA`, df$`LC507-normal-RNA`, 
              df$`LC508-normal-RNA`, df$`LC509-normal-RNA`)

RNA_tumor <-c(df$`LC001-tumor-RNA`, df$`LC033-tumor-RNA`, 
              df$`LC034-tumor-RNA`, df$`LC501-tumor-RNA`,
              df$`LC502-tumor-RNA`, df$`LC505-tumor-RNA`,
              df$`LC506-tumor-RNA`, df$`LC507-tumor-RNA`, 
              df$`LC508-tumor-RNA`, df$`LC509-tumor-RNA`)

cor(RPF_norm, RNA_norm, method = c("pearson")) # 0.9008479
cor(RPF_norm, RNA_norm, method = c("spearman")) # 0.8357499
plot(x = RPF_norm, y = RNA_norm, xlab = 'Ribo-seq', ylab = 'RNA-seq', col='#0c4c8a', 
     main = 'Корреляция числа каунтов между экспериментами RNA-Seq и Ribo-Seq в нормальной ткани')

cor(RPF_tumor, RNA_tumor, method = c("pearson")) # 0.6774872
cor(RPF_tumor, RNA_tumor, method = c("spearman")) # 0.8806405
plot(x = RPF_tumor, y = RNA_tumor, xlab = 'Ribo-seq', ylab = 'RNA-seq', col='#0c4c8a', 
     main = 'Корреляция числа каунтов между экспериментами RNA-Seq и Ribo-Seq в опухолевой ткани')

