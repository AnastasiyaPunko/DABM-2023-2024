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

#####--------------- Ribo-seq count distribution analysis ---------------#####
df <-  df %>% dplyr::select(`LC001-normal-RPF`, `LC001-tumor-RPF`, `LC033-normal-RPF`, `LC033-tumor-RPF`,
                     `LC034-normal-RPF`, `LC034-tumor-RPF`, `LC501-normal-RPF`, `LC501-tumor-RPF`,
                     `LC502-normal-RPF`, `LC502-tumor-RPF`, `LC505-normal-RPF`, `LC505-tumor-RPF`,
                     `LC506-normal-RPF`, `LC506-tumor-RPF`, `LC507-normal-RPF`, `LC507-tumor-RPF`, 
                     `LC508-normal-RPF`, `LC508-tumor-RPF`, `LC509-normal-RPF`, `LC509-tumor-RPF`)


mean_counts <- apply(df[,1:20], 1, mean)        #The second argument '1' of 'apply' function indicates the function being applied to rows. Use '2' if applied to columns 
variance_counts <- apply(df[,1:20], 1, var)
df_stat <- data.frame(mean_counts, variance_counts)

ggplot(df_stat) +                                    #Each data point represents a gene and the red line represents x = y.
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

#####--------------- Distribution ---------------#####

install.packages('fitdistrplus')
library(fitdistrplus)

test <-  df[3,1:20]
hist(as.numeric(test))

data <- as.vector(as.numeric(test))
fit.norm <- fitdist(data, "norm")
summary(fit.norm)
plot(fit.norm)

fit_n  <- fitdist(data, "norm")
fit_p  <- fitdist(data, "pois")
fit_nb  <- fitdist(data, "nbinom")
fit_ln  <- fitdist(data, "lnorm")
cdfcomp(list(fit_n, fit_p, fit_ln, fit_nb), legendtext = c("Normal", "Poisson", "Log Normal", "Negative Binomial"))

summary(fit_nb)
plot(fit_nb)
summary(fit_n)

fit_ln.boot <- bootdist(fit_ln, "param")
CIcdfplot(fit_ln.boot, CI.output = "probability")

descdist(data, discrete=TRUE, boot=500)
