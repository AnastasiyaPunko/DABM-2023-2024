## Author: Anastasiya Punko
## First created: 08-10-2023
## Description: HW1
## Keywords: PK modeling, Two-Compartment Model, Human Ethanol Metabolism
rm(list = ls())

#####--------------- Load libraries and functions ---------------#####
install.packages("tidyverse")
install.packages("readxl")
install.packages("gridExtra")
library(tidyverse)
library(readxl)
library(gridExtra)
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())


#####--------------- Parameters used in the model  ---------------#####
FHV <- 1.5 #liters/min
VB <- 48 #liters
VL <- 0.61 #liters
Vmax <- 2.75 #mmol/min/80 kg 
Km <- 0.1 #mM 
CL = 0
CB = 0
QIV = 720 * 0.08 / (2 * 60) * 0.789 / 46.069

N = 1000000

t <- seq(from = 0, to = 480,length.out=1000000)

CB <- rep(0, length(t))
CL <- rep(0, length(t))

range <- 2:N
CB <- c()
for (i in range) {
  CB[i] = CB[i - 1] + (t[i] - t[i - 1]) * (FHV * (CL[i - 1] - CB[i - 1]) + QIV) / VB
}
lapply(1:1000000,function(i) {
  CB[i] = CB[i - 1] + (t[i] - t[i - 1]) * (FHV * (CL[i - 1] - CB[i - 1]) + QIV) / VB}
)


for (i in range) {
  CB[i] = CB[i - 1] + (t[i] - t[i - 1]) * (FHV * (CL[i - 1] - CB[i - 1]) + (if (t[i - 1] <= 2 * 60) (QIV) else (0))) / VB
  CL[i] = CL[i - 1] + (t[i] - t[i - 1]) * (FHV * (CB[i - 1] - CL[i - 1]) - VMax * CL[i - 1] / (KM + CL[i - 1])) / VL
}

data_str = """0.083 0.063 
0.25 0.16 
0.5 0.26 
0.75 0.37
1.0 0.46
1.5 0.59
2.0 0.74
2.083 0.70
2.167 0.66
2.25 0.63
2.5 0.58
2.75 0.55
3.0 0.50
3.5 0.43
4.0 0.35
4.5 0.28
5.0 0.18
5.25 0.15
5.5 0.11
5.75 0.076
6.0 0.052
6.25 0.036
6.5 0.025
6.75 0.013
7.0 0.010
7.25 0.0066
7.5 0.0040
7.75 0.0028"""



install.packages("rxode2")
library(rxode2)



#> rxode2 2.0.13.9000 using 8 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`

mod1 <- function() {
  ini({
    FHV = 1.5 #liters/min
    VB = 48 #liters
    VL = 0.61 #liters
    Vmax = 0.00275 #2.75 mmol/min/80 kg 
    Km = 0.0001 #0.1 mM 
    QIV = 720 * 0.08 / (2 * 60) * 0.789 / 46.069
  })
  model({
    CB <- centr/VB
    CL <- peri/VL
    d/dt(centr) <- FHV*(CL - CB) + QIV
    d/dt(peri)  <- FHV*(CB - CL) - Vmax*CL/(Km + CL)
  })
}



et <- eventTable();
et$add.sampling(seq(0, 7, length.out=200));
et$add.dosing(20, start.time=0);

transit <- rxSolve(mod1, et)

plot(transit, cen, ylab="Central Concentration")


mod1 <- mod1() # create the ui object (can also use `rxode2(mod1)`)
mod1

summary(mod1$simulationModel)
m1 <- rxode2(mod1)
print(m1)


mod <- RxODE({
  ## Table 3 from Savic 2007
  cl = 17.2 # (L/hr)
  vc = 45.1 # L
  ka = 0.38 # 1/hr
  mtt = 0.37 # hr
  bio=1
  n = 20.1
  k = cl/vc
  ktr = (n+1)/mtt
  d/dt(depot) = transit(n,mtt,bio)-ka*depot
  d/dt(cen) = ka*depot-k*cen
  
  
})

et <- eventTable();
et$add.sampling(seq(0, 7, length.out=200));
et$add.dosing(20, start.time=0);

transit <- rxSolve(mod, et)