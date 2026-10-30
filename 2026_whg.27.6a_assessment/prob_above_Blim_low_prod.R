### Calculate prob SSB > Blim in 2026 and add to forecast outputs (Table 3)

rm(list=ls()) 

library(stockassessment)
library(icesAdvice)

# Last year of data (= year prior assessment year)
data.year <- 2025 
assessment.year <- 2026

#  Reference points - WKNSEA 2021 Benchmark
F_MSY <- 0.21
F_MSY_lower <- 0.173
F_MSY_upper <- 0.21
Fpa <- 0.21
Flim <- 0.31
Blim <- 17286
Bpa <- 25597
MSY_Btrigger <- 25597

# Load stock weights data 
sw <- read.ices("../data/sw.dat")

# Load model results
load("../SAM_output/fit.RData")

# Load forecast results
load("./results/forecast_lowprod.RData")

# Load Table 3
Table3 <- read.csv("Table3_lowprod.csv")

# Final year in the assessment
final.yr <- max(fit$data$years)

# Resample R from previous 10 years Rec
R.range <-(final.yr-6):(final.yr-1)

# Use three year average for mean weights & L-D partition 
# (see above/below for age 1 & 2 fix, WGCSE 2022)
av.yrs <- (final.yr-3):(final.yr-1)

# Use 5 year average selection patter
sel.yrs <- (final.yr-5):(final.yr-1)

# Get stock weights used in forecast
swa <- colMeans(sw[as.character(av.yrs),])

# Get prob SSB > Blim in 2028
perc <-rep(NA,length(FC))
for (i in 1:length(FC)){
  perc[i] <-100*sum(rowSums(sweep(exp(FC[[i]][[4]]$sim[,1:8]), MARGIN = 2, FUN = "*", STATS = (swa * tail(fit$data$propMat,1))))>Blim)/dim(FC[[i]][[4]]$sim)[1]
}

# Get prob SSB < Blim in 2028
perc <- 100-perc

Table3$Prob.less.Blim <-icesRound(perc)
write.csv(Table3, file = "Table3_lowprod_final.csv")

