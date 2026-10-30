
### Forecast run with recruitment resampled only from last 6 years to account for lower productivity

rm(list=ls())
library(stockassessment)
load("../SAM_output/fit.RData")

# Model output summary
fitsum <- as.data.frame(cbind(summary(fit), rbind(catchtable(fit), NA)))
fitsum$year <- as.numeric(rownames(fitsum))
names(fitsum) <- c("Rec","RecLow", "RecHigh",
                   "SSB", "SSBLow", "SSBHigh",
                   "Fbar", "FLow", "FHigh",
                   "Catch", "CatchLow", "CatchHigh", "year")


### Input Parameters

#  Reference points - WKNSEA 2021 Benchmark
F_MSY <- 0.21
F_MSY_lower <- 0.173
F_MSY_upper <- 0.21
Fpa <- 0.21
Flim <- 0.31
Blim <- 17286
Bpa <- 25597
MSY_Btrigger <- 25597


# TAC & advice for 2026
prev.advice <- 5364 # 5116 in 2025
prev.TAC <- 5364 # 4952 in 2025

# Current catch
#cur.cat <- 710 # 751 in 2022

# Number of simulations
numberOfSims <- 1e3


###  FORECAST ASSUMPTIONS
# Modified WGCSE 2022
# - No intermediate year in the assessment
# - exclusion of age 1 & age 2 discard weight info (catch mean wts, stock mean wts)
#   in final assessment year from 3 year average - fix is to replace
#   final year with average of previous 2

# Final year in the assessment
final.yr <- max(fit$data$years)

#  Start year of forecast - needs to be year before intermediate year
start.yr <- final.yr-1

# Resample R from previous 6 years Rec to match current low productivity
R.range <-(final.yr-6):(final.yr-1)

# Use three year average for mean weights & L-D partition 
# (see above/below for age 1 & 2 fix, WGCSE 2022)
av.yrs <- (final.yr-3):(final.yr-1)

# Use 5 year average selection patter
sel.yrs <- (final.yr-5):(final.yr-1)

# F status quo is last observed value
f_sq <- fbartable(fit)[as.character(final.yr-1), "Estimate"]
# f_mean <- mean(fbartable(fit)[as.character(sel.yrs),"Estimate"])

#----- Final assessment year F - needed to run the forecast through 2021 to get numbers/SSB at start of intermediate year
f_final <- f_sq
# 
# fit.save <- fit

## WG agreed forecast fixes
# Replace mean weight-at-age 1 for 2020 with 2019 value
fit$data$catchMeanWeight["2020", "1",1] <- fit$data$catchMeanWeight["2019", "1",1]


## 2021 catchMeanWeight & stockMeanWeight ages 0 & 1 replaced with 
# mean of previous; WGCSE 2022
ags <- as.character(0:1)
rep.yrs <- c(2019, 2020)

# No need to adjust discMeanWeight as not used in the forecast 
# (Total disc = Catch - Land)
fit$data$catchMeanWeight["2021", ags,1] <- colMeans(fit$data$catchMeanWeight[as.character(rep.yrs),ags,1])
fit$data$stockMeanWeight["2021", ags] <- colMeans(fit$data$stockMeanWeight[as.character(rep.yrs),ags])


## 2022 catchMeanWeight & landMeanWeight replaced with 
# mean of previous 3; WGCSE 2023
ags <- as.character(0:7)
rep.yrs <- c(2019:2021)

# No need to adjust discMeanWeight as not used in the forecast 
# (Total disc = Catch - Land)
fit$data$catchMeanWeight["2022", ags,1] <- colMeans(fit$data$catchMeanWeight[as.character(rep.yrs),ags,1])
fit$data$landMeanWeight["2022", ags,1] <- colMeans(fit$data$landMeanWeight[as.character(rep.yrs),ags,1])


### FORECAST

scen <- list(
  "MSY approach: F_MSY" = list(fval = c(f_sq, f_sq, F_MSY, F_MSY)),
  "Precautionary approach: Fpa" = list(fval = c(f_sq, f_sq, Fpa, Fpa)),
  "FMSY upper" = list(fval = c(f_sq, f_sq, F_MSY_upper, F_MSY_upper)),
  "FMSY lower" = list(fval = c(f_sq, f_sq, F_MSY_lower, F_MSY_lower)),
  "F = 0" = list(fval = c(f_sq, f_sq, NA, NA), fscale=c(NA,NA,0.0,0.0)),
  "F = Flim" = list(fval = c(f_sq, f_sq, Flim, Flim)),
  "Fsq"  = list(fval = c(f_sq, f_sq, f_sq, f_sq)),
  "0.05*Fsq" = list(fval = c(f_sq, f_sq, 0.05*f_sq, 0.05*f_sq)),
  "0.25*Fsq" = list(fval = c(f_sq, f_sq, 0.25*f_sq, 0.25*f_sq)),
  "0.5*Fsq" = list(fval = c(f_sq, f_sq, 0.5*f_sq, 0.5*f_sq)),
  "0.75*Fsq" = list(fval = c(f_sq, f_sq, 0.75*f_sq, 0.75*f_sq))
)


# Forecast storage list
FC <- list()
#FC <- vector("list", length(scen))
#names(FC) <- names(scen)

# Run F scenarios
for(i in seq(scen)){ # i<-1
  RNGkind(sample.kind = "Rounding")
  set.seed(12345)
  ARGS <- scen[[i]] # fval = ARGS$fval
  ARGS <- c(ARGS,
            list(fit = fit, 
                 ave.years = av.yrs, 
                 rec.years = R.range, 
                 year.base = start.yr,
                 label = names(scen)[i], 
                 overwriteSelYears = sel.yrs, 
                 splitLD = TRUE,
                 nosim = numberOfSims,
                 savesim = TRUE))  
  FC[[i]] <- do.call(forecast, ARGS) ###
  
  print(paste0("forecast : ", "'", names(scen)[i], "'", " is complete"))
}

# load forecast function which saves the naytable as an attribute
source("./forecast_sam_with_natage_v3.r") 
#source("N:/STOCK_ASSESSMENT/West of Scotland whiting/2026/SAM/forecast/forecast_sam_with_natage_v3.r")

# load custom function to debug
source("./processNoiseF_v2.r")
#source("N:/STOCK_ASSESSMENT/West of Scotland whiting/2026/SAM/forecast/processNoiseF_v2.r")

# ICES advice rule for MSY
MSYappr <- function(fit, 
                    Fmsy = NULL, # Fmsy = 0.21
                    Btrig = NULL, # Btrig = 25597
                    fscale_init = 1, 
                    fval_init =NA, 
                    catchval_init = NA, 
                    label = NA){
  
  fscale <- rep(fscale_init,2)
  fval <- rep(fval_init,2)
  catchval <- rep(catchval_init,2)
  
  for (i in 3:4){
    
    fscale[i] <- NA
    fval[i] <- Fmsy
    catchval[i] <- NA
    
    set.seed(12345)
    f1 <- forecast_sam_with_natage(fit, fscale = fscale, fval = fval, catchval = catchval, year.base = start.yr, ave.years = av.yrs, rec.years = R.range, overwriteSelYears = sel.yrs)
    SSB <- attr(f1, "shorttab")[3,i]
    fval[i] <- min(Fmsy*SSB/Btrig, Fmsy)
    
  }
  
  set.seed(12345)
  f_final <- forecast_sam_with_natage(fit, 
                      fscale = fscale, 
                      fval = fval, 
                      catchval = catchval, 
                      year.base = start.yr,
                      ave.years = av.yrs, 
                      rec.years = R.range, 
                      label = label, 
                      overwriteSelYears = sel.yrs, 
                      splitLD = TRUE, 
                      nosim = numberOfSims,
                      savesim = TRUE)
  return(f_final)
}

FC[[length(FC)+1]] <- MSYappr(fit, Fmsy = F_MSY, Btrig = MSY_Btrigger, fscale_init = NA, fval_init = f_sq, catchval_init = NA, label = paste0(final.yr, "F=Fsq then Fmsy rule"))

FC[[length(FC)+1]] <- MSYappr(fit, Fmsy = F_MSY_lower, Btrig = MSY_Btrigger, fscale_init = NA, fval_init = f_sq, catchval_init = NA, label = paste0(final.yr, "F=Fsq then Fmsy rule lower"))

FC[[length(FC)+1]] <- MSYappr(fit, Fmsy = F_MSY_upper, Btrig = MSY_Btrigger, fscale_init = NA, fval_init = f_sq, catchval_init = NA, label = paste0(final.yr, "F=Fsq then Fmsy rule upper"))


# SSB target
# X% increase
# incr <- 10
# inc <- seq(0, 50, incr)/100
# inc <- 0 # doesn't run for increases
# SSBnext <- attributes(FC[[1]])$shorttab[3,2]

# for (i in 1:length(inc)){
#   SSBtarget <- SSBnext*(1+inc[i])
#   RNGkind(sample.kind = "Rounding")
#   set.seed(12345)
#   FC[[length(FC)+1]] <- forecast(fit, 
#                                  fscale = c(NA,NA,NA), 
#                                  fval = c(f_sq,NA,f_sq),
#                                  nextssb = c(NA,SSBtarget,NA),
#                                  year.base = final.yr,
#                                  label = paste0("2024F=Fsq then ",(i-1)*incr,"% SSB increase"),
#                                  rec.years = R.range,
#                                  ave.years = av.yrs,
#                                  overwriteSelYears = sel.yrs,
#                                  splitLD = TRUE,
#                                  nosim = numberOfSims,
#                                  #nosim = 1e5,
#                                  savesim = TRUE)
# }


# SSB target
# SSB 2028 = SSB 2027
  SSBtarget <- attributes(FC[[1]])$shorttab[3,3]
  
  ### CHECK IF TARGET ATTAINABLE ###
  
  # Maximum attainable SSB assumed to be the F=0 scenario
  maxSSB <- attr(FC[[5]], "shorttab")[3,4]
  
  if (maxSSB >= SSBtarget) {
  
  RNGkind(sample.kind = "Rounding")
  set.seed(12345)
  FC[[length(FC)+1]] <- forecast(fit,
                                 fscale = c(NA,NA,NA,NA),
                                 fval = c(f_sq,f_sq,NA,f_sq),
                                 nextssb = c(NA,NA,SSBtarget,NA),
                                 year.base = start.yr,
                                 label = paste0(start.yr, "F=Fsq then SSB(2028)=SSB(2027)"),
                                 rec.years = R.range,
                                 ave.years = av.yrs,
                                 overwriteSelYears = sel.yrs,
                                 splitLD = TRUE,
                                 nosim = numberOfSims,
                                 #nosim = 1e5,
                                 savesim = TRUE)
  } else {
    message("SSB target not attainable: target = ",
      round(SSBtarget, 0),
      ", maximum SSB at F=0 = ",
      round(maxSSB, 0)) 
  }

# IntermediateF=Fsq then SSB=Blim
RNGkind(sample.kind = "Rounding")
set.seed(12345)
FC[[length(FC)+1]] <- forecast_sam_with_natage(fit,
                               fscale = c(NA,NA,NA,NA),
                               fval = c(f_sq,f_sq,NA,f_sq),
                               nextssb = c(NA,NA,Blim,NA),
                               year.base = start.yr,
                               label = paste0(start.yr, "F=Fsq then SSB=Blim"), 
                               rec.years = R.range,
                               ave.years = av.yrs,
                               overwriteSelYears = sel.yrs,
                               nosim = numberOfSims,
                               splitLD = TRUE,
                               savesim = TRUE)
FC[[length(FC)]]
attr(FC[[length(FC)]], "shorttab")


# IntermediateF=Fsq then SSB=Bpa
RNGkind(sample.kind = "Rounding")
set.seed(12345)
  FC[[length(FC)+1]] <- forecast_sam_with_natage(
    fit,
    fscale = c(NA,NA,NA,NA),
    fval = c(f_sq,f_sq, NA,f_sq),
    nextssb = c(NA,NA,Bpa,NA),
    year.base = start.yr,
    label = paste0(final.yr, "F=Fsq then SSB=Bpa"),
    rec.years = R.range,
    ave.years = av.yrs,
    overwriteSelYears = sel.yrs,
    nosim = numberOfSims,
    splitLD = TRUE,
    savesim = TRUE
  )
FC[[length(FC)]]
attr(FC[[length(FC)]], "shorttab")


# IntermediateF=Fsq then SSB=MSYBtrigger
# RNGkind(sample.kind = "Rounding")
# set.seed(12345)
# FC[[length(FC)+1]] <- forecast(fit,
#                                fscale = c(NA,NA,NA),
#                                fval = c(f_sq,NA,f_sq),
#                                nextssb = c(NA,MSY_Btrigger,NA),
#                                year.base = final.yr,
#                                label = paste0(final.yr, "F=Fsq then SSB=MSYBtrigger"), 
#                                rec.years = R.range,
#                                ave.years = av.yrs,
#                                overwriteSelYears = sel.yrs,
#                                nosim = numberOfSims,
#                                splitLD = TRUE,
#                                savesim = TRUE) # Currently cannot be achieved
# FC[[length(FC)]]
# attr(FC[[length(FC)]], "shorttab")


# TAC reduction
# X% increase

inc <- seq(0.25,1,0.25)

for (i in 1:length(inc)){
  Catchtarget <- prev.TAC*inc[i]
  RNGkind(sample.kind = "Rounding")
  set.seed(12345)
  FC[[length(FC)+1]] <- forecast_sam_with_natage(fit, 
                                 fval = c(f_sq,f_sq,NA,NA), 
                                 catchval = c(NA,NA,prev.TAC*inc[i],prev.TAC*inc[i]), 
                                 year.base = start.yr,
                                 label = paste0("2025F=Fsq then prev.TAC*",inc[i]), 
                                 rec.years = R.range,
                                 ave.years = av.yrs,
                                 overwriteSelYears = sel.yrs,
                                 splitLD = TRUE,
                                 nosim = numberOfSims,
                                 savesim = TRUE)
}  



# nms <- names(scen)
# LM <-length(scen)
# for (i in (LM+1):length(FC)){
#   nms <-c(nms,attr(FC[[i]],"label"))
# }
# names(FC) <-nms


# Check list object indices
lapply(FC, function(x){return(attributes(x)$label)})
attr(FC[[13]], "shorttab")

# Change in advice values
# Numbers-at-age
adviceChangeN <- as.data.frame(attr(FC[[12]], "naytable"))
# This will have several rows with "median.50%", "low.2.5%", "high.97.5%". 
# Each set of these rows is for a different forecast year
adviceChangeN <- adviceChangeN[grep("med", rownames(adviceChangeN)),]

# Get numbers-at-age from data years
datYears <- as.data.frame(ntable(fit))[rownames(ntable(fit)) %in% c(final.yr-3, final.yr-2),] ### CHANGE BECAUSE WE CHANGED START YEAR
names(datYears) <- names(adviceChangeN)

# Combine and format
adviceChangeN <- rbind(datYears, adviceChangeN)
adviceChangeN$year <- c((final.yr-3):(final.yr+2))
adviceChangeN$WG <- paste("WGCSE", final.yr)
adviceChangeN$type <- c("data", "data", "data",
                        "intermediate year", 
                        "forecast","forecast")

# Mean stock weights-at-age years
meanWts <- colMeans(fit$data$stockMeanWeight[rownames(fit$data$stockMeanWeight) %in% as.character(av.yrs),])
adviceChangeB <- t(t(as.matrix(adviceChangeN[,c(1:8)])) * meanWts)
dataB <- adviceChangeN[c(1,2),c(1:8)]
dataMeanWts <- fit$data$stockMeanWeight[rownames(dataB),]
dataB <- dataB * dataMeanWts
adviceChangeB <- rbind(dataB, adviceChangeB[3:6,])
adviceChangeB <- cbind(adviceChangeB, adviceChangeN[,9:11])


# Go wide-to-long
adviceChangeN <- reshape(adviceChangeN, 
                         direction = "long",
                         varying = list(names(adviceChangeN)[1:8]),
                         v.names = "N",
                         idvar = c("year", "WG", "type"),
                         timevar = "age",
                         times = 0:7)
rownames(adviceChangeN) <- NULL
adviceChangeN <- adviceChangeN[with(adviceChangeN, order(year, age)), ]

adviceChangeB <- reshape(adviceChangeB, 
                         direction = "long",
                         varying = list(names(adviceChangeB)[1:8]),
                         v.names = "N",
                         idvar = c("year", "WG", "type"),
                         timevar = "age",
                         times = 0:7)
rownames(adviceChangeB) <- NULL
adviceChangeB <- adviceChangeB[with(adviceChangeB, order(year, age)), ]

# Catch
adviceChangeC <- as.data.frame(attr(FC[[12]], "caytable"))
adviceChangeC <- adviceChangeC[grep("med", rownames(adviceChangeC)),]

# Get catch-at-age from data years
datYears <- as.data.frame(caytable(fit))[rownames(ntable(fit)) %in% c(final.yr-3, final.yr-2),]
#names(datYears) <- names(adviceChangeN)

# Load object produced by stockassessment.org (using caytable function which doesn't work in stockassessment v9.0)
# This is because everything above was produced with stockassessment v9.0 in order to get the residuals (no clue why it doesn't work with new version)
# Here we load the one fit produced by stockassessment.org with stockassessment v12.0 (most recent) as it contains stuff needed for code below
# load("N:/STOCK_ASSESSMENT/West of Scotland whiting/2024/SAM/SAM_output/catage.Rdata")
# datYears <- as.data.frame(catage)[rownames(ntable(fit)) %in% c(final.yr-2, final.yr-1),]
# names(datYears) <- names(adviceChangeC)


# Combine and format
adviceChangeC <- rbind(datYears, adviceChangeC)
adviceChangeC$year <- c((final.yr-3):(final.yr+2))
adviceChangeC$WG <- paste("WGCSE", final.yr)
adviceChangeC$type <- c("data","data", "data",
                        "intermediate year", 
                        "forecast","forecast")
adviceChangeC <- reshape(adviceChangeC, 
                         direction = "long",
                         varying = list(names(adviceChangeC)[1:8]),
                         v.names = "N",
                         idvar = c("year", "WG", "type"),
                         timevar = "age",
                         times = 0:7)
rownames(adviceChangeC) <- NULL
adviceChangeC <- adviceChangeC[with(adviceChangeC, order(year, age)), ]


# Save change in advice output
write.csv(adviceChangeN, paste("./results/n-at-age_WGCSE_",final.yr,"_lowprod.csv",sep=''), 
          row.names = FALSE, quote = FALSE)
write.csv(adviceChangeB, paste("./results/B-at-age_WGCSE_",final.yr,"_lowprod.csv",sep=''), 
          row.names = FALSE, quote = FALSE)
write.csv(adviceChangeC, paste("./results/catch-at-age_WGCSE_",final.yr,"_lowprod.csv",sep=''), 
          row.names = FALSE, quote = FALSE)

# Save list of forecasts
save(FC, file = "./results/forecast_lowprod.RData")
