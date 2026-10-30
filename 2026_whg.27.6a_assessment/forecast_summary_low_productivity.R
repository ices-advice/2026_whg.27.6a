# Process SAM forecast outputs ####

rm(list=ls()) 

library(stockassessment)
library(TMB)
library(icesAdvice)
library(ggplot2)
library(S7)

root.dir <- getwd()

#sourceDir("utils")

#Resdir <- "results"

# load SAM outputs and FLStock object
load("../SAM_output/fit.RData")
#load("N:/STOCK_ASSESSMENT/West of Scotland whiting/2024/SAM/user348-whiting_27.6a_2024_postADG/whiting_27.6a_2024_postADG/run/forecast.Rdata")
load("./results/forecast_lowprod.RData")
#load("whiting_27.6a_2024_postADG/run/model.RData")
#load("whiting_27.6a_2024_postADG/run/forecast.RData")

# Load change in advice data
adviceChange26 <- read.csv("./results/n-at-age_WGCSE_2026_lowprod.csv")
adviceChange25 <- read.csv("./results/n-at-age_WGCSE_2025.csv")
adviceChangeN <- rbind(adviceChange25, adviceChange26)
adviceChangeN$type <- as.factor(adviceChangeN$type)

adviceChange26 <- read.csv("./results/B-at-age_WGCSE_2026_lowprod.csv")
adviceChange25 <- read.csv("./results/B-at-age_WGCSE_2025.csv")
adviceChangeB <- rbind(adviceChange25, adviceChange26)
adviceChangeB$type <- as.factor(adviceChangeB$type)

adviceChange26 <- read.csv("./results/catch-at-age_WGCSE_2026_lowprod.csv")
adviceChange25 <- read.csv("./results/catch-at-age_WGCSE_2025.csv")
adviceChangeC <- rbind(adviceChange25, adviceChange26)
adviceChangeC$type <- as.factor(adviceChangeC$type)

# Check list object indices
lapply(FC, function(x){return(attributes(x)$label)})
attr(FC[[7]], "shorttab")

## ===============================================================
## Helpers: extract historical + forecast time series
## ===============================================================

get_hist_ts <- function(fit) {
  
  ssb_hist   <- ssbtable(fit)
  fbar_hist  <- fbartable(fit)
  catch_hist <- catchtable(fit)
  
  years <- Reduce(
    intersect,
    list(
      rownames(ssb_hist),
      rownames(fbar_hist),
      rownames(catch_hist)
    )
  )
  years <- as.numeric(years)
  
  list(
    ssb = data.frame(
      year = years,
      median = ssb_hist[as.character(years), "Estimate"],
      low    = ssb_hist[as.character(years), "Low"],
      high   = ssb_hist[as.character(years), "High"]
    ),
    fbar = data.frame(
      year = years,
      median = fbar_hist[as.character(years), "Estimate"],
      low    = fbar_hist[as.character(years), "Low"],
      high   = fbar_hist[as.character(years), "High"]
    ),
    catch = data.frame(
      year = years,
      median = catch_hist[as.character(years), "Estimate"],
      low    = catch_hist[as.character(years), "Low"],
      high   = catch_hist[as.character(years), "High"]
    )
  )
}

get_fc_ts <- function(fc) {
  
  tab <- as.data.frame(attr(fc, "tab"))
  tab$year <- as.numeric(rownames(tab))
  
  list(
    ssb = data.frame(
      year = tab$year,
      median = tab$`ssb:median`,
      low    = tab$`ssb:low`,
      high   = tab$`ssb:high`
    ),
    fbar = data.frame(
      year = tab$year,
      median = tab$`fbar:median`,
      low    = tab$`fbar:low`,
      high   = tab$`fbar:high`
    ),
    catch = data.frame(
      year = tab$year,
      median = tab$`catch:median`,
      low    = tab$`catch:low`,
      high   = tab$`catch:high`
    )
  )
}

### One ggplot function for every scenario

plot_ts_gg <- function(hist, fc, ylab, show_x = FALSE) {
  
  ggplot() +
    geom_ribbon(
      data = hist,
      aes(year, ymin = low, ymax = high),
      fill = "grey80"
    ) +
    geom_line(
      data = hist,
      aes(year, median),
      linetype = "dashed",
      linewidth = 0.8
    ) +
    geom_errorbar(
      data = fc,
      aes(year, ymin = low, ymax = high),
      width = 0
    ) +
    geom_line(
      data = fc,
      aes(year, median),
      linewidth = 0.9
    ) +
    geom_point(
      data = fc,
      aes(year, median),
      size = 2
    ) +
    labs(
      y = ylab,
      x = if (show_x) "Year" else NULL
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 10)
    )
}

### Wrapper to plot any FC scenario

plot_forecast_scenario <- function(fc, fit, filename) {
  
  hist <- get_hist_ts(fit)
  fc_ts <- get_fc_ts(fc)
  
  p1 <- plot_ts_gg(hist$ssb,   fc_ts$ssb,   "SSB")
  p2 <- plot_ts_gg(hist$fbar,  fc_ts$fbar,  expression(F[1-3]))
  p3 <- plot_ts_gg(hist$catch, fc_ts$catch, "Catch", show_x = TRUE)
  
  jpeg(
    filename,
    height = 210, width = 210, units = "mm", res = 350
  )
  
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 1)))
  print(p1, vp = grid::viewport(layout.pos.row = 1))
  print(p2, vp = grid::viewport(layout.pos.row = 2))
  print(p3, vp = grid::viewport(layout.pos.row = 3))
  
  dev.off()
}

## ===============================================================
## Unified ggplot-style forecast plots
## ===============================================================

# FMSY
plot_forecast_scenario(
  FC[[1]],
  fit,
  "FMSY_forecast_gg_lowprod.jpeg"
)

# Fsq
plot_forecast_scenario(
  FC[[7]],
  fit,
  "FSQ_forecast_gg_lowprod.jpeg"
)

# HCR / MSY rule
plot_forecast_scenario(
  FC[[12]],
  fit,
  "HCR_forecast_gg_lowprod.jpeg"
)


# Some summary values
assessment_summary <- as.data.frame(summary(fit))
assessment_summary  #For Standard Graph input


# TAC & advice for 2026
prev.advice <- 5364 # 5116 in 2025
prev.TAC <- 5364 # 4952 in 2025

# Create the tables ####

# Final year of assessment
final.yr <- max(fit$data$years)
# If the assessment includes the int year (and this is the 1st year in forecast 
# then set int.yr=1)
int.yr <- 2

Table2 <- NULL
Table3 <- NULL
for (i in 1:length(FC)){#i=1
  Table2 <- rbind(Table2, data.frame(Basis = attributes(FC[[i]])$label,
                                     Ftot = attributes(FC[[i]])$shorttab[1,int.yr],
                                     SSB = attributes(FC[[i]])$shorttab[3,int.yr+1],
                                     R1 = attributes(FC[[i]])$shorttab[2,int.yr],
                                     R2 = attributes(FC[[i]])$shorttab[2,int.yr+1],
                                     TotCatch = attributes(FC[[i]])$shorttab[7,int.yr]+attributes(FC[[i]])$shorttab[8,int.yr],
                                     Landings = attributes(FC[[i]])$shorttab[7,int.yr],
                                     Discards = attributes(FC[[i]])$shorttab[8,int.yr]
  ))
  
  Table3 <- rbind(Table3, data.frame(Basis = attributes(FC[[i]])$label,
                                     TotCatch = attributes(FC[[i]])$shorttab[7,int.yr+1]+attributes(FC[[i]])$shorttab[8,int.yr+1],
                                     WCatch = attributes(FC[[i]])$shorttab[7,int.yr+1],
                                     UCatch = attributes(FC[[i]])$shorttab[8,int.yr+1],
                                     Ftot = attributes(FC[[i]])$shorttab[1,int.yr+1],
                                     Fw = attributes(FC[[i]])$shorttab[5,int.yr+1],
                                     Fu = attributes(FC[[i]])$shorttab[6,int.yr+1],
                                     SSB = attributes(FC[[i]])$shorttab[3,int.yr+2],
                                     SSBpc = NA,
                                     TACpc = NA,
                                     Advicepc = NA
                                     
  ))
}

# Calculate Discards
Table2$Discards <- Table2$TotCatch - Table2$Landings
Table2 # check rows are all the same

# Catch scenario summary table
Table3
Table3$Fu <- Table3$Ftot - Table3$Fw
Table3$UCatch <-Table3$TotCatch - Table3$WCatch
Table3$SSBpc <- 100 * (Table3$SSB/Table2$SSB -1)
Table3$TACpc <- 100 * (Table3$TotCatch/prev.TAC - 1)
Table3$Advicepc <- 100 * (Table3$TotCatch/prev.advice - 1)
write.csv(Table2, file = "Table2_lowprod.csv")

# Create ICES rounded table
Tab3 <- Table3
Tab3$Ftot <- icesRound(Table3$Ftot)
Tab3$Fw <- icesRound(Table3$Fw)
Tab3$Fu <- icesRound(Table3$Fu)
Tab3$TotCatch <- round(Table3$TotCatch)
Tab3$WCatch <- round(Table3$WCatch)
Tab3$UCatch <- round(Table3$UCatch)
Tab3$SSB <- round(Table3$SSB)
Tab3$SSBpc <- icesRound(Table3$SSBpc,percent=TRUE,sign=FALSE)
Tab3$TACpc <- icesRound(Table3$TACpc,percent=TRUE,sign=FALSE)
Tab3$Advicepc <- icesRound(Table3$Advicepc,percent=TRUE,sign=FALSE)
write.csv(Tab3, file = "Table3_lowprod.csv")

# Getting median values where rounding to 4 places is required
round(median(FC[[7]][[3]]$fbarL), 4)
round(median(FC[[12]][[3]]$fbarL), 4)

# Use three year average for mean weights & L-D partition 
# (see above/below for age 1 & 2 fix, WGCSE 2022)
av.yrs <- (final.yr-3):(final.yr-1)

## WG agreed forecast fixes
# Replace mean weight-at-age 1 for 2020 with 2019 value
fit$data$catchMeanWeight["2020", "1",1] <- fit$data$catchMeanWeight["2019", "1",1]

# 2021 catchMeanWeight & stockMeanWeight ages 0 & 1 replaced with mean of previous 
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

# Get catch/stock weights used in forecast
cw <-fit$data$catchMeanWeight
cwa <- colMeans(cw[rownames(cw) %in% av.yrs,,1])
sw <-fit$data$stockMeanWeight
swa <- colMeans(sw[rownames(sw) %in% av.yrs,])


## Contribution of recruitment assumption to catch and SSB
# FMSY approach
lapply(FC, function(x){return(attributes(x)$label)})
i <- 12
attr(FC[[i]], "shorttab")

# attributes(FC[[i]])$label
FC[[i]][[3]]$year


ssb_sim <-sweep(exp(FC[[i]][[4]]$sim[,1:8]), MARGIN = 2, FUN = "*", STATS = (swa * tail(fit$data$propMat,1)))
ssb <- ssb_sim[which.min(abs(rowSums(ssb_sim) - median(rowSums(ssb_sim)))),]
sum(ssb) # check ok

#C by age Intplus1
catch_sim <- sweep(FC[[i]][[3]]$catchatage, MARGIN = 1, FUN = "*", STATS = cwa)
catch <-  catch_sim[,which.min(abs(colSums(catch_sim) - median(colSums(catch_sim))))] 
sum(catch)


#  exp(apply(FC[[i]][[4]]$sim, 2, median)) [1:8]  * colMeans(tail(fit$data$catchMeanWeight * fit$data$propMat[-dim(fit$data$propMat)[1],],1))   #stockn


year <- 2026
plot.data <- rbind(
  data.frame(recruitment = year + 1:-6, val = catch, type = paste(year+1, "Catch")),
  data.frame(recruitment = year + 2:-5, val = ssb, type = paste(year+2, "SSB"))
)

plot.data$est <- plot.data$recruitment >= year

# ggplot(plot.data, aes(x=factor(recruitment), y=val, fill=est)) + 
#   facet_wrap(~type, scales="free_y") + geom_col() + 
#   xlab("Recruitment year") +ylab("Contribution (tonnes)") +
#   scale_fill_discrete(name="Forecast assumption")
# 
# ggsave("Forecast makeup.png", 
#        width = 20, height = 15, units = "cm", dpi = 300, type = "cairo-png")

jpeg("forecast_contribution_lowprod.jpeg", 
     height = 100, width = 250, units = "mm", res = 350)
par(mfrow = c(1, 2), mar = c(0, 4.1, 1.1, 1.1), oma = c(5.1, 0, 0, 0))
barplot(height = rev(plot.data$val[plot.data$type == "2027 Catch"]),
        names.arg = rev(plot.data$recruitment[plot.data$type == "2027 Catch"]),
        col = c("#56B4E9", "#E69F00")[as.factor(rev(plot.data$est[plot.data$type == "2027 Catch"]))],
        ylab = "Contribution (tonnes)")
legend("topright",
       legend = "2027 Catch",
       cex = 1.5,
       text.font = 2,
       bty = "n")
barplot(height = rev(plot.data$val[plot.data$type == "2028 SSB"]),
        names.arg = rev(plot.data$recruitment[plot.data$type == "2028 SSB"]),
        col = c("#56B4E9", "#E69F00")[as.factor(rev(plot.data$est[plot.data$type == "2028 SSB"]))])
legend("topleft",
       legend = "2028 SSB",
       cex = 1.5,
       text.font = 2,
       bty = "n")
mtext("Recruitment Year", side = 1, outer = TRUE, line = 4, cex = 1.5)
dev.off()

# % contributions to SSB
#2027 contribution
rev(plot.data$val[plot.data$type == "2028 SSB"])[6]/sum(rev(plot.data$val[plot.data$type == "2028 SSB"]))
#2026 contribution
rev(plot.data$val[plot.data$type == "2028 SSB"])[7]/sum(rev(plot.data$val[plot.data$type == "2028 SSB"]))
#2025 contribution
rev(plot.data$val[plot.data$type == "2028 SSB"])[5]/sum(rev(plot.data$val[plot.data$type == "2028 SSB"]))

# % contributions to Catch
rev(plot.data$val[plot.data$type == "2027 Catch"])[7]/sum(rev(plot.data$val[plot.data$type == "2027 Catch"]))
rev(plot.data$val[plot.data$type == "2027 Catch"])[8]/sum(rev(plot.data$val[plot.data$type == "2027 Catch"]))

## Change in advice plot
pointType = c(19, 0, 2)

# Numbers
jpeg("../figures/compare_stock_n_at_age_lowprod.jpeg", 
     height = 150, width = 150, units = "mm", res = 350)
#windows()
par(mfrow = c(2, 2), mar = c(0, 0, 0, 0), oma = c(5.1, 4.1, 1.1, 1.1))
for(i in 2024:2027){
  #i = 2023
  xx <- adviceChangeN[adviceChangeN$year == i,]
  plot(xx$age[xx$WG == "WGCSE 2025"], xx$N[xx$WG == "WGCSE 2025"],
       pch = pointType[xx$type[xx$WG == "WGCSE 2025"]],
       xlab = "",
       ylab = "",
       # xlim = c(-0.5, 7.5),
       ylim = c(0, 3.5e5),
       col = "#E69F00",
       axes = FALSE,
       type = "b")
  lines(xx$age[xx$WG == "WGCSE 2026"], xx$N[xx$WG == "WGCSE 2026"],
        pch = pointType[xx$type[xx$WG == "WGCSE 2026"]],
        col = "#56B4E9",
        type = "b")
  text(0, 3.2e5, labels = i, cex = 1.5, font = 2, pos = 4)
  box()
  if(i %in% c(2024, 2026)){axis(2)}
  if(i %in% c(2026, 2027)){axis(1)}
  if(i == 2025){legend("topright",
                       legend = c("WGCSE 25", "WGCSE 26", 
                                  "Data", "Intermediate", "Forecast"),
                       bty = "n",
                       y.intersp = 1.2,
                       col = c("#E69F00", "#56B4E9", rep("#000000", 3)),
                       lty = c(1, 1, NA, NA, NA),
                       pch = c(NA, NA, 19, 2, 0))}
}
mtext("Numbers-at-age", outer = TRUE, side = 2, line = 2.5)
mtext("Age Class", outer = TRUE, side = 1, line = 3)
dev.off()

# Biomass
jpeg("../figures/compare_stock_b_at_age_lowprod.jpeg", 
     height = 150, width = 150, units = "mm", res = 350)
par(mfrow = c(2, 2), mar = c(0, 0, 0, 0), oma = c(5.1, 4.1, 1.1, 1.1))
for(i in 2024:2027){
  #i = 2023
  xx <- adviceChangeB[adviceChangeB$year == i,]
  plot(xx$age[xx$WG == "WGCSE 2025"], xx$N[xx$WG == "WGCSE 2025"],
       pch = pointType[xx$type[xx$WG == "WGCSE 2025"]],
       xlab = "",
       ylab = "",
       # xlim = c(-0.5, 7.5),
       ylim = c(0, 1.2e4),
       col = "#E69F00",
       axes = FALSE,
       type = "b")
  lines(xx$age[xx$WG == "WGCSE 2026"], xx$N[xx$WG == "WGCSE 2026"],
        pch = pointType[xx$type[xx$WG == "WGCSE 2026"]],
        col = "#56B4E9",
        type = "b")
  text(0, 1.1e4, labels = i, cex = 1.5, font = 2, pos = 4)
  box()
  if(i %in% c(2024, 2026)){axis(2)}
  if(i %in% c(2026, 2027)){axis(1)}
  if(i == 2025){legend("topright",
                       legend = c("WGCSE 25", "WGCSE 26", 
                                  "Data", "Intermediate", "Forecast"),
                       bty = "n",
                       y.intersp = 1.2,
                       col = c("#E69F00", "#56B4E9", rep("#000000", 3)),
                       lty = c(1, 1, NA, NA, NA),
                       pch = c(NA, NA, 19, 2, 0))}
}
mtext("Biomass-at-age", outer = TRUE, side = 2, line = 2.5)
mtext("Age Class", outer = TRUE, side = 1, line = 3)
dev.off()

# Catch
jpeg("../figures/compare_catch_at_age_lowprod.jpeg", 
     height = 150, width = 150, units = "mm", res = 350)
par(mfrow = c(2, 2), mar = c(0, 0, 0, 0), oma = c(5.1, 4.1, 1.1, 1.1))
for(i in 2024:2027){
  #i = 2023
  xx <- adviceChangeC[adviceChangeC$year == i,]
  plot(xx$age[xx$WG == "WGCSE 2025"], xx$N[xx$WG == "WGCSE 2025"],
       pch = pointType[xx$type[xx$WG == "WGCSE 2025"]],
       xlab = "",
       ylab = "",
       # xlim = c(-0.5, 7.5),
       ylim = c(0, 2.1e4),
       col = "#E69F00",
       axes = FALSE,
       type = "b")
  lines(xx$age[xx$WG == "WGCSE 2026"], xx$N[xx$WG == "WGCSE 2026"],
        pch = pointType[xx$type[xx$WG == "WGCSE 2026"]],
        col = "#56B4E9",
        type = "b")
  text(0, 2.1e4, labels = i, cex = 1.5, font = 2, pos = 4)
  box()
  if(i %in% c(2024, 2026)){axis(2)}
  if(i %in% c(2026, 2027)){axis(1)}
  if(i == 2025){legend("topright",
                       legend = c("WGCSE 25", "WGCSE 26", 
                                  "Data", "Intermediate", "Forecast"),
                       bty = "n",
                       y.intersp = 1.2,
                       col = c("#E69F00", "#56B4E9", rep("#000000", 3)),
                       lty = c(1, 1, NA, NA, NA),
                       pch = c(NA, NA, 19, 2, 0))}
}
mtext("Catch-at-age", outer = TRUE, side = 2, line = 2.5)
mtext("Age Class", outer = TRUE, side = 1, line = 3)
dev.off()
