library(stockassessment)
library(gtools)

# Last year of data (= year prior assessment year)
data.year <- 2025 
assessment.year <- 2026

# Reference points
FMSY <- 0.21
MSYBtrigger <- 25597

## Create directory for master dataset
if(dir.exists("configuration_files") == FALSE){dir.create("configuration_files")}
if(dir.exists("figures") == FALSE){dir.create("figures")}
if(dir.exists("SAM_output") == FALSE){dir.create("SAM_output")}




# Import data files ####
cn <- read.ices("./data/cn.dat")
cw <- read.ices("./data/cw.dat")
dw <- read.ices("./data/dw.dat")
lf <- read.ices("./data/lf.dat")
lw <- read.ices("./data/lw.dat")
mo <- read.ices("./data/mo.dat")
nm <- read.ices("./data/nm.dat")
pf <- read.ices("./data/pf.dat")
pm <- read.ices("./data/pm.dat")
sw <- read.ices("./data/sw.dat")
surveys <- read.ices("./data/survey.dat")

# Remove catch numbers-at-age value, considered unreliable
# As per WGCSE 2022
cn["2021", c("0", "1")] <- NA

## Save Catch weights
totalCatchWeights <- data.frame(
  catch = rowSums(cn * lf * lw, na.rm = TRUE) + rowSums(cn * (1-lf) * dw, na.rm = TRUE),
  landings = rowSums(cn * lf * lw, na.rm = TRUE),
  discards = rowSums(cn * (1-lf) * dw, na.rm = TRUE)
)
write.csv(totalCatchWeights, "./SAM_output/catchWeights.csv")


## Survey variances as weights
sv <- stockassessment:::read.surveys("./data/survey_variance.txt")
W1 <- (surveys[[2]]^2)/sv[[1]]  # rationale sqrt(sv)/obs ~ cv ~ sd on log scale, hence (obs^2)/sv approx weight  
W2 <- (surveys[[3]]^2)/sv[[2]]
attr(surveys[[2]], "weight") <- W1
attr(surveys[[3]], "weight") <- W2


## Create collected data object for SAM
dat <- setup.sam.data(surveys = surveys,
                      residual.fleet = cn, 
                      prop.mature = mo, 
                      stock.mean.weight = sw, 
                      catch.mean.weight = cw, 
                      dis.mean.weight = dw, 
                      land.mean.weight = lw,
                      prop.f = pf, 
                      prop.m = pm, 
                      natural.mortality = nm, 
                      land.frac = lf)




# Model Configuration ####

# Generate a default model configuration
# This configuration can be changed by modifying/overwriting the elements in the list. 
conf <- defcon(dat)

# Change lowest and highest age included in Fbar
conf$fbarRange <- c(1, 3)

# Coupling of the fishing mortality states (normally only first row is used).
conf$keyLogFsta <- matrix(c(0, 1, 2, 3, 4, 5, 6, 7,
                            -1,-1,-1,-1,-1,-1,-1,-1,
                            -1,-1,-1,-1,-1,-1,-1,-1,
                            -1,-1,-1,-1,-1,-1,-1,-1), nrow = 4, byrow = TRUE)


# Observation variance coupling
conf$keyVarObs <- matrix(c(0,  1,  1,  1,  1,  1,  1,  2,
                           -1,  3,  3,  3,  3,  3,  3, -1,
                           -1,  4,  4,  4,  4,  4,  4, -1,
                           5,  6,  6,  6,  6,  7,  7,  7), nrow = 4, byrow = TRUE)

# Fleet covariance structure
conf$obsCorStruct <- factor(c("ID", "ID", "ID", "AR"), levels = c("ID", "AR", "US"))
cns <- colnames(conf$keyCorObs)
conf$keyCorObs <- matrix(c(NA, NA, NA, NA, NA, NA, NA,
                           -1, NA, NA, NA, NA, NA, -1,
                           -1, NA, NA, NA, NA, NA, -1,
                           0,  1,  1,  1,  1,  1,  2), nrow = 4, byrow = TRUE)
colnames(conf$keyCorObs) <- cns


## Save model configurations to files
saveConf(conf, "./configuration_files/model.cfg", overwrite = TRUE)




# Generate default initial values for all our model parameters ####
# These default initial can be modified but it is rarely necessary. 
par <- defpar(dat, conf)




# Load model output, or fit model ####
if(file.exists("./SAM_output/output.RData")){
  
  load("./SAM_output/output.RData")
  
} else {
  
  # Optimize the models
  set.seed(1337)
  fit <- sam.fit(dat, conf, par)
  
  # Get residuals
  res <- residuals(fit)
  resp <- procres(fit)
  
  # Retros
  Retmat <- cbind(data.year-0:4, 3000, assessment.year-0:4, data.year-0:4)
  mretro <- retro(fit, year = Retmat)
  lo <- leaveout(fit)
  
  # Save Model Output & Results
  save(fit, res, resp, mretro, lo,
       file = "./SAM_output/output.RData")
  
}




# Results & Model Diagnostics ####


## Model summaries
modeltable(fit)


## Parameter Tables
partable(fit)
#                                        par      sd(par)     exp(par)          Low         High
# Note in fit1 - transfIRARdist_0  3.4165138 1.382003e+03 30.463030467 0.0000000000          Inf
# Indicates weak to non-existent level of autocorrelation between age groups 0 and 1
# i.e. Model is still acceptable to the WG
write.csv(partable(fit), "./SAM_output/parameters.csv")


## Save model summaries
fitsum <- as.data.frame(cbind(summary(fit), rbind(catchtable(fit), NA)))
fitsum$year <- as.numeric(rownames(fitsum))
names(fitsum) <- c("Rec","RecLow", "RecHigh",
                   "SSB", "SSBLow", "SSBHigh",
                   "Fbar", "FLow", "FHigh",
                   "Catch", "CatchLow", "CatchHigh", "year")
write.table(fitsum, "./SAM_output/summary.txt", row.names = FALSE)

# ICES model summary
write.csv(cbind(fitsum[,c(2,1,3,5,4,6)],
                rbind(totalCatchWeights[,c(2,3)], c(NA, NA)),
                fitsum[,c(8,7,9)]), "./SAM_output/summary.csv")

# Numbers-at-age
write.csv(ntable(fit), "./SAM_output/ntable.csv")

# F-at-age
write.csv(faytable(fit), "./SAM_output/ftable.csv")


## Data fits
jpeg("./figures/fitplot_conf.jpeg", height = 210, width = 210, units = "mm", res = 350)
fitplot(fit)
dev.off()


## One-observation-ahead residuals
jpeg("./figures/res.jpeg", height = 210, width = 210, units = "mm", res = 350)
plot(res)
dev.off()


## Process Residuals
jpeg("./figures/procres.jpeg", height = 210, width = 210, units = "mm", res = 350)
plot(resp)
dev.off()


## Residual between-age correlation matrices by fleet
jpeg("./figures/empirobscorr.jpeg", height = 210, width = 210, units = "mm", res = 350)
empirobscorrplot(res)
dev.off()


## Retrospective analyses
jpeg("./figures/retro.jpeg", height = 210, width = 210, units = "mm", res = 350)

par(mfrow=c(3,1))
ssbplot(mretro)
recplot(mretro, drop=1)
fbarplot(mretro, drop=1)

dev.off()


## Mohn's retro values
mohn(mretro)["SSB"] # 2026 value is 0.07600905
mohn(mretro, lag = 1)[c("Fbar(1-3)", "R(age 0)")] # 2026 values are Fbar(1-3) 0.2816111 & R(age 0) 0.2859491 


## Leave-one-out analysis
jpeg("./figures/loo.jpeg", height = 210, width = 210, units = "mm", res = 350)
plot(lo)
dev.off()




# Save best model for forecasting ####
save(fit, file = "./SAM_output/fit.RData")




##### Some additional data visualizations

# Transparent colour function
t_col <- function(color, percent = 50, name = NULL){
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

# CVD-friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                        "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                        

## Summary plot of model estimates - Updated model
jpeg("./figures/summary.jpeg", height = 160, width = 210, units = "mm", res = 350)
par(mfrow = c(2, 2), mar = c(0, 4.1, 0, 1.1), oma = c(5.1, 0, 1.1, 0))

# Catch
plot(1981:data.year, catchtable(fit)[,1]/1e3,
     xlab = "Year",
     ylab = "Catch ('000 t)",
     ylim = c(-0.4636151, max(catchtable(fit))/1e3),
     axes = FALSE,
     lwd = 2,
     type = "n")
polygon(c(1981:data.year, rev(1981:data.year)), 
        c(catchtable(fit)[,2]/1e3, rev(catchtable(fit)[,3]/1e3)),
        border = "#ffebbf",
        col = "#ffebbf")
lines(fitsum$year, fitsum$Catch/1e3,
      lwd = 2,
      col = "#E69F00")
points(1981:data.year, rowSums(cw*cn, na.rm = TRUE)/1e3,
       cex = 0.5,
       col = "#000000",
       pch = 8)
axis(2)
box()

# Recruitment
plot(fitsum$year, fitsum$Rec,
     xlab = "Year",
     ylab = "Recruitment (Age 0)",
     ylim = c(8812.92, 2200000),
     axes = FALSE,
     col = "#000000",
     lwd = 2,
     type = "n")
polygon(c(1981:assessment.year, rev(1981:assessment.year)), 
        c(fitsum$RecLow, rev(fitsum$RecHigh)),
        border = "#ffebbf",
        col = "#ffebbf")
lines(fitsum$year, fitsum$Rec,
      lwd = 2,
      col = "#E69F00")
axis(2)
box()

# SSB
plot(fitsum$year, fitsum$SSB/1e3,
     xlab = "Year",
     ylab = "SSB ('000 t)",
     ylim = c(3, 150),
     lwd = 2,
     col = "#000000",
     type = "n")
polygon(c(1981:assessment.year, rev(1981:assessment.year)), 
        c(fitsum$SSBLow, rev(fitsum$SSBHigh))/1e3,
        border = "#ffebbf",
        col = "#ffebbf")
lines(fitsum$year, fitsum$SSB/1e3,
      lwd = 2,
      col = "#E69F00")
segments(x0 = 1981, x1 = assessment.year, 
         y0 = MSYBtrigger/1e3, y1 = MSYBtrigger/1e3, lty = 3)

# Fbar
plot(fitsum$year, fitsum$Fbar,
     xlab = "Year",
     ylab = expression("F"["1-3"]),
     lwd = 2,
     ylim = c(0, 0.85),
     col = "#000000",
     type = "n")
polygon(c(1981:assessment.year, rev(1981:assessment.year)), 
        c(fitsum$FLow, rev(fitsum$FHigh)),
        border = "#ffebbf",
        col = "#ffebbf")
lines(fitsum$year, fitsum$Fbar,
      lwd = 2,
      col = "#E69F00")
segments(x0 = 1981, x1 = assessment.year, y0 = FMSY, y1 = FMSY, lty = 3)

# Add axis label
mtext("Year", side = 1, outer = TRUE, line = 3)
dev.off()


## Partial F
# F-at-age
jpeg("./figures/f_at_age.jpeg", height = 210, width = 210, units = "mm", res = 350)
matplot(fit$data$years[-length(fit$data$years)], 
        faytable(fit)[-nrow(faytable(fit)),], 
        type = "l", 
        # col = 1:ncol(faytable(fit)), 
        col = cbPalette, 
        lty = 1:ncol(faytable(fit)), 
        lwd = 2, 
        xlab = "Year", 
        ylab = "F")
legend("topright", 
       col = cbPalette, 
       lty = 1:ncol(faytable(fit)), 
       legend = c(colnames(faytable(fit))[1:7],"7+"), 
       bty = 'n', 
       lwd = 2)
dev.off()


## SRR with CIs
library(plotrix)
colfunc <- colorRampPalette(c("#E69F00", "#0072B2"))


## CI for recruitment
jpeg("./figures/srr_random_walk_1981.jpeg", width = 15, height = 7.5, res = 300, units = "in")
par(mar = c(5.1, 4.1, 1.1, 1.1))
plot(fitsum$SSB/1e3, fitsum$Rec/1e3,
     type = "n",
     xlab = "SSB ('000 t)",
     ylab = "Recruitment Age 0 ('000s)",
     col = "#000000",
     xlim = c(2.38596, 115),
     ylim = c(0, 2300),
     lwd = 2)
segments(x0 = head(fitsum$SSB/1e3,-1),
         y0 = head(fitsum$Rec/1e3,-1),
         x1 = tail(fitsum$SSB/1e3,-1),
         y1 = tail(fitsum$Rec/1e3,-1),
         lwd = 2,
         col = colfunc(nrow(fitsum)-1))
points(fitsum$SSB/1e3, fitsum$Rec/1e3,
       col = colfunc(nrow(fitsum)),
       pch = 19)
plotCI(fitsum$SSB/1e3, fitsum$Rec/1e3, 
       ui = fitsum$RecHigh/1e3, 
       li = fitsum$RecLow/1e3,
       col = colfunc(nrow(fitsum)),
       add = TRUE)
text(fitsum$SSB/1e3, fitsum$Rec/1e3, 
     substr(as.character(fitsum$year), 3, 4), 
     col = colfunc(nrow(fitsum)),
     pos = 4)
dev.off()


## Retro zoomed in on last ten years
jpeg("./figures/retro_short.jpeg", 
     height = 160, width = 160, units = "mm", res = 350)
par(mfrow = c(3, 1), mar = c(0, 4.1, 0, 1.1), oma = c(5.1, 0, 1.1, 0))

# SSB
plot(fitsum$year, fitsum$SSB/1e3,
     xlab = "Year",
     ylab = "SSB ('000 t)",
     xlim = c((assessment.year-9), assessment.year),
     ylim = c((min(fitsum[as.character((assessment.year-9):assessment.year),"SSBLow"])/1e3),(max(fitsum[as.character((assessment.year-9):assessment.year),"SSBHigh"])/1e3)),
     axes = F,
     type = "n")
polygon(c(1981:assessment.year, rev(1981:assessment.year)), 
        c(fitsum$SSBLow, rev(fitsum$SSBHigh))/1e3,
        border = "#e7e7e7",
        col = "#e7e7e7")
lines(fitsum$year, fitsum$SSB/1e3,
      lwd = 2,
      col = "#000000")
# Peel 1
lines(rownames(summary(mretro[[1]])), summary(mretro[[1]])[,"SSB"]/1e3,
      lwd = 2,
      col = "#E69F00")
# Peel 2
lines(rownames(summary(mretro[[2]])), summary(mretro[[2]])[,"SSB"]/1e3,
      lwd = 2,
      col = "#56B4E9")
# Peel 3
lines(rownames(summary(mretro[[3]])), summary(mretro[[3]])[,"SSB"]/1e3,
      lwd = 2,
      col = "#009E73")
# Peel 4
lines(rownames(summary(mretro[[4]])), summary(mretro[[4]])[,"SSB"]/1e3,
      lwd = 2,
      col = "#F0E442")
# Peel 5
lines(rownames(summary(mretro[[5]])), summary(mretro[[5]])[,"SSB"]/1e3,
      lwd = 2,
      col = "#0072B2")
rhossb <- round(mohn(mretro)["SSB"],2)
legend("topleft", bty = "n", cex = 2, text.font = 2,
       legend = bquote(rho ~ "=" ~ .(rhossb)))
axis(2)
box()

# F
plot(fitsum$year, fitsum$Fbar,
     xlab = "Year",
     ylab = expression("F"["1-3"]),
     xlim = c((assessment.year-9), assessment.year),
     ylim = c((min(fitsum[as.character((data.year-9):data.year),"FLow"])),(max(fitsum[as.character((data.year-9):data.year),"FHigh"]))),
     axes = F,
     type = "n")
polygon(c(1981:data.year, rev(1981:data.year)), 
        c(fitsum[fitsum$year<assessment.year,"FLow"], rev(fitsum[fitsum$year<assessment.year,"FHigh"])),
        border = "#e7e7e7",
        col = "#e7e7e7")
lines(fitsum[fitsum$year<assessment.year,"year"], fitsum[fitsum$year<assessment.year,"Fbar"],
      lwd = 2,
      col = "#000000")
# Peel 1
lines(rownames(summary(mretro[[1]]))[-length(rownames(summary(mretro[[1]])))],
      summary(mretro[[1]])[-length(rownames(summary(mretro[[1]]))),"Fbar(1-3)"],
      lwd = 2,
      col = "#E69F00")
# Peel 2
lines(rownames(summary(mretro[[2]]))[-length(rownames(summary(mretro[[2]])))],
      summary(mretro[[2]])[-length(rownames(summary(mretro[[2]]))),"Fbar(1-3)"],
      lwd = 2,
      col = "#56B4E9")
# Peel 3
lines(rownames(summary(mretro[[3]]))[-length(rownames(summary(mretro[[3]])))],
      summary(mretro[[3]])[-length(rownames(summary(mretro[[3]]))),"Fbar(1-3)"],
      lwd = 2,
      col = "#009E73")
# Peel 4
lines(rownames(summary(mretro[[4]]))[-length(rownames(summary(mretro[[4]])))],
      summary(mretro[[4]])[-length(rownames(summary(mretro[[4]]))),"Fbar(1-3)"],
      lwd = 2,
      col = "#F0E442")
# Peel 5
lines(rownames(summary(mretro[[5]]))[-length(rownames(summary(mretro[[5]])))],
      summary(mretro[[5]])[-length(rownames(summary(mretro[[5]]))),"Fbar(1-3)"],
      lwd = 2,
      col = "#0072B2")
rhoF <- round(mohn(mretro, lag = 1)["Fbar(1-3)"],2)
legend("topleft", bty = "n", cex = 2, text.font = 2,
       legend = bquote(rho ~ "=" ~ .(rhoF)))
axis(2)
box()

# Recruitment
plot(fitsum$year, fitsum$Rec,
     xlab = "Year",
     ylab = "Recruitment (Age 0)",
     ylim = c((min(fitsum[as.character((data.year-9):data.year),"RecLow"])),(max(fitsum[as.character((data.year-9):data.year),"RecHigh"]))),
     xlim = c((assessment.year-9), assessment.year),
     type = "n")
polygon(c(1981:data.year, rev(1981:data.year)), 
        c(fitsum[fitsum$year<assessment.year,"RecLow"], rev(fitsum[fitsum$year<assessment.year,"RecHigh"])),
        border = "#e7e7e7",
        col = "#e7e7e7")
lines(fitsum[fitsum$year<assessment.year,"year"], fitsum[fitsum$year<assessment.year,"Rec"], 
      lwd = 2,
      col = "#000000")
# Peel 1
lines(rownames(summary(mretro[[1]]))[-length(rownames(summary(mretro[[1]])))],
      summary(mretro[[1]])[-length(rownames(summary(mretro[[1]]))),"R(age 0)"],
      lwd = 2,
      col = "#E69F00")
# Peel 2
lines(rownames(summary(mretro[[2]]))[-length(rownames(summary(mretro[[2]])))],
      summary(mretro[[2]])[-length(rownames(summary(mretro[[2]]))),"R(age 0)"],
      lwd = 2,
      col = "#56B4E9")
# Peel 3
lines(rownames(summary(mretro[[3]]))[-length(rownames(summary(mretro[[3]])))],
      summary(mretro[[3]])[-length(rownames(summary(mretro[[3]]))),"R(age 0)"],
      lwd = 2,
      col = "#009E73")
# Peel 4
lines(rownames(summary(mretro[[4]]))[-length(rownames(summary(mretro[[4]])))],
      summary(mretro[[4]])[-length(rownames(summary(mretro[[4]]))),"R(age 0)"],
      lwd = 2,
      col = "#F0E442")
# Peel 5
lines(rownames(summary(mretro[[5]]))[-length(rownames(summary(mretro[[5]])))],
      summary(mretro[[5]])[-length(rownames(summary(mretro[[5]]))),"R(age 0)"],
      lwd = 2,
      col = "#0072B2")
rhorec <- round(mohn(mretro, lag = 1)["R(age 0)"],2)
legend("topleft", bty = "n", cex = 2, text.font = 2,
       legend = bquote(rho ~ "=" ~ .(rhorec)))
box()
mtext("Year", side = 1, outer = TRUE, line = 3)
dev.off()
