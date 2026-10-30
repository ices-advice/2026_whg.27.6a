catch.year <- 2025
assessment.year <- 2026

# CVD-friendly palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", 
                         "#009E73", "#0072B2", 
                         "#D55E00", "#CC79A7")
                         
## Discards proportions-at-age
dat <- read.table("../data/lf.dat", skip = 5)
dat <- 1-dat


# Plot discards-at-age
jpeg("../figures/discard_proportions_at_age_lowprod.jpeg", 
     height = 150, width = 300, units = "mm", res = 350)
plot(1981:catch.year, dat$V2,
     ylim = c(0, 1),
     xlab = "Year",
     ylab = "Discard Proportion (numbers-at-age)",
     type = "l")
points(1981:catch.year, dat$V2,
       pch = "1")
for(i in 3:ncol(dat)){
  lines(1981:catch.year, dat[,i],
        col = cbbPalette[i-1])
  points(1981:catch.year, dat[,i],
         pch = as.character(i-1),
         col = cbbPalette[i-1])
}
dev.off()


## Recruitment
list.files("../SAM_output/")
summ <- read.table("../SAM_output/summary.txt", header = TRUE)

jpeg("../figures/recruitment_assumption_lowprod.jpeg", 
     height = 150, width = 210, units = "mm", res = 350)
plot(summ$year, summ$Rec, 
     type = "n",
     ylim = c(8812.92, 2200000),
     xlab = "Year",
     ylab = "Recruitment (Age 0)")
polygon(c(summ$year, rev(summ$year)), 
        c(summ$RecLow, rev(summ$RecHigh)),
        border = "#ffebbf",
        col = "#ffebbf")
lines(summ$year, summ$Rec,
      lwd = 2,
      col = "#E69F00")
segments(x0 = assessment.year-6, x1 = assessment.year-1,
         y1 = median(summ$Rec[summ$year %in% (assessment.year-6):(assessment.year-1)]),
         y0 = median(summ$Rec[summ$year %in% (assessment.year-6):(assessment.year-1)]),
         lwd = 3,
         lty = 3)
legend("topright",
       bty = "n",
       lty = 3,
       lwd = 3,
       legend = paste0("Median Recruitment ", 
                       assessment.year-6, "-", assessment.year-1))
dev.off()


## Survey Indices
library(stockassessment)
surveys <- read.ices("../data/survey.dat")
cols <- colorRampPalette(c("#ee1408", "#635143"))
range01 <- function(x){
  (x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))
}

# Visualise
jpeg("../figures/indices_at_age_lowprod.jpeg", height = 120, width = 300, units = "mm", res = 350)
par(mar = c(0, 0, 0, 0), oma = c(5.1, 4.1, 1.1, 1.1), mfrow = c(4, 2))
plot(1996:catch.year, range01(surveys$SWC.Q4[,1]),
     xlim = c(1985, assessment.year),
     ylim = c(0, 1),
     col = "#E69F00",
     xlab = "",
     ylab = "",
     lwd = 2,
     pch = "0",
     axes = FALSE,
     type = "b")
legend("left",
       legend = c("Old Q1", "New Q1", "Q4"),
       col = c("#000000", "#56B4E9", "#E69F00"),
       bty = "n",
       cex = 1.5,
       lwd = 2)
axis(2)
box()
for(i in 2:ncol(surveys$SWC.Q4)){
  plot(1996:catch.year, range01(surveys$SWC.Q4[,i]), 
       xlim = c(1985, assessment.year),
       type = "b", 
       pch = as.character(i-1),
       axes = FALSE,
       lwd = 2,
       col = "#E69F00")
  if(i %in% c(2:7)){
    # xx <- cbind(2011:assessment.year, range01(surveys$SCO.Q1[1:nrow(surveys$SCO.Q1),i-1]))
    # xx <- xx[complete.cases(xx),]
    lines(2011:assessment.year, range01(surveys$SCO.Q1[1:nrow(surveys$SCO.Q1),i-1]),
          type = "b",
          col = "#56B4E9",
          lwd = 2,
          pch = as.character(i-1))
  }
  if(i %in% c(2:7)){
    lines(1985:2010, range01(surveys$WCIBTS.Q1[,i-1]),
          type = "b",
          col = "#000000",
          lwd = 2,
          pch = as.character(i-1))
  }
  if(i %in% c(3, 5, 7)){
    axis(2)
  }
  if(i %in% c(7, 8)){
    axis(1)
  }
  box()
}
mtext("Year", side = 1, outer = TRUE, line = 3)
mtext("Standardised Index", side = 2, outer = TRUE, line = 3)
dev.off()


## Catch weights-at-age by gear/country
cdat <- read.csv(paste("catch_weights_2011-",catch.year,".csv",sep=''))
names(cdat) <- tolower(names(cdat))
cdat[cdat == 0] <- NA

def120_lnd <- cdat[cdat$fleets == "OTB_DEF_>=120_0_0_all" & cdat$catch.cat. == "Landings",]
def120_lnd_IRE <- def120_lnd[def120_lnd$country == "Ireland",]
def100 <- cdat[cdat$fleets == "OTB_DEF_100-119_0_0_all",]
defCRU <- cdat[cdat$fleets == "OTB_CRU_70-99_0_0_all",]

def120_lnd_IRE <- aggregate(def120_lnd_IRE[,7:14], list(def120_lnd_IRE$year), 
                            mean, na.rm = TRUE)


jpeg("../figures/TR1_LND_WTS_lowprod.jpeg", height = 150, width = 250, 
     units = "mm", res = 350)
par(mar = c(5.1, 4.1, 1.1, 1.1))

plot(unique(def120_lnd$year[def120_lnd$country == "UK(Scotland)"]),
     def120_lnd$age.0[def120_lnd$country == "UK(Scotland)"],
     col = cbbPalette[5],
     ylim = c(range(def120_lnd[,7:17], na.rm = TRUE)[1], 1250),
     pch = "0",
     xlab = "Year",
     ylab = "Mean Weight (g)")

for(i in 8:14){
  lines(unique(def120_lnd$year[def120_lnd$country == "UK(Scotland)"]),
        def120_lnd[def120_lnd$country == "UK(Scotland)", i],
        col = cbbPalette[5],
        lwd = 2,
        type = "b",
        pch = as.character(i-7))
}


for(i in 2:9){
  lines(def120_lnd_IRE$Group.1+0.25, def120_lnd_IRE[,i],
        col = cbbPalette[6],
        pch = as.character(i-2),
        type = "b",
        lty = 3)
  
}
dev.off()
