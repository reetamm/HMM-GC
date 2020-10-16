library(ggplot2)
library(forecast)
library(lubridate)
library(gdata)
library(tidyr)
library(psych)
theme_set(theme_bw())
library(sf)

simFile <- 'sim_hmm_cbay_s6g3'              #MVNHMM simulation output
copulaSimFile <- 'gammacopula'              #File generated from gammacorrels.R
obsFile <- 'cbayJulSep0019'                 #Observation data in MVNHMM format
date <- rep(seq(ymd('2000-07-01'), ymd('2000-09-30'), by = 'days'), 20)
month <- month(date)
day <- day(date)
year <- rep(2000:2019, each = 92)
date <- make_date(year, month, day)
date <- substr(date, 1, 7)
head(date)

obsdata <- read.table(obsFile)
hmmoutput <- read.table(simFile)
copulaoutput <- read.table(copulaSimFile, header = T)

month <- rep(c(7,8,9),20)

tmp1 <- aggregate(obsdata, list(date), sum)
tmp2 <- aggregate(copulaoutput, list(date), sum)
tmp1 <- tmp1[,-1]
tmp2 <- tmp2[,-1]
tmp1 <- aggregate(tmp1, list(month), mean)
tmp2 <- aggregate(tmp2, list(month), mean)
input.monthly.mean <- tmp1[,-1]
copula.monthly.mean <- tmp2[,-1]
tmp2 <- c(t(as.matrix(copula.monthly.mean)))
tmp1 <- c(t(as.matrix(input.monthly.mean)))


tmp3 <- aggregate(obsdata, list(date), function(x)mean(x == 0))
tmp4 <- aggregate(copulaoutput, list(date), function(x)mean(x == 0))
tmp3 <- tmp3[, -1]
tmp4 <- tmp4[, -1]
tmp3 <- aggregate(tmp3, list(month), mean)
tmp4 <- aggregate(tmp4, list(month), mean)
input.monthly.prop <- tmp3[,-1]
copula.monthly.prop <- tmp4[,-1]
tmp4 <- c(t(as.matrix(copula.monthly.prop)))
tmp3 <- c(t(as.matrix(input.monthly.prop)))

summaryplot <- data.frame(month = rep(month.abb[7:9], each = 1927), mean.obs = tmp1, mean.sim = tmp2, prop.obs = tmp3, prop.sim = tmp4)
summaryplot$month <- relevel(summaryplot$month, ref = 'Jul')
RMSE(summaryplot$mean.sim,summaryplot$mean.obs) #plug this value into the next plot
xrng <- range(summaryplot$mean.obs)
yrng <- range(summaryplot$mean.sim)
ggplot(summaryplot,aes(x = mean.obs, y = mean.sim, col = month)) + geom_point() + geom_abline() +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 15), 
      axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
      axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
      axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = 2, face = "plain"),
      legend.text = element_text(size=12), legend.title = element_text(size=15)) +
      guides(colour = guide_legend(override.aes = list(size=3))) +
      labs(x = 'IMERG precipitation (mm)', y = 'Synthetic precipitation from HMM-GC (mm)', col = 'Month') +
      annotate(geom = "text", x = xrng[1], y = yrng[2], label = 'RMSE = 11.69 mm', hjust = 0, vjust = 1, size = 7)

RMSE(summaryplot$prop.sim,summaryplot$prop.obs) #plug this value into the next plot
xrng <- range(summaryplot$prop.obs)
yrng <- range(summaryplot$prop.sim)
ggplot(summaryplot, aes(x = prop.obs,y = prop.sim,col = month)) + geom_point() + geom_abline() +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 15), 
      axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
      axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
      axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
      axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = 2, face = "plain"),
      legend.text = element_text(size=12), legend.title = element_text(size=15)) +
      guides(colour = guide_legend(override.aes = list(size=3))) +
      labs(x = 'Proportion of dry days based on IMERG', y = 'Proportion of dry days based on HMM', col = 'Month') +
      annotate(geom = "text", x = xrng[1], y = yrng[2], label = 'RMSE = 0.03', hjust = 0, vjust = 1, size = 7)

summaryplot2 <- gather(summaryplot[, 1:3], source, mean, mean.obs:mean.sim)
summaryplot3 <- gather(summaryplot[, c(1,4,5)], source, prop, prop.obs:prop.sim)

ggplot(summaryplot2, aes(x = month, y = mean, fill = source)) + geom_boxplot()
ggplot(summaryplot3, aes(x = month, y = prop, fill = source)) + geom_boxplot()