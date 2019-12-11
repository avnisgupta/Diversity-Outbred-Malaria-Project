# Figure 1A Survival Curve

library(survminer)
library(RTCGA.clinical)
library(survival)
library(tidyverse)
library(readxl)
library(kableExtra)
library(ggpubr)


library(readr)
DO  <- read_csv("DO Experiment Data Analysis Version.csv")


Measurements <- tibble(Number = as.data.frame(table(DO$Number))[,1], Count = as.data.frame(table(DO$Number))[,2])

names(Measurements) <- c("Number", "Measurements")

uDO <- distinct(DO,Number, .keep_all = TRUE)


batchDO <- dplyr::select(uDO, Batch, Number, Description, `Lived?`)

batchMeasurements <-  merge(Measurements, batchDO)
#We need columns Day, Mouse ID, Mouse.Genome, Status
#Made a vector of the number of mice
Mouse.ID <- batchMeasurements$Number
batchMeasurements$Measurements[batchMeasurements$Measurements==14] <- 13
batchMeasurements$Measurements[batchMeasurements$Measurements==20] <- 13
batchMeasurements$Measurements[batchMeasurements$Measurements==23] <- 13

#Apply and adjust for 0 to 5 day sampling gap
Death.Day <- batchMeasurements$Measurements + 3
#Bring the genotype to the mouse ID
Mouse.Genotype <- batchMeasurements$Description
Mouse.Status <- ifelse(Death.Day <15, 1,0)

ACODsurv <- data.frame(Death.Day, Mouse.ID, Mouse.Status, Mouse.Genotype)

fit <- survfit(Surv(Death.Day, Mouse.Status) ~ Mouse.Genotype,
               data = ACODsurv)
# Visualize with survminer
ggsurvplot(fit, data = ACODsurv, size = 1.5, xlim = c(0,15), ylim = c(0,1.01), legend = "none", axes.offset = F, palette = c("#f4b5c0", "#f7db14", "#696868", "#2c9d36", "lightgrey", "#1861a9", "#82cfe2", "#f94136", "#b41183"), ggtheme = theme_classic(base_size = 18), size =1) + ylab("Percent Survival") +xlab("Days") 
