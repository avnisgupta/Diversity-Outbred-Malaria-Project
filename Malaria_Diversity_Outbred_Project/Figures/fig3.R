# Figure 2 Survival Curve

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

library(qtl2)
library(tidyverse)

setwd("/Volumes/Adam/DO Project")
load("axiom8ex.Rda")
load("axiom7ex_apr.Rda")
load("axiom7ex_kinship.Rda")
load("axiom7ex_sex.Rda")
load("axiom7ex_Xcovar.Rda")


axiom8ex_apr <-axiom7ex_apr
axiom8ex_kinship <-axiom7ex_kinship
axiom8ex_sex <- axiom7ex_sex
axiom8ex_Xcovar <-axiom7ex_Xcovar


cross <- axiom8ex
cross_name <- "axiom8ex"

# Lived

## QTL

pheno <- "Died"
model <- "binary" #c("normal", "binary")
perms <- 1000


load(paste(paste(cross_name, sep="_", pheno), sep=".", "Rda"))

load(paste(paste("perm", cross_name, sep="_", pheno), sep=".", "Rda"))

summary(perm_cross_outpheno, alpha= c(0.2, 0.1, 0.05, 0.01))
sig99_cross_outpheno <- summary(perm_cross_outpheno, alpha= 0.01)
sig95_cross_outpheno <- summary(perm_cross_outpheno, alpha= 0.05)
sig90_cross_outpheno <- summary(perm_cross_outpheno, alpha= 0.1)
sig80_cross_outpheno <- summary(perm_cross_outpheno, alpha= 0.2)
bestmarker_cross_outpheno <- rownames (max(cross_outpheno, cross$pmap))
markerpos_bestmarker_cross_outpheno <- find_markerpos(cross, bestmarker_cross_outpheno)
chr_cross_outpheno <- as.numeric(markerpos_bestmarker_cross_outpheno[1,1])

markerpos_bestmarker_cross_outpheno

threshold_cross_outpheno <-maxlod(cross_outpheno, map = NULL, chr = NULL)-1

#Find peaks above significance thresholds
peaks_cross_outpheno <- find_peaks(scan1_output = cross_outpheno,
                                   map= cross$pmap,
                                   threshold = threshold_cross_outpheno,
                                   peakdrop = Inf,
                                   drop = NULL,
                                   prob = NULL,
                                   thresholdX = NULL,
                                   peakdropX = NULL,
                                   dropX = NULL,
                                   probX = NULL,
                                   expand2markers = TRUE,
                                   sort_by = "lod",
                                   cores = 1)


##To estimate the QTL effects from each of the 8 founder strains on the chromosome with the highest QTL
load(paste(paste("coef", cross_name, sep="_", pheno), sep=".", "Rda"))


##Plot qtl graph
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(cross_outpheno,
     cross$pmap,col="grey61", altcol = "grey61", bgcol = "white", altbgcol = "white", gap = 0)



##Add horizontal threshold lines to plot
abline(h = sig99_cross_outpheno, col="#FF0000FF")
abline(h = sig95_cross_outpheno, col="#FFAA00FF")
abline(h = sig90_cross_outpheno, col="#FFFF80FF")



## Parental Contribution 

plot_coefCC(coef_cross_outpheno,
            map = cross$pmap[chr_cross_outpheno],
            columns = 1:8,
            scan1_output = cross_outpheno,
            add = TRUE,
            xlim =NULL,
            ylim = NULL,
            bgcolor = "white", legend = "topright")


##Calculate LOD support intervals using lod_int()-(qtl2scan)
lod_int_bestmarker_cross_outpheno <- lod_int(scan1_output = cross_outpheno,
                                             map= cross$pmap,
                                             chr= chr_cross_outpheno,
                                             lodcolumn = 1,
                                             threshold = threshold_cross_outpheno,
                                             peakdrop = Inf,
                                             drop = 1.5,
                                             expand2markers = TRUE)

lod_int_bestmarker_cross_outpheno


## Genes of interest

peak_Mbp <- lod_int_bestmarker_cross_outpheno[2]
k <- calc_kinship(axiom7ex_apr, "loco")
query_variants <-create_variant_query_func("../Desktop/cc_variants.sqlite")
query_genes <- create_gene_query_func("../Desktop/mouse_genes_mgi.sqlite")
# This works but the LOD doesn't match the previous scan so obviously I'm feeding in something wrong. Actually I'm not sure I'm feeding it anything from the actual phenotype.
out_snps <- scan1snps(axiom8ex_pr, axiom8ex$pmap, axiom8ex$pheno, k[["10"]], axiom7ex_sex, query_func=query_variants,chr= 10, start=41.39297, end=43.77759, keep_all_snps=TRUE)


genes <- query_genes(chr = 10,start=41.39297, end=43.77759)



## Plotting Interval

par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out_snps$lod, out_snps$snpinfo, drop_hilit=1.5, genes=genes)

top <- top_snps(out_snps$lod, out_snps$snpinfo)
