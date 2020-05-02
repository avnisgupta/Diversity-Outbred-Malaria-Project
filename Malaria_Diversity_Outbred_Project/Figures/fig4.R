

library(readr)
Data <- read_csv("GitHub/Diversity-Outbred-Malaria-Project/Malaria_Diversity_Outbred_Project/Setup/Creating_Phenotype_File/")

Data$Sex <- "F"


# Making an initial weight and temperature dataframe

#This converts columns we care about to numeric
id <- c("Temperature","Weight", "Accuri Count")
Data[id] = data.matrix(Data[id])

Day0 <- Data[Data$Day==0,]

Initial <- subset(Day0, select = c("Number", "Temperature", "Weight"))
names(Initial) <- c("Number", "Initial Temperature", "Initial Weight")
Datai <-  merge(Data,Initial)


# Computing Mouse physiology parameters

Dataic <- Datai %>%
  # calculate percent weight relative to day 0
  mutate(`Percent Bodyweight` = ((Datai$Weight) / (Datai$`Initial Weight`) * 100))  %>%
  # calculate percent weight loss relative to day 0
  mutate(`Percent Change in Weight` = (((Datai$Weight)-(Datai$`Initial Weight`)) / (Datai$`Initial Weight`)) * 100)%>%
  # calculate weight loss relative to day 0
  mutate(`Change in Weight` = (Datai$Weight)-(Datai$`Initial Weight`)) %>%
  # calculate RBCs from raw accuri score
  mutate(RBC = (Datai$`Accuri Count`) / 2000) %>%  
  # calculate delta temperature relative to day 0
  mutate(`Change in Temperature` = ((Datai$Temperature) - (Datai$`Initial Temperature`))) %>%
  # calculate parasite density 
  mutate(`Parasite Density` = Parasitemia * RBC)

# Removing a measurement error
Dataic$RBC[2710] <- NA
Dataic %>%
  subset(Number == 210 & Day == 0) %>%
  select(RBC) 

# Summarizing Statistics for QTL Analysis

Phenotypes <- Dataic %>% 
  subset(Day %in% c(0,5:15)) %>%
  group_by(Number) %>%
  summarize(White = as.numeric(`Coat Color`[Day==0]== "White"),
            Black = as.numeric(`Coat Color`[Day==0]== "Black"),
            Star = `Star?`[Day==0],
            Lived = `Lived?`[Day==0],
            MinDeltaTemp = min(`Change in Temperature`, na.rm = TRUE),
            Day0Weight= `Initial Weight`[Day==0],
            MinPercWeight = min(`Percent Change in Weight`, na.rm=TRUE),
            MinDeltaWeight = min(`Change in Weight`, na.rm = TRUE),
            #This is likely just a measure of age since only young mice are this phenotype
            #Gainer = as.numeric(MinDeltaWeight==0), 
            MinRBC= min(RBC, na.rm=TRUE),
            MaxPD = max(`Parasite Density`, na.rm=TRUE),
            MaxParasitemia = max(Parasitemia, na.rm = TRUE),
            Description = Description[1])

# Changing experimental measurements from a mouse (415) that died before the experiment 
Phenotypes[397,6:13] <- NA

Phenotypes$MinDeltaTemp[Phenotypes$Lived==0] <- NA


axiom8_phenotypes <- read_csv("axiom7/axiom8_phenotypes.csv", skip = 3)

axiom8_Par <- axiom8_phenotypes %>% select(Number,ParasitemiaClusters) 
Together <-merge(Dataic, axiom8_Par, by = "Number")
InterestTogether <- subset(Together,is.na(Together$ParasitemiaClusters)==F)
InterestTogether$ParasitemiaClusters <- as.factor(InterestTogether$ParasitemiaClusters)
levels(InterestTogether$ParasitemiaClusters) <- c("Moderate", "Mild", "Severe")
InterestTogether$ParasitemiaClusters <-factor(InterestTogether$ParasitemiaClusters, levels = c("Mild", "Moderate", "Severe"))
ggplot(InterestTogether, aes(x = Day, y = Parasitemia, color = ParasitemiaClusters, group = Number)) + scale_x_continuous(limits=c(0,15), expand = c(0,0))+ scale_y_continuous(expand = c(0,0))+ theme_classic(base_size = 17) + geom_line(size = 1, alpha =0.7) + facet_wrap(~ParasitemiaClusters) + theme(legend.position = "none") + scale_color_manual(values = c("#9dab86","#cca176","#dd7162")) +labs(y = expression(atop("Parasitemia",paste("(Infected RBC per Total RBC)"))), x = "Day Post Infection")


sum(axiom8_phenotypes$UndectedParasitemia,na.rm = T)
Eachmouse <-subset(InterestTogether, InterestTogether$Day==0)
table(Eachmouse$ParasitemiaClusters)
83+379+27
subset(Together,is.na(Together$ParasitemiaClusters)==T)
library(qtl2)
library(tidyverse)


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

pheno <- "ParasitemiaCluster2"
model <- "normal" #c("normal", "binary")
perms <- 1


load(paste(paste(cross_name, sep="_", pheno), sep=".", "Rda"))

paste(paste(cross_name, sep="_", pheno), sep=".", "Rda")
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
     cross$pmap,col="black", altcol = "black", bgcol = "white", altbgcol = "white", gap = 0)



##Add horizontal threshold lines to plot
abline(h = sig99_cross_outpheno, col="#FF0000FF")
abline(h = sig95_cross_outpheno, col="#FFAA00FF")
abline(h = sig90_cross_outpheno, col="#FFAA00FF", lty=5)


## Parental Contribution 


plot_coefCC(coef_cross_outpheno,
            map = cross$pmap[chr_cross_outpheno],
            columns = 1:8,
            scan1_output = cross_outpheno,
            add = TRUE,
            xlim =NULL,
            ylim = NULL,
            bgcolor = "white", legend = "topright")

chr_cross_outpheno
cross_outpheno
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

k <- calc_kinship(axiom7ex_apr, "loco")
query_variants <-create_variant_query_func("../Desktop/cc_variants.sqlite")
query_genes <- create_gene_query_func("../Desktop/mouse_genes_mgi.sqlite")
# This works but the LOD doesn't match the previous scan so obviously I'm fe
# This works but the LOD doesn't match the previous scan so obviously I'm feeding in something wrong. Actually I'm not sure I'm feeding it anything from the actual phenotype.
out_snps <- scan1snps(axiom8ex_pr, axiom8ex$pmap, axiom8ex$pheno, k[["13"]], axiom7ex_sex, query_func=query_variants,chr= 13, start = 54.07724, end = 54.17411, keep_all_snps=TRUE)

genes <- query_genes(chr= 13, start = 54.07724, end = 54.17411)

## Plotting Interval

par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out_snps$lod, out_snps$snpinfo, drop_hilit=1.5, genes=genes)

top <- top_snps(out_snps$lod, out_snps$snpinfo)
print(top[,c(1, 8:15, 20)], row.names=FALSE)
