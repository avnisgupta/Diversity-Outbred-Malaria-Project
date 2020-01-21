# Set the directory in Malaria_Diversity_Project
setwd("C:/Users/Schneider Lab/Desktop/Malaria_Diversity_Outbred_Project")

# Load Packages
library(readr)
library(tidyverse)
library(ggpubr)
Data <- read_csv("Setup/Raw_Data/axiom8_RawData.csv")

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

Lived <- Dataic %>%
  subset(`Lived?`==1)
Died <- Dataic %>%
  subset(`Lived?`==0)

dt <- c(-16,4)
pa <- c(0,0.8)
dw <- c(-40,35)
rb <- c(0, 12)

a <- ggplot(Lived, aes(y = `Change in Temperature`, x = Parasitemia, group = Number)) + geom_path(color = "royalblue", alpha=0.7) + geom_path(data = Died, color = "salmon") + theme_classic()+ coord_cartesian(xlim = pa, ylim = dt, expand = FALSE)

b <- ggplot(Lived, aes(y = `Change in Temperature`, x = RBC, group = Number)) + geom_path(color = "royalblue", alpha=0.7) + geom_path(data = Died, color = "salmon") + theme_classic()+ coord_cartesian(xlim = rb, ylim = dt, expand = FALSE)

c <- ggplot(Lived, aes(y = `Change in Temperature`, x = `Percent Change in Weight`, group = Number)) + geom_path(color = "royalblue", alpha=0.7) + geom_path(data = Died, color = "salmon") + theme_classic()+ coord_cartesian(xlim = dw, ylim = dt, expand = FALSE)

d <- ggplot(Lived, aes(y = Parasitemia, x = RBC, group = Number)) + geom_path(color = "royalblue", alpha=0.7) + geom_path(data = Died, color = "salmon") + theme_classic()+ coord_cartesian(xlim = rb, ylim = pa, expand = FALSE)

e <- ggplot(Lived, aes(y = Parasitemia, x = `Change in Temperature`, group = Number)) + geom_path(color = "royalblue", alpha=0.7) + geom_path(data = Died, color = "salmon") + theme_classic()+ coord_cartesian(xlim = dt, ylim = pa, expand = FALSE)

f <- ggplot(Lived, aes(y = Parasitemia, x = `Percent Change in Weight`, group = Number)) + geom_path(color = "royalblue", alpha=0.7) + geom_path(data = Died, color = "salmon") + theme_classic()+ coord_cartesian(xlim = dw, ylim = pa, expand = FALSE)

g <- ggplot(Lived, aes(y = RBC, x = Parasitemia, group = Number)) + geom_path(color = "royalblue", alpha=0.7) + geom_path(data = Died, color = "salmon") + theme_classic()+ coord_cartesian(xlim = pa, ylim = rb, expand = FALSE)

h <- ggplot(Lived, aes(y = RBC, x = `Change in Temperature`, group = Number)) + geom_path(color = "royalblue", alpha=0.7) + geom_path(data = Died, color = "salmon") + theme_classic()+ coord_cartesian(xlim = dt, ylim = rb, expand = FALSE)

i <- ggplot(Lived, aes(y = RBC, x = `Percent Change in Weight`, group = Number)) + geom_path(color = "royalblue", alpha=0.7) + geom_path(data = Died, color = "salmon") + theme_classic()+ coord_cartesian(xlim = dw, ylim = rb, expand = FALSE)


ggarrange(a,b,c,d,e,f,g,h,i, ncol = 3, nrow = 3)

OrderedPhenotypes <- Phenotypes %>%
  arrange(-MinDeltaTemp) %>%
  mutate(Index = as.numeric(rownames(OrderedPhenotypes)))


ggplot(OrderedPhenotypes, aes(x = Index, y = MinDeltaTemp, color = as.factor(Lived))) + geom_point() + theme_classic(base_size = 18) + scale_color_manual(values =c("salmon","royalblue")) + theme(legend.position = "none") 
