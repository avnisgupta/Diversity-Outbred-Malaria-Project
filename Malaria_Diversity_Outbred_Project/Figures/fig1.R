#install.packages("plotrix")
library(plotrix)
library(tidyverse)
# 
head(Dataic)
ggplot(Dataic, aes( x= Day, y = RBC, group = Number, color = Description))+ geom_line() + xlim(c(0,15))
head(Dataic)
Grouped <- Dataic %>%
  subset(Description != "DO") %>%
  group_by(Description, Day) %>%
  summarize(RBCavg = mean(RBC,na.rm = T), RBCse = std.error(RBC, na.rm = T), Parasitemiaavg = mean(Parasitemia, na.rm = T), Parasitemiase = std.error(Parasitemia, na.rm = T), PercentChangeinWeightavg = mean(`Percent Change in Weight`,na.rm=T), PercentChangeinWeightse = std.error(`Percent Change in Weight`,na.rm = T), ChangeinTemperatureavg = mean(`Change in Temperature`, na.rm=T), ChangeinTemperaturese = std.error(`Change in Temperature`, na.rm=T), ParasiteDensityavg = mean(`Parasite Density`, na.rm=T), ParasiteDensityse = std.error(`Parasite Density`, na.rm = T) )

Grouped$Description <- factor(Grouped$Description, levels = c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HILtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ","DO"))

Dataic$Description <- factor(Dataic$Description, levels = c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ", "NZO/HILtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ","DO"))

strain_colors <- c("#FFDC00", "DarkGrey", "LightPink", "#0064C9", "#7FDBFF", "#2ECC40", "#FF4136", "#B10DC9", "LightGrey") 

#RBC
A <- ggplot(Grouped, aes(x = Day, y = RBCavg, color = Description))  + geom_ribbon(aes(ymax = RBCavg + RBCse, ymin = RBCavg- RBCse, fill = Description), color = NA, alpha =0.3) + geom_line( size = 0.8) + scale_color_manual(values = strain_colors) +  scale_fill_manual(values = strain_colors)+ labs(y = expression(atop(paste("10"^"6"," RBC"),"per ??L of Blood")), x = "")+ theme(legend.position = "none", axis.title.y = element_text(size = 12)) + scale_x_continuous(limits = c(0,15), expand = c(0,0)) +scale_y_continuous(limits = c(0,11), expand = c(0,0))


B <- ggplot(Grouped, aes(x = Day, y = RBCavg, color = Description)) + geom_line(data = subset(Dataic,Dataic$Description=="DO"),aes(x=Day,y=RBC, group = Number), color = "Light gray", alpha = 0.5) + geom_line( size = 0.8)+ scale_color_manual(values = strain_colors)  + labs(y = expression(atop(paste("10"^"6"," RBC"),"per ??L of Blood")), x = "")+ theme(legend.position = "none", axis.title.y = element_text(size = 12)) + scale_x_continuous(limits = c(0,15), expand = c(0,0)) +scale_y_continuous(limits = c(0,11), expand = c(0,0))


#PercentChangeinWeight
C <- ggplot(Grouped, aes(x = Day, y = PercentChangeinWeightavg, color = Description))  + geom_ribbon(aes(ymax = PercentChangeinWeightavg + PercentChangeinWeightse, ymin = PercentChangeinWeightavg- PercentChangeinWeightse, fill = Description), color = NA, alpha =0.3) + geom_line( size = 0.8) + scale_color_manual(values = strain_colors) +  scale_fill_manual(values = strain_colors) + labs(y = "Percent Change in Weight", x = "")+ theme(legend.position = "none", axis.title.y = element_text(size = 12))+ scale_x_continuous(limits = c(0,15), expand = c(0,0)) +scale_y_continuous(limits = c(-30,35), expand = c(0,0))

D <- ggplot(Grouped, aes(x = Day, y = PercentChangeinWeightavg, color = Description))  + geom_line(data = subset(Dataic,Dataic$Description=="DO"),aes(x=Day,y=`Percent Change in Weight`, group = Number),color = "Light gray", alpha = 0.5) + geom_line( size = 0.8)+ scale_color_manual(values = strain_colors)  + labs(y = "Percent Change in Weight", x = "")+ theme(legend.position = "none", axis.title.y = element_text(size = 12)) + scale_x_continuous(limits = c(0,15), expand = c(0,0)) +scale_y_continuous(limits = c(-30,35), expand = c(0,0))



#Change in Temperature
E <-  ggplot(Grouped, aes(x = Day, y = ChangeinTemperatureavg, color = Description))  + geom_ribbon(aes(ymax = ChangeinTemperatureavg + ChangeinTemperaturese, ymin = ChangeinTemperatureavg- ChangeinTemperaturese, fill = Description), color = NA, alpha =0.3) + geom_line( size = 0.8) + scale_color_manual(values = strain_colors) +  scale_fill_manual(values = strain_colors) + labs(y = "Change in Temperature (°C)", x = "")+theme(legend.position = "none", axis.title.y = element_text(size = 12)) + scale_x_continuous(limits = c(0,15), expand = c(0,0)) +scale_y_continuous(limits = c(-16,4), expand = c(0,0))

F <- ggplot(Grouped, aes(x = Day, y = ChangeinTemperatureavg, color = Description))  + geom_line(data = subset(Dataic,Dataic$Description=="DO"),aes(x=Day,y=`Change in Temperature`, group = Number),color = "Light gray", alpha = 0.5) + geom_line( size = 0.8)+ scale_color_manual(values = strain_colors)   + labs(y = "Change in Temperature (°C)", x = "")+ theme(legend.position = "none", axis.title.y = element_text(size = 12)) + scale_x_continuous(limits = c(0,15), expand = c(0,0)) +scale_y_continuous(limits = c(-16,4), expand = c(0,0))


#Parasitemia
G <- ggplot(Grouped, aes(x = Day, y = Parasitemiaavg, color = Description))  + geom_ribbon(aes(ymax = Parasitemiaavg + Parasitemiase, ymin = Parasitemiaavg- Parasitemiase, fill = Description), color = NA, alpha =0.3) + geom_line( size = 0.8) + scale_color_manual(values = strain_colors) +  scale_fill_manual(values = expression(atop("Parasitemia",paste("(Infected RBC per Total RBC)"))), x = "") + theme(legend.position = "none", axis.title.y = element_text(size = 12))  + scale_x_continuous(limits = c(0,15), expand = c(0,0)) +scale_y_continuous(limits = c(0,0.6), expand = c(0,0))

H <- ggplot(Grouped, aes(x = Day, y = Parasitemiaavg, color = Description)) + geom_line(data = subset(Dataic,Dataic$Description=="DO"),aes(x=Day,y=Parasitemia, group = Number), color = "Light gray", alpha = 0.5) + geom_line( size = 0.8)+ scale_color_manual(values = strain_colors)  + labs(y = expression(atop("Parasitemia",paste("(Infected RBC per Total RBC)"))), x = "")+ theme(legend.position = "none", axis.title.y = element_text(size = 12)) + scale_x_continuous(limits = c(0,15), expand = c(0,0)) +scale_y_continuous(limits = c(0,0.6), expand = c(0,0))


#Parasite Density 
I <- ggplot(Grouped, aes(x = Day, y = ParasiteDensityavg, color = Description)) + geom_ribbon(aes(ymax = ParasiteDensityavg + ParasiteDensityse, ymin = ParasiteDensityavg- ParasiteDensityse, fill = Description), color = NA, alpha =0.3) + geom_line( size = 0.8) + scale_color_manual(values = strain_colors) +  scale_fill_manual(values = strain_colors) + labs(y = expression(atop("Parasite Density",paste("(10"^"6"," Parasites per ??L of Blood)"))), x = "Day Post Infection")  + theme(legend.position = "none", axis.title.y = element_text(size = 12)) + scale_x_continuous(limits = c(0,15), expand = c(0,0)) +scale_y_continuous(limits = c(0,2.5), expand = c(0,0))


J <- ggplot(Grouped, aes(x = Day, y = ParasiteDensityavg, color = Description)) + geom_line(data = subset(Dataic,Dataic$Description=="DO"),aes(x=Day,y=`Parasite Density`, group = Number), color = "Light gray", alpha = 0.5) + geom_line( size = 0.8)+ scale_color_manual(values = strain_colors, name = "Mouse Strain")  + labs(y = expression(atop("Parasite Density",paste("(10"^"6"," Parasites per ??L of Blood)"))), x = "Day Post Infection")+ scale_x_continuous(limits = c(0,15), expand = c(0,0)) +scale_y_continuous(limits = c(0,2.5), expand = c(0,0))+ theme(axis.title.y = element_text(size = 12), legend.position = "none") 

library(ggpubr)

p <- ggarrange(A,B,C,D,E,F,G,H,I,J, labels = c("A","B","C","D","E","F","G","H","I","J"), nrow = 5, ncol = 2, align = "v")


#install.packages("svglite")
library(svglite)
ggsave("Figure2", device = "pdf", plot = p, height = 14, width = 12,units = "in" )

plegend <- ggplot(Dataic, aes(x = Day, y = RBC, group = Number , color = Description)) + geom_line(size = 5) + scale_color_manual(values = strain_colors, name = "Mouse Strains") + theme(legend.position = "bottom")

ggsave("Figure2Legend", device = "pdf", plot = plegend, height = 14, width = 12,units = "in")

