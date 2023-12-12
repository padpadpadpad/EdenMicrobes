#eden_microclimate
#Rainforest
microclim<-read.csv("/Users/jd768/Documents/PROJECTS/eDNA_PROJECT/microclimate/Full_Data_RBMB.csv")

library(ggplot2)
microclim_plot_temp <- 
  
  ggplot(microclim, aes(x=Habitat, y=temperature.C., color=Rep))+theme_classic()+
  geom_boxplot()

##ANOVA on DIVERSITY
anova1 <- aov(temperature.C.~Habitat, data = microclim)
summary(anova1)
TukeyHSD(anova1)

microclim_plot_moist <- 
  
  ggplot(microclim, aes(x=Habitat, y=volumetricwatercontent..., color=Rep))+theme_classic()+
  geom_boxplot()

microclim_plot_temp <- 
  
  ggplot(microclim, aes(x=Habitat, y=temperature.C., color=Rep))+theme_classic()+
  geom_point()+
  geom_jitter()


microclim_plot_moist <- 
  
  ggplot(microclim, aes(x=Habitat, y=volumetricwatercontent..., color=Rep))+theme_classic()+
  geom_point()


microclim_plot_temp <- 
  
  ggplot(microclim, aes(x=Habitat, y=temperature.C., color=Rep))+theme_classic()+
  geom_point()+
  geom_jitter()


microclim_plot_moist <- 
  
  ggplot(microclim, aes(x=Habitat, y=volumetricwatercontent..., color=Rep))+theme_classic()+
  geom_point()


microclim_plot_temp <- 
  
  ggplot(microclim, aes(x=Sample_Code, y=temperature.C., color=Rep))+theme_classic()+
  geom_boxplot()

microclim_plot_temp+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


microclim_plot_moist <- 
  
  ggplot(microclim, aes(x=Sample_Code, y=volumetricwatercontent..., color=Rep))+theme_classic()+
  geom_boxplot()

microclim_plot_moist+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

microclim2<- microclim
microclim2$Plot_rep <-paste(microclim2$Sample_Code,"-",microclim2$Rep)

microclim_plot_temp <- 
  
  ggplot(microclim2, aes(x=Plot_rep, y=temperature.C., color=Rep))+theme_classic()+
  geom_boxplot()

microclim_plot_temp+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


microclim_plot_moist <- 
  
  ggplot(microclim2, aes(x=Plot_rep, y=volumetricwatercontent..., color=Rep))+theme_classic()+
  geom_boxplot()

microclim_plot_moist+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# generate hourly averages per plot

library(dplyr)
microclim$DT<- paste(microclim$Date, " ", microclim$Time)

microclim$DT2 <- strptime(microclim$DT, "%d/%m/%Y %H:%M")

microclim2<-microclim

microclim2$group <- cut(microclim2$DT2, breaks="hour")
microclim2$group2<- paste(microclim2$Sample_Code, "-", microclim2$group)
setwd("/Users/jd768/Documents/PROJECTS/eDNA_PROJECT/microclimate/data/data_for_Euclidean/")
write.csv(microclim2,"microclim_hourlyC.csv")
microclim<-read.csv("microclim_hourlyB.csv")



library(data.table)
microclim3<-as.data.table(cbind(microclim2$group2,microclim2$temperature.C.,microclim2$volumetricwatercontent...))
microclim3$V2<-as.numeric(microclim3$V2)
microclim3$V3<-as.numeric(microclim3$V3)

microclim4<-aggregate(microclim3[, 2:3], list(microclim3$V1), mean)
#this option avoids the NA issues, and calculates mean of remaining groups
microclim4<-aggregate(microclim3[, 2:3], list(microclim3$V1), mean,na.rm = TRUE)
plot(microclim4$V2,microclim4$V3)

write.csv(microclim4,"microclim_hourly3.csv")
microclim5<-read.csv("microclim_hourly3b.csv")
table(microclim5$SampleSample)

setwd("/Users/jd768/Documents/PROJECTS/eDNA_PROJECT/Fungi_Bacteria/Analysis/postmetabaR/")

edenbc<- read.csv("Eden_soil_chem_data2.csv")
edenbc3<-edenbc[1:4]

microclim6<-merge(microclim5,edenbc3,by="SampleSample")

microclim_plot_moist <- 
  
  ggplot(microclim6, aes(x=SampleSample, y=Moisture, color=Ecosystem))+theme_classic()+
  geom_boxplot()

microclim_plot_moist+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

microclim_plot_temp <- 
  
  ggplot(microclim6, aes(x=SampleSample, y=Temperature, color=Ecosystem))+theme_classic()+
  geom_boxplot()

microclim_plot_temp+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

