# This is a script that will produce a bar-plot
# representing biomasss contributions of different phytoplankton
# assemblages, organised by cruise, latitude and time-of-sampling

setwd('D:/Documents/')
# You will need to set the path to the location
# that you saved the following files
# 'Biomasses_40.csv'
# If you do not have a package, run the code 'install.packages('package name')'

Biomasses<-read.csv('Biomasses_40.csv')
# Read the data

cruises<-Biomasses$cruises
# define a vector of cruises

J<-Biomasses$Julian
J<-(J/365)*12
J<-order(J)
Julian<-Biomasses$Julian[J]
# define a vector of time-of-sampling 

diatom_masses<-(rowSums(Biomasses[,5:55]))
dinoflagellate_masses<-(rowSums(Biomasses[,56:98]))
ciliate_masses<-(rowSums(Biomasses[,99:110]))
PNAN_masses<-(rowSums(Biomasses[,c(111,112)]))
HNAN_masses<-(rowSums(Biomasses[,c(113,114)]))
# Key groups
# HNAN: heterotrophic nanoflagellates
# PNAN: phototrophic nanoflagellates

biomasses<-as.data.frame(cbind(1:40,diatom_masses[J]/10^6,
dinoflagellate_masses[J]/10^6,
ciliate_masses[J]/10^6,PNAN_masses[J]/10^6,
HNAN_masses[J]/10^6,Julian))
colnames(biomasses)[2]<-'diatoms'
colnames(biomasses)[3]<-'dinoflagellates'
colnames(biomasses)[4]<-'ciliates'
colnames(biomasses)[5]<-'PNAN'
colnames(biomasses)[6]<-'HNAN'
colnames(biomasses)[7]<-'Julian'
# Define a data frame

for(i in 1:length(unique(Biomasses$Julian)) ) {
	rows<-which(biomasses$Julian==unique(biomasses$Julian)[i])
	if(i==1){
		bmas<-t(as.matrix(colSums(biomasses[rows,2:6])/length(rows)))
		jfix<-t(as.matrix(unique(biomasses$Julian)[i]))
		pooled<-cbind(bmas,jfix)
	} else  {
		bmas<-t(as.matrix(colSums(biomasses[rows,2:6])/length(rows)))
		jfix<-t(as.matrix(unique(biomasses$Julian)[i]))
		new_row<-cbind(bmas,jfix)
		pooled<-rbind(pooled,new_row)
	}

}
pooled<-as.data.frame(pooled)
colnames(pooled)[6]<-'Julian'
# Pool the data by unique time-of-sampling

library('reshape2')
# Package for organising data frames

biomasses2<-melt(pooled,id.var='Julian')
# Organise data frame

library(ggplot2)
# Package for producing plots

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Colour scheme (colourblind friendly)

#Produce a bar-plot organised by time-of-sampling
biomass_bars<-
ggplot(data=biomasses2,aes(fill=variable,y=value,x=Julian))+
geom_bar(position='stack',stat='identity',width=1)+
scale_fill_manual("Biomass", values = c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5],
"PNAN"=colorBlindGrey8[6],"HNAN"=colorBlindGrey8[7]))+
xlab('Ordinal time (by Julian day)')+
ylim(c(0,350))+
ylab( expression(paste('Biomass',,' (',mu,'gC L'^-1,')')) )+
theme(legend.position='right',
axis.text.x=element_text(angle=90),
axis.line.x=element_line(),
axis.line.y=element_line(),
	legend.key.width= unit (1/3,'cm'),
	legend.key.height= unit (1/3,'cm'))
# Produce a bar-plot of phytoplankton assemblage contributions across the pooled dataset

library(cowplot)
# Package for manipulating plots

legend<-cowplot::get_legend(biomass_bars) #extract legend

block_1<-biomasses2[which(biomasses2$Julian <125),]
block_2<-biomasses2[which(biomasses2$Julian <187 & biomasses2$Julian >150),]
block_3<-biomasses2[which(biomasses2$Julian >187),]
# Define durations of individual cruises for linear interpolations

# Bar-plot organised by time with linear interpolations within cruises
biomass_bars<-
ggplot(data=biomasses2,aes(fill=variable,y=value,x=Julian))+
geom_area(data=block_1, position = 'stack', alpha=0.4) +
geom_area(data=block_2, position = 'stack', alpha=0.4) +
geom_area(data=block_3, position = 'stack', alpha=0.4) +
geom_bar(position='stack',stat='identity',width=1)+
scale_fill_manual("Biomass", values = c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5],
"PNAN"=colorBlindGrey8[6],"HNAN"=colorBlindGrey8[7]))+
xlab('Julian day')+
ylim(c(0,350))+
ylab( expression(paste('Biomass',,' (',mu,'gC L'^-1,')')) )+
theme(legend.position='none',
axis.text.x=element_text(angle=90),
axis.line.x=element_line(),
axis.line.y=element_line(),
	legend.key.width= unit (1/3,'cm'),
	legend.key.height= unit (1/3,'cm'))

select<-which(cruises=='Spring_2018')
lat<-Biomasses$Lat[select]
lat_Spring_2018<-lat
lat<-order(lat)
lat_Spring_2018<-lat_Spring_2018[lat]
# Extract latitudes of samples for spring 2018 cruise

biomasses<-as.data.frame(cbind(1:length(lat),diatom_masses[select][lat]/10^6,
dinoflagellate_masses[select][lat]/10^6,
ciliate_masses[select][lat]/10^6,PNAN_masses[select][lat]/10^6,
HNAN_masses[select][lat]/10^6))
colnames(biomasses)[2]<-'diatoms'
colnames(biomasses)[3]<-'dinoflagellates'
colnames(biomasses)[4]<-'ciliates'
colnames(biomasses)[5]<-'PNAN'
colnames(biomasses)[6]<-'HNAN'
biomasses2<-melt(biomasses,id.var='V1')
# Produce a data frame organised by individual cruises and ordinal latitude

# Plot the spring 2018 samples by ordinal latitude
spring_2018_biomass<-
ggplot(data=biomasses2,aes(fill=variable,y=value,x=V1))+
geom_bar(position='stack',stat='identity',width=1)+
scale_fill_manual("Biomass", values = c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5],
"PNAN"=colorBlindGrey8[6],"HNAN"=colorBlindGrey8[7]))+
xlab( expression(paste('Ordinal latitude',,' ('^o,'N)')) )+
ylim(c(0,350))+
scale_x_continuous(breaks=c(1:length(lat))-0.25, labels=format(round(lat_Spring_2018,digits=1),nsmall=1) )+
ylab( expression(paste('Biomass',,' (',mu,'gC L'^-1,')')) )+
theme(legend.position='none',
axis.ticks.x=element_blank(),
axis.text.x=element_text(angle=90),
axis.line.x=element_line(),
axis.line.y=element_line(),
	legend.key.width= unit (1/3,'cm'),
	legend.key.height= unit (1/3,'cm'))

select<-which(cruises=='Summer_2018')
lat<-Biomasses$Lat[select]
lat_Summer_2018<-lat
lat<-order(lat)
lat_Summer_2018<-lat_Summer_2018[lat]
# Extract latitudes of samples for summer 2018 cruise

biomasses<-as.data.frame(cbind(1:length(lat),diatom_masses[select][lat]/10^6,
dinoflagellate_masses[select][lat]/10^6,
ciliate_masses[select][lat]/10^6,PNAN_masses[select][lat]/10^6,
HNAN_masses[select][lat]/10^6))
colnames(biomasses)[2]<-'diatoms'
colnames(biomasses)[3]<-'dinoflagellates'
colnames(biomasses)[4]<-'ciliates'
colnames(biomasses)[5]<-'PNAN'
colnames(biomasses)[6]<-'HNAN'
biomasses2<-melt(biomasses,id.var='V1')
# Produce a data frame organised by individual cruises and ordinal latitude

# Plot the summer 2018 samples by ordinal latitude
summer_2018_biomass<-
ggplot(data=biomasses2,aes(fill=variable,y=value,x=V1))+
geom_bar(position='stack',stat='identity',width=1)+
scale_fill_manual("Biomass", values = c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5],
"PNAN"=colorBlindGrey8[6],"HNAN"=colorBlindGrey8[7]))+
xlab( expression(paste('Ordinal latitude',,' ('^o,'N)')) )+
ylim(c(0,350))+
scale_x_continuous(breaks=c(1:length(lat))-0.25, labels=format(round(lat_Summer_2018,digits=1),nsmall=1) )+
ylab( expression(paste('Biomass',,' (',mu,'gC L'^-1,')')) )+
theme(legend.position='none',
axis.ticks.x=element_blank(),
axis.text.x=element_text(angle=90),
axis.line.x=element_line(),
axis.line.y=element_line(),
	legend.key.width= unit (1/3,'cm'),
	legend.key.height= unit (1/3,'cm'))


select<-which(cruises=='Summer_2017')
lat<-Biomasses$Lat[select]
lat_Summer_2017<-lat
lat<-order(lat)
lat_Summer_2017<-lat_Summer_2017[lat]
# Extract latitudes of samples for summer 2017 cruise

biomasses<-as.data.frame(cbind(1:length(lat),diatom_masses[select][lat]/10^6,
dinoflagellate_masses[select][lat]/10^6,
ciliate_masses[select][lat]/10^6,PNAN_masses[select][lat]/10^6,
HNAN_masses[select][lat]/10^6))
colnames(biomasses)[2]<-'diatoms'
colnames(biomasses)[3]<-'dinoflagellates'
colnames(biomasses)[4]<-'ciliates'
colnames(biomasses)[5]<-'PNAN'
colnames(biomasses)[6]<-'HNAN'
biomasses2<-melt(biomasses,id.var='V1')
# Produce a data frame organised by individual cruises and ordinal latitude

# Plot the summer_2017 samples by ordinal latitude
summer_2017_biomass<-
ggplot(data=biomasses2,aes(fill=variable,y=value,x=V1))+
geom_bar(position='stack',stat='identity',width=1)+
scale_fill_manual("Biomass", values = c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5],
"PNAN"=colorBlindGrey8[6],"HNAN"=colorBlindGrey8[7]))+
xlab( expression(paste('Ordinal latitude',,' ('^o,'N)')) )+
ylim(c(0,350))+
scale_x_continuous(breaks=c(1:length(lat))-0.25, labels=format(round(lat_Summer_2017,digits=1),nsmall=1) )+
ylab( expression(paste('Biomass',,' (',mu,'gC L'^-1,')')) )+
theme(legend.position='none',
axis.ticks.x=element_blank(),
axis.text.x=element_text(angle=90),
axis.line.x=element_line(),
axis.line.y=element_line(),
	legend.key.width= unit (1/3,'cm'),
	legend.key.height= unit (1/3,'cm'))

library(ggpubr)
# Package for combing plots

plot_1<-ggarrange(spring_2018_biomass,summer_2018_biomass,summer_2017_biomass,
biomass_bars,labels=c('(A)','(B)', '(C)','(D)'),hjust=.02,vjust=1.2)
ggarrange(plot_1,legend,nrow=1,widths=c(1,0.2))
# combine the plots

plot_1
# show the plot with all the sub-plots included 



