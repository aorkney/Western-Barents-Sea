# This script will produce pie charts of major phytoplankton biomasses, distributed 
# across maps of sampling expeditions.
# The data will then be represented as a stacked bar chart, organised by seasonality.

setwd('D:/Documents/Final_nutrients_sian')
# You will need to set the path to the location
# that you saved the following file
# 'Biomasses_40.csv'
# If you do not have a package, run the code 'install.packages('package name')'

Biomasses<-read.csv('Biomasses_40.csv')

library(ggplot2)
# This package is required for plotting
library(reshape2)
# This package is required to manage data
library(scatterpie)
# This package is required to plot pie charts
library(ggrepel)
# This package is required to plot labels
library(ggpubr)
# This package is required to combine subplots


colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# This is a colourblind palette, sourced from StackExchange

cruises<-Biomasses$cruises
# This is a list of the cruises that the biomass samples originate from

Julian<-Biomasses$Julian
# Julian day of the year that biomass samples were taken

diatom_masses<-rowSums(Biomasses[,5:55])
# Sum the diatom masses
dinoflagellate_masses<-rowSums(Biomasses[,56:98])
# Sum the dinoflagellate masses
ciliate_masses<-rowSums(Biomasses[,99:110])
# Sum the ciliate masses
PNAN_masses<-rowSums(Biomasses[,111:112])
# Sum the Photo-Nanoflagellate masses
HNAN_masses<-rowSums(Biomasses[,113:114])
# Sum the Hetero-Nanoflagellate masses

biomasses<-as.data.frame(cbind(Biomasses$Lon[order(Julian)],Biomasses$Lat[order(Julian)],1:40,diatom_masses[order(Julian)],
dinoflagellate_masses[order(Julian)],ciliate_masses[order(Julian)],PNAN_masses[order(Julian)],
HNAN_masses[order(Julian)]))
colnames(biomasses)[1]<-'Lon'
colnames(biomasses)[2]<-'Lat'
colnames(biomasses)[3]<-'V1'
colnames(biomasses)[4]<-'diatoms'
colnames(biomasses)[5]<-'dinoflagellates'
colnames(biomasses)[6]<-'ciliates'
colnames(biomasses)[7]<-'PNAN'
colnames(biomasses)[8]<-'HNAN'
# We have arranged the masses into a dataframe

radius<-((rowSums(biomasses[,4:8]))/(10^8.5) )^.5
# Define piechart radii

biomasses<-cbind(biomasses,radius)
# Add this to our dataframe

rownames(biomasses)<-c(1:40)
# Define rownames

biomasses$Lon<-biomasses$Lon+c(1:40)/1000
# Disrupt non-unique geographical locations by a small amount (necessary for plotting functions to work)


world <- map_data('world')
# Load world map data
world<-world[which(world$long>10 & world$long <40 & world$lat>69 & world$lat <85),]
# Restrict our map data to the Barents Sea region

p<-3 # Font size

spring_map<-
ggplot(world, aes(long, lat)) + # Create plot with world map data
geom_map(map=world, aes(map_id=region), color="black", fill='black') + # Plot the world map data
    coord_quickmap()+ 
ylim(69,85)+ # Set latitude limits
xlim(10,35)+ # Set longitude limits
 geom_scatterpie(aes(x=Lon, y=Lat, group=V1, r=radius ), data=biomasses[which( cruises[order(Julian)]=='Spring_2018'),],
cols= colnames(biomasses)[-c(1:3,9)],color=NA )+  # Plot piecharts of biomass
scale_fill_manual(values= c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5],
"PNAN"=colorBlindGrey8[6],"HNAN"=colorBlindGrey8[7]) ) + # Use colourblind palette 
 coord_equal()+ # Set coordinates to equal 
geom_scatterpie_legend(radius=c(0.3976354,0.9740037), x=10, y=70, n=3, labeller= function(x) x=c(50,300) )+ # Make legend for mass
ylab(expression(paste('Latitude',,' ('^o,'N)')) )+ # Label latitude
xlab(expression(paste('Longitude',,' ('^o,'E)')) )+ # Label longitude 
geom_label_repel(data=biomasses[which( cruises[order(Julian)]=='Spring_2018'),], aes(x=Lon,y=Lat), label=order(Julian)[which( cruises[order(Julian)]=='Spring_2018')],alpha=.5,size=p)+
# Unique labels for samples
theme(legend.position='none',
axis.line.x=element_line(),
axis.line.y=element_line()) # Tidy the plotting area

# You can call this plot by typing 'spring_map' in the terminal

# I will now desist with commenting the code for the remaining map subplots


summer_2017_map<-
ggplot(world, aes(long, lat)) +
geom_map(map=world, aes(map_id=region), color="black", fill='black') +
    coord_quickmap()+
ylim(69,85)+
xlim(10,35)+
 geom_scatterpie(aes(x=Lon, y=Lat, group=V1, r=radius ), data=biomasses[which( cruises[order(Julian)]=='Summer_2017'),],
cols= colnames(biomasses)[-c(1:3,9)],color=NA )+ 
scale_fill_manual(values= c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5],
"PNAN"=colorBlindGrey8[6],"HNAN"=colorBlindGrey8[7]) ) + coord_equal()+
geom_label_repel(data=biomasses[which( cruises[order(Julian)]=='Summer_2017'),], aes(x=Lon,y=Lat), label=order(Julian)[which( cruises[order(Julian)]=='Summer_2017')],alpha=.5,size=p)+
ylab(expression(paste('Latitude',,' ('^o,'N)')) )+
xlab(expression(paste('Longitude',,' ('^o,'E)')) )+
theme(legend.position='none',
axis.line.x=element_line(),
axis.line.y=element_line())


summer_2018_map<-
ggplot(world, aes(long, lat)) +
geom_map(map=world, aes(map_id=region), color="black", fill='black') +
    coord_quickmap()+
ylim(69,85)+
xlim(10,35)+
 geom_scatterpie(aes(x=Lon, y=Lat, group=V1, r=radius ), data=biomasses[which( cruises[order(Julian)]=='Summer_2018'),],
cols= colnames(biomasses)[-c(1:3,9)],color=NA )+ 
scale_fill_manual(values= c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5],
"PNAN"=colorBlindGrey8[6],"HNAN"=colorBlindGrey8[7]) ) + coord_equal()+
ylab(expression(paste('Latitude',,' ('^o,'N)')) )+
xlab(expression(paste('Longitude',,' ('^o,'E)')) )+
geom_label_repel(data=biomasses[which( cruises[order(Julian)]=='Summer_2018'),], aes(x=Lon,y=Lat), label=order(Julian)[which( cruises[order(Julian)]=='Summer_2018')],alpha=.5,size=p)+
theme(legend.position='none',
axis.line.x=element_line(),
axis.line.y=element_line())


# We will now produce a stacked barplot, with biomasses organised by Julian day along the x-axis.

biomasses<-as.data.frame(cbind(1:40,diatom_masses[order(Julian)]/10^6,
dinoflagellate_masses[order(Julian)]/10^6,
ciliate_masses[order(Julian)]/10^6,PNAN_masses[order(Julian)]/10^6,
HNAN_masses[order(Julian)]/10^6))
colnames(biomasses)[2]<-'diatoms'
colnames(biomasses)[3]<-'dinoflagellates'
colnames(biomasses)[4]<-'ciliates'
colnames(biomasses)[5]<-'PNAN'
colnames(biomasses)[6]<-'HNAN'

biomasses2<-melt(biomasses,id.var='V1')
# Organise by sample


biomass_bars<-ggplot(data=biomasses2,aes(fill=variable,y=value,x=V1))+ # Produce a plot with the biomass data
geom_bar(position='stack',stat='identity',width=1)+ # Stacked bar plot
scale_fill_manual("Biomass", values = c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5],
"PNAN"=colorBlindGrey8[6],"HNAN"=colorBlindGrey8[7]))+ # Colourblind palette
xlab('Ordinal time (by Julian day)')+ # x-axis title
ylab( expression(paste('Biomass',,' (',mu,'gC L'^-1,')')) )+ # y-axis title
theme(legend.position='right',
axis.ticks.x=element_blank(),
axis.text.x=element_blank(),
axis.line.x=element_line(),
axis.line.y=element_line(),
	legend.key.width= unit (1/3,'cm'),
	legend.key.height= unit (1/3,'cm'))

# You can call this plot by typing 'biomass_bars' in the terminal

ggarrange(spring_map,summer_2018_map,summer_2017_map,
biomass_bars,labels=c('(A)','(B)', '(C)','(D)'),hjust=.02,vjust=1.2)
# Product a combined plot



