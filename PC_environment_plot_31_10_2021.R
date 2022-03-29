# This script produces biplots of the leading important 
# Principal Components of the dataset of western Barents Sea
# Environmental conditions associated with the phytoplankton samples.

setwd('D:/Documents/')
# You will need to set the path to the location
# that you saved the following file

# 'Environment_40.csv'


# If you do not have a package, run the code 'install.packages('package name')'
library(ggplot2)
# Plotting package
library(ggpubr)
# Combine subplots
library(cowplot)
# Extracting plot elements


Environment<-read.csv('Environment_40.csv')
# Load the data





Julian<-Environment$Julian
# The Julian day of sample collection
J<-(Julian/365)*12 
# Quick conversion of Julian day to months (rough)

explan<-prcomp(Environment[,-c(11:14)],scale=T)
# PCA decomposition 

hydrog_vector<-rep('Barents',dim(Environment)[1])
hydrog_vector[which(Environment$temp <0 & Environment$salinity <34.8)]<- 'Arctic' #Arw
hydrog_vector[which(Environment$temp >0 & Environment$temp <3 & Environment$salinity <34.4)]<- 'Melt' #Mw
hydrog_vector[which(Environment$temp >3 & Environment$salinity <34.8)]<- 'Coastal' #Cw
hydrog_vector[which(Environment$temp >3 & Environment$salinity >34.8)]<- 'Atlantic' #AtW
# Define Watermasses

d<-data.frame(cbind(explan$x,hydrog_vector))
colnames(d)[11]<-'Watermass'
# Make dataframe

d<-cbind(d,J)
colnames(d)[12]<-'Month'
d[,1:10]<-as.matrix(sapply(d[,1:10], as.character ))  
d[,1:10]<-as.matrix(sapply(d[,1:10], as.numeric ))  
# Combine month-of-sampling with dataframe


PC1_2<-ggplot(data=d,aes(x=PC1,y=PC2,col=Watermass))+ # Produce a plot with the dataframe; PC1-2 space
scale_x_continuous(expand = c(0.15, 0.15)) + scale_y_continuous(expand = c(0.15, 0.15))+ # Margin limits
geom_point(aes(shape=Watermass,col=Watermass),size=4)+ # Define colours as watermass
scale_colour_manual("Watermass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Watermass", values = c("Arctic" = 15, "Atlantic" = 16, "Barents" = 17, "Coastal" = 18, "Melt" = 19))+
geom_text(label=c(1:40),col='black',size=2)+
theme(legend.position='right',
axis.line.x=element_line(),
axis.line.y=element_line(),
panel.background=element_blank())

# Type 'PC1_2' directly in the terminal to view this plot

legend<-cowplot::get_legend(PC1_2)
# Extract the plot legend for later use

PC1_2<-ggplot(data=d,aes(x=PC1,y=PC2,col=Watermass))+
scale_x_continuous(expand = c(0.15, 0.15)) + scale_y_continuous(expand = c(0.15, 0.15))+
geom_point(aes(shape=Watermass,col=Watermass),size=4)+
scale_colour_manual("Watermass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Watermass", values = c("Arctic" = 15, "Atlantic" = 16, "Barents" = 17, "Coastal" = 18, "Melt" = 19))+
geom_text(label=c(1:40),col='black',size=2)+
theme(legend.position='none',
axis.line.x=element_line(),
axis.line.y=element_line(),
panel.background=element_blank())
# Remove plot legend from figure

PC1_3<-ggplot(data=d,aes(x=PC1,y=PC3,col=Watermass))+
scale_x_continuous(expand = c(0.15, 0.15)) + scale_y_continuous(expand = c(0.15, 0.15))+
geom_point(aes(shape=Watermass,col=Watermass),size=4)+
scale_colour_manual("Watermass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Watermass", values = c("Arctic" = 15, "Atlantic" = 16, "Barents" = 17, "Coastal" = 18, "Melt" = 19))+
geom_text(label=c(1:40),col='black',size=2)+
theme(legend.position='none',
axis.line.x=element_line(),
axis.line.y=element_line(),
panel.background=element_blank())

# PC1-3 space plot

# We will now produce a plot showing that PC5-6 are predictive of the sampling seasonality
d2<-as.data.frame(cbind(J,predict(lm(J~explan$x[,6]+explan$x[,5])),hydrog_vector))
# Linear additive model relates PC5,6 and month-of-sampling
colnames(d2)<-c('Month','Prediction','Watermass')
d2$Month<-as.numeric(as.character(d2$Month))
d2$Prediction<-as.numeric(as.character(d2$Prediction))

Predict_plot<-ggplot(data=d2,aes(x=Month,y=Prediction,col=Watermass))+
scale_x_continuous(expand = c(0.15, 0.15)) + scale_y_continuous(expand = c(0.15, 0.15))+
geom_point(aes(shape=Watermass,col=Watermass),size=4)+
scale_colour_manual("Watermass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Watermass", values = c("Arctic" = 15, "Atlantic" = 16, "Barents" = 17, "Coastal" = 18, "Melt" = 19))+
geom_text(label=c(1:40),col='black',size=2)+
geom_abline(intercept=0,slope=1,size=2,linetype='dashed')+
theme(legend.position='none',
axis.line.x=element_line(),
axis.line.y=element_line(),
panel.background=element_blank())

# Produce a plot of the prediction verses the actual time of sampling


ggarrange(ncol=2,nrow=2,PC1_2,PC1_3,Predict_plot,
legend,labels=c('(A)','(B)','(C)'),font.label=c(size=12),hjust=0.01)

# Combine the plots into a single image


