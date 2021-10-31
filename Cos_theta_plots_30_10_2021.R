# This script will perform Redundancy Analyses, 
# relating phytoplankton community composition and environmental conditions.
# The results will be represented as cos-theta matrices, 
# which represent the angle between samples and environmental vectors-
# a measure for the distribution of environmental affinities across the samples.

setwd('D:/Documents/')
# You will need to set the path to the location
# that you saved the following files
# 'Biomasses_40.csv'
# 'Environment_40.csv'
# 'Absorption_spectra_40.csv'
# 'Cell_counts_40.csv'
# 'Pigments_40.csv'
# If you do not have a package, run the code 'install.packages('package name')'

library('vegan')
# Code for multivariate statistical analyses
library(ggplot2)
# Plotting package
library(reshape2)
# Data tidying
library(cowplot)
# Extracting plot elements
library(ggpubr)
# Combine subplots

colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Colour blind palette for plotting


Biomasses<-read.csv('Biomasses_40.csv')
Environment<-read.csv('Environment_40.csv')
Absorption_spectra<-read.csv('Absorption_spectra_40.csv')
Cell_counts<-read.csv('Cell_counts_40.csv')
Pigments<-read.csv('Pigments_40.csv')
# Load the requisite datasets

Julian<-Biomasses$Julian
# Julian day of the year that biomass samples were taken
diatom_masses<-rowSums(Biomasses[,5:55])/rowSums(Biomasses[,5:110])
dinoflagellate_masses<-(rowSums(Biomasses[,56:98])/(rowSums(Biomasses[,5:110])))
ciliate_masses<-(rowSums(Biomasses[,99:110])/(rowSums(Biomasses[,5:110])))
# Relative contribution of diatoms/different groups to microphytoplankton

Order<-order(Julian) # you can choose to order by diatom_masses instead if you wish
J<-(Julian/365)*12 # Quick conversion of Julian day to months (rough)


Absorption_data<-Absorption_spectra[,5:355]/rowSums(Absorption_spectra[,5:355])
# Normalise to integral
Cell_data<- cbind(log10(Cell_counts[,5:110]+1),log10((Cell_counts[,111:114]*1000)+1))
Cell_data<-Cell_data/rowSums(Cell_data)
# Log10 + 1 transform, then normalise to integral
Pigment_data<-Pigments[,5:32]/rowSums(Pigments[,5:32])
# Normalise to integral

cell_rda<-rda( Cell_data ~ . ,Environment[,-c(11:14)])
optics_rda<-rda( Absorption_data ~ . ,Environment[,-c(11:14)])
pigment_rda<-rda( Pigment_data ~ . ,Environment[,-c(11:14)])
# Perform Redundancy Analyses 


biomasses<-as.data.frame(cbind(1:40,diatom_masses[Order],
dinoflagellate_masses[Order],
ciliate_masses[Order]))
colnames(biomasses)[2]<-'diatoms'
colnames(biomasses)[3]<-'dinoflagellates'
colnames(biomasses)[4]<-'ciliates'
biomasses<-melt(biomasses,id.var='V1')
# Make a data frame and arrange by sample

p<-8 # font size

biomass_bars<-ggplot(data=biomasses,aes(fill=variable,x=value,y=V1))+ # Plot based on biomass
geom_bar(position='stack',stat='identity',width=1)+ # Stacked bar chart
scale_fill_manual("Biomass", values= c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5] ))+
# With a colourblind palette
theme(legend.position='left',
legend.key.width=unit(0.3,'cm'),
legend.key.height=unit(0.3,'cm'),
legend.title=element_text(size=p),
legend.text=element_text(size=p),
axis.title.x=element_blank(),
axis.title.y=element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank(),
axis.text.y=element_blank(),
axis.text.x=element_blank(),
panel.background=element_blank())

# This is a stacked bar plot of the relative dominance of diatoms, dinoflagellates and ciliates with time
# type 'biomass_bars' into the terminal to view

biomass_legend<-cowplot::get_legend(biomass_bars) # Extract legend

biomass_bars<-ggplot(data=biomasses,aes(fill=variable,x=value,y=V1))+
geom_bar(position='stack',stat='identity',width=1)+
scale_fill_manual("Biomass", values= c("diatoms" = colorBlindGrey8[4], "dinoflagellates" = colorBlindGrey8[8], "ciliates" = colorBlindGrey8[5] ))+
	labs(fill = 'Biomass' )+
theme(legend.position='none',
axis.title.x=element_blank(),
axis.title.y=element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank(),
axis.text.y=element_blank(),
axis.text.x=element_blank(),
panel.background=element_blank())

# Legend removed

month_data<-as.data.frame(J[Order])
colnames(month_data)<-'Month'
J_bars<-ggplot(data=month_data,aes(fill=Month,x=1,y=c(1:length(J)) )  )+
geom_tile()+
scale_fill_gradientn(colors=c(low='blue',mid='yellow',high='red'))+ # Adjust if colour deficiency poses a problem
theme(legend.position='left',
legend.key.width=unit(0.3,'cm'),
legend.key.height=unit(0.3,'cm'),
legend.title=element_text(size=p),
legend.text=element_text(size=p),
axis.title.x=element_blank(),
axis.title.y=element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank(),
axis.text.y=element_blank(),
axis.text.x=element_blank(),
panel.background=element_blank())

# This is a colour gradient plot of the month-of-sampling. A red-yellow-blue gradient is chosen 

J_legend<-cowplot::get_legend(J_bars) # Extract legend

J_bars<-ggplot(data=month_data,aes(fill=Month,x=1,y=c(1:length(J)) )  )+
geom_tile()+
scale_fill_gradientn(colors=c(low='blue',mid='yellow',high='red'))+
theme(legend.position='none',
axis.title.x=element_blank(),
axis.title.y=element_blank(),
axis.ticks.x=element_blank(),
axis.ticks.y=element_blank(),
axis.text.y=element_blank(),
axis.text.x=element_blank(),
panel.background=element_blank())

angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}
# This is a function to calculate the angle between samples and environmental vectors across RDA-space

x<-scores(choices=c(1:optics_rda$CCA$rank),optics_rda,scaling=3,display='sites')
# Sample vectors for optics dataset
y<-scores(choices=c(1:optics_rda$CCA$rank),optics_rda,scaling=3,display='bp')
# Environmental vectors for optics dataset

store1<-matrix(NaN,dim(x)[1],dim(y)[1])
for(i in 1: (dim(x)[1]) ){
	for(j in 1: (dim(y)[1]) ){
		store1[i,j]<-cos(angle(t(as.matrix(x[i,])),as.matrix(y[j,])))
	}
}
# Angles between vectors computed

t<-melt(store1[Order,]) 
# Organise by Julian day

image_raw_optics<-ggplot(data=t,aes(x=Var2,y=Var1,fill=value))+ # Plot of cos-theta matrix for absorption spectral data
geom_tile()+ # tile plot
scale_fill_gradient2(low="blue", mid="white", high="red",midpoint=0,limits=c(-1,1))+ # Red-blue colour scale
	labs(fill = expression(bold("cos"[theta])) )+ # Labelling 
scale_x_continuous("",breaks=c(1:10),labels=parse(text =  c(expression(bold("O"[2])) , expression(bold(Salinity)),
 expression(bold(sigma[theta])),expression(bold(Temperature)),expression(bold(Stratification)),expression(bold(Ammonium)),
expression(bold(Silicic~acid)),expression(bold(Nitrate)),
expression(bold(PAR)),
expression(bold(Depth)) ) ) ,position='top')+ # Labelling
theme_bw()+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0),
axis.text.y=element_blank(),
axis.title.y=element_blank(),
legend.position='right',
legend.key.width=unit(0.3,'cm'),
legend.key.height=unit(0.3,'cm'),
legend.title=element_text(size=p),
legend.text=element_text(size=p))

# Image of cos-theta matrix for RDA based on optics dataset
# I will desist with further comments for plots; we will produce subplots for each data type from this point

legend<-cowplot::get_legend(image_raw_optics)
# Extract legend

image_raw_optics<-ggplot(data=t,aes(x=Var2,y=Var1,fill=value))+
geom_tile()+
scale_fill_gradient2(low="blue", mid="white", high="red",midpoint=0,limits=c(-1,1))+
	labs(fill = expression(bold("cos"[theta])) )+
scale_x_continuous("",breaks=c(1:10),labels=parse(text =  c(expression(bold("O"[2])) , expression(bold(Salinity)),
 expression(bold(sigma[theta])),expression(bold(Temperature)),expression(bold(Stratification)),expression(bold(Ammonium)),
expression(bold(Silicic~acid)),expression(bold(Nitrate)),
expression(bold(PAR)),
expression(bold(Depth)) ) ) ,position='top')+
theme_bw()+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0, size=p-2),
axis.text.y=element_blank(),
axis.title.y=element_blank(),
legend.position='none')

x<-scores(choices=c(1:cell_rda$CCA$rank),cell_rda,scaling=3,display='sites')
y<-scores(choices=c(1:cell_rda$CCA$rank),cell_rda,scaling=3,display='bp')

store2<-matrix(NaN,dim(x)[1],dim(y)[1])
for(i in 1: (dim(x)[1]) ){
	for(j in 1: (dim(y)[1]) ){
		store2[i,j]<-cos(angle(t(as.matrix(x[i,])),as.matrix(y[j,])))
	}
}

t<-melt(store2[Order,])

image_cell<-ggplot(data=t,aes(x=Var2,y=Var1,fill=value))+
geom_tile()+
scale_fill_gradient2(low="blue", mid="white", high="red",midpoint=0,limits=c(-1,1))+
	labs(fill = expression(bold("cos"[theta])) )+
scale_x_continuous("",breaks=c(1:10),labels=parse(text =  c(expression(bold("O"[2])) , expression(bold(Salinity)),
 expression(bold(sigma[theta])),expression(bold(Temperature)),expression(bold(Stratification)),expression(bold(Ammonium)),
expression(bold(Silicic~acid)),expression(bold(Nitrate)),
expression(bold(PAR)),
expression(bold(Depth)) ) ) ,position='top')+
theme_bw()+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0, size=p-2),
axis.text.y=element_blank(),
axis.title.y=element_blank(),
legend.position='none')


x<-scores(choices=c(1:pigment_rda$CCA$rank),pigment_rda,scaling=3,display='sites')
y<-scores(choices=c(1:pigment_rda$CCA$rank),pigment_rda,scaling=3,display='bp')

store3<-matrix(NaN,dim(x)[1],dim(y)[1])
for(i in 1: (dim(x)[1]) ){
	for(j in 1: (dim(y)[1]) ){
		store3[i,j]<-cos(angle(t(as.matrix(x[i,])),as.matrix(y[j,])))
	}
}

t<-melt(store3[Order,])

image_pigment<-ggplot(data=t,aes(x=Var2,y=Var1,fill=value))+
geom_tile()+
scale_fill_gradient2(low="blue", mid="white", high="red",midpoint=0,limits=c(-1,1))+
	labs(fill = expression(bold("cos"[theta])) )+
scale_x_continuous("",breaks=c(1:10),labels=parse(text =  c(expression(bold("O"[2])) , expression(bold(Salinity)),
 expression(bold(sigma[theta])),expression(bold(Temperature)),expression(bold(Stratification)),expression(bold(Ammonium)),
expression(bold(Silicic~acid)),expression(bold(Nitrate)),
expression(bold(PAR)),
expression(bold(Depth)) ) ) ,position='top')+
theme_bw()+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0, size=p-2),
axis.text.y=element_blank(),
axis.title.y=element_blank(),
legend.position='none')


# Now we perform the same analysis again, but rotate the environmental data by principal component decomposition

explan<-prcomp(Environment[,-c(11:14)],scale=T) 
# PCA decomposition

cell_rda<-rda( Cell_data ~ . ,as.data.frame(explan$x))
optics_rda<-rda( Absorption_data ~ . ,as.data.frame(explan$x))
pigment_rda<-rda( Pigment_data ~ . ,as.data.frame(explan$x))
# Run analyses


x<-scores(choices=c(1:optics_rda$CCA$rank),optics_rda,scaling=3,display='sites')
y<-scores(choices=c(1:optics_rda$CCA$rank),optics_rda,scaling=3,display='bp')

store4<-matrix(NaN,dim(x)[1],dim(y)[1])
for(i in 1: (dim(x)[1]) ){
	for(j in 1: (dim(y)[1]) ){
		store4[i,j]<-cos(angle(t(as.matrix(x[i,])),as.matrix(y[j,])))
	}
}

t<-melt(store4[Order,])

image_raw_optics_pc<-ggplot(data=t,aes(x=Var2,y=Var1,fill=value))+
geom_tile()+
scale_fill_gradient2(low="blue", mid="white", high="red",midpoint=0,limits=c(-1,1))+
	labs(fill = expression(bold("cos"[theta])) )+
scale_x_continuous("",breaks=c(1:10),labels=c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10') ,position='top')+
theme_bw()+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0 , size=p-2),
axis.text.y=element_blank(),
axis.title.y=element_blank(),
legend.position='none')


x<-scores(choices=c(1:cell_rda$CCA$rank),cell_rda,scaling=3,display='sites')
y<-scores(choices=c(1:cell_rda$CCA$rank),cell_rda,scaling=3,display='bp')

store5<-matrix(NaN,dim(x)[1],dim(y)[1])
for(i in 1: (dim(x)[1]) ){
	for(j in 1: (dim(y)[1]) ){
		store5[i,j]<-cos(angle(t(as.matrix(x[i,])),as.matrix(y[j,])))
	}
}

t<-melt(store5[Order,])

image_cell_pc<-ggplot(data=t,aes(x=Var2,y=Var1,fill=value))+
geom_tile()+
scale_fill_gradient2(low="blue", mid="white", high="red",midpoint=0,limits=c(-1,1))+
	labs(fill = expression(bold("cos"[theta])) )+
scale_x_continuous("",breaks=c(1:10),labels=c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10') ,position='top')+
theme_bw()+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0 , size=p-2),
axis.text.y=element_blank(),
axis.title.y=element_blank(),
legend.position='none')

x<-scores(choices=c(1:pigment_rda$CCA$rank),pigment_rda,scaling=3,display='sites')
y<-scores(choices=c(1:pigment_rda$CCA$rank),pigment_rda,scaling=3,display='bp')

store6<-matrix(NaN,dim(x)[1],dim(y)[1])
for(i in 1: (dim(x)[1]) ){
	for(j in 1: (dim(y)[1]) ){
		store6[i,j]<-cos(angle(t(as.matrix(x[i,])),as.matrix(y[j,])))
	}
}

t<-melt(store6[Order,])


image_pigment_pc<-ggplot(data=t,aes(x=Var2,y=Var1,fill=value))+
geom_tile()+
scale_fill_gradient2(low="blue", mid="white", high="red",midpoint=0,limits=c(-1,1))+
	labs(fill = expression(bold("cos"[theta])) )+
scale_x_continuous("",breaks=c(1:10),labels=c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10') ,position='top')+
theme_bw()+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=0 , size=p-2),
axis.text.y=element_blank(),
axis.title.y=element_blank(),
legend.position='none')


# It is now time to combine the subplots, margin plots and legends into a larger figure. 

data_plot1<-ggarrange(nrow=1,ncol=5,J_bars,biomass_bars,image_cell,image_pigment,image_raw_optics,align='h',widths=c(1/3,2/3,2,2,2),
labels=c('','','(A) Cell counts',' (B) Pigments','(C) Absorption spectra'),font.label = list(size = 8, color = "black"),hjust=-.1)

data_plot2<-ggarrange(nrow=1,ncol=5,J_bars,biomass_bars,image_cell_pc,image_pigment_pc,image_raw_optics_pc,align='h',widths=c(1/3,2/3,2,2,2),
labels=c('','','(D) Cell counts','(E) Pigments','(F) Absorption spectra'),font.label = list(size = 8, color = "black"),hjust=-.1)

combined_data<-ggarrange(nrow=2,data_plot1,data_plot2)

legends<-ggarrange(J_legend,biomass_legend,legend,ncol=1)

ggarrange(combined_data,legends,ncol=2,widths=c(8,1.5))
# Produce final plot



