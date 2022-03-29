# This script will produce triplots of the leading axes that 
# define the multivariate dependencies of
# phytoplankton community structure on environmental conditions.
# The spaces are illustrated with representative species, 
# taxon-indicative pigments and absorption spectra, respectively. 

setwd('D:/Documents/')
# You will need to set the path to the location
# that you saved the following files
# 'Environment_40.csv'
# 'Absorption_spectra_40.csv'
# 'Cell_counts_40.csv'
# 'Pigments_40.csv'
# If you do not have a package, run the code 'install.packages('package name')'

library('vegan')
# Code for multivariate statistical analyses
library(ggplot2)
# Plotting package
library(cowplot)
# Extracting plot elements
library(ggpubr)
# Combine subplots

Environment<-read.csv('Environment_40.csv')
Absorption_spectra<-read.csv('Absorption_spectra_40.csv')
Cell_counts<-read.csv('Cell_counts_40.csv')
Pigments<-read.csv('Pigments_40.csv')
# Load the requisite datasets

explan<-prcomp(Environment[,-c(11:14)],scale=T)
# Perform PCA decomposition on Environmental data, in order
# to eliminate multiple co-linearity

Absorption_data<-Absorption_spectra[,5:355]/rowSums(Absorption_spectra[,5:355])
# Normalise to integral
Cell_data<- cbind(log10(Cell_counts[,5:110]+1),log10((Cell_counts[,111:114]*1000)+1))
Cell_data<-Cell_data/rowSums(Cell_data)
# Log10 + 1 transform, then normalise to integral
Pigment_data<-Pigments[,5:32]/rowSums(Pigments[,5:32])
# Normalise to integral

cell_rda<-rda( Cell_data ~ . ,as.data.frame(explan$x))
optics_rda<-rda( Absorption_data ~ . ,as.data.frame(explan$x))
pigment_rda<-rda( Pigment_data ~ . ,as.data.frame(explan$x))
# Perform redundancy analysis

hydrog_vector<-rep('Barents',dim(Environment)[1] )
hydrog_vector[which(Environment$temp <0 & Environment$salinity <34.8)]<- 'Arctic' #Arw
hydrog_vector[which(Environment$temp >0 & Environment$temp <3 & Environment$salinity <34.4)]<- 'Melt' #Mw
hydrog_vector[which(Environment$temp >3 & Environment$salinity <34.8)]<- 'Coastal' #Cw
hydrog_vector[which(Environment$temp >3 & Environment$salinity >34.8)]<- 'Atlantic' #AtW
# Define a vector of the Watermass classifications of each sample

d<-as.data.frame(cbind(scores(cell_rda,display='sites'),
hydrog_vector))
d[,1:2]<-as.matrix(sapply(d[,1:2], as.character ))  
d[,1:2]<-as.matrix(sapply(d[,1:2], as.numeric ))  
d[,1]<-d[,1]/max(abs(d[,1]))
d[,2]<-d[,2]/max(abs(d[,2]))
colnames(d)[3]<-'Watermass'
# Make a dataframe for the cell counts analysis

cell_segments<-as.data.frame(cbind(rep(0,dim(Environment)[2]-4),rep(0,dim(Environment)[2]-4),
summary(cell_rda)$biplot[,1],summary(cell_rda)$biplot[,2]))/max(abs(summary(cell_rda)$biplot))
# Triplot vectors for environmental variance axes

cell_segments<-cell_segments/max(abs(cell_segments))
# Normalise

cell_var<-c(1,2,3,5,6)
# The chosen axes for plotting (those identified as significant in the significance testing script)

cell_picks<-c(18,14,42,36,97,40,1,72,80,9)
# Select the species and genera to plot in the triplot
dcells<-as.data.frame(cbind(rep(0,10),rep(0,10),
scores(cell_rda,display='species')[cell_picks,1],scores(cell_rda,display='species')[cell_picks,2])/ (max(abs(scores(cell_rda,display='species')[cell_picks,]))*1.2) )
dcells<-dcells/max(abs(dcells))

cell_names<-c(expression(italic('C.soc')),
expression(italic('T.nor')),expression(italic('N.fri')),
expression(italic('Fra')),expression(italic('Lab')),expression(italic('Att')),
expression(italic('P.del')),expression(italic('Scr')),expression(italic('G.spi')),
expression(italic('R.heb')))
# Abbreviate the names of the species and genera

cell_rda_plot<-
ggplot(data=d,aes(x=RDA1,y=RDA2,col=Watermass))+ # Plot with cell RDA data and watermass properties
ylim(c(-1,1))+
xlim(c(-1,1))+ # Set axis limits
scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.15, 0.15))+ 
geom_segment(data=cbind(cell_segments[cell_var,]),aes(x=V1,y=V2,xend=V3,yend=V4),size=2, 
alpha=1,color='gray50')+ # Plot vectors for environmental axes
geom_segment(data=dcells,aes(x=V1,y=V2,xend=V3,yend=V4),size=2,
alpha=1,color='grey')+  # Plot vectors for representative species
geom_point(aes(shape=Watermass,col=Watermass),size=4,alpha=1)+ # Plot sample points 
scale_colour_manual("Watermass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Watermass", values = c("Arctic" = 15, "Atlantic" = 16, "Barents" = 17, "Coastal" = 18, "Melt" = 19))+ # Colours and shapes
geom_text(label=c(1:40),col='black',size=2)+ # Label points with unique numbers 
geom_label(data=dcells,aes(x=V3,y=V4),label=cell_names,color='black',
size=2,alpha=.5)+ # Add species and genera abbreviations to the plot space 
geom_text(data=cell_segments[cell_var,]*1,aes(x=V3,y=V4),label=c(colnames(explan$x)[cell_var]),color='black',
size=4)+ # Add environmental axis labels 
theme(legend.position='right',
axis.line.x=element_line(),
axis.line.y=element_line(),
panel.background=element_blank())
# Completed plot; type 'cell_rda_plot' into terminal to view. 

legend<-cowplot::get_legend(cell_rda_plot)
# Obtain legend.

cell_rda_plot<-
ggplot(data=d,aes(x=RDA1,y=RDA2,col=Watermass))+
ylim(c(-1,1))+
xlim(c(-1,1))+
scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.15, 0.15))+
geom_segment(data=cbind(cell_segments[cell_var,]),aes(x=V1,y=V2,xend=V3,yend=V4),size=2,
alpha=1,color='gray50')+
geom_segment(data=dcells,aes(x=V1,y=V2,xend=V3,yend=V4),size=2,
alpha=1,color='grey')+
geom_point(aes(shape=Watermass,col=Watermass),size=4,alpha=1)+
scale_colour_manual("Watermass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Watermass", values = c("Arctic" = 15, "Atlantic" = 16, "Barents" = 17, "Coastal" = 18, "Melt" = 19))+
geom_text(label=c(1:40),col='black',size=2)+
geom_label(data=dcells,aes(x=V3,y=V4),label=cell_names,color='black',
size=2,fontface='bold',alpha=.5,seed=1234)+
geom_label(data=dcells,aes(x=V3,y=V4),label=cell_names,color='black',
size=2,fontface='bold',seed=1234,fill=NA)+
geom_label(data=cell_segments[cell_var,]*1,aes(x=V3,y=V4),label=c(colnames(explan$x)[cell_var]),color='black',
size=2,fontface='bold',alpha=.5,)+
theme(legend.position='none',
axis.line.x=element_line(),
axis.line.y=element_line(),
panel.background=element_blank())
# Plot with legend removed

# Now we will repeat this exercise for the pigment dataset

d2<-as.data.frame(cbind(scores(pigment_rda,display='sites'),
hydrog_vector))
d2[,1:2]<-as.matrix(sapply(d2[,1:2], as.character ))  
d2[,1:2]<-as.matrix(sapply(d2[,1:2], as.numeric ))  
d2[,1]<-d2[,1]/max(abs(d2[,1]))
d2[,2]<-d2[,2]/max(abs(d2[,2]))
colnames(d2)[3]<-'Watermass'

pigment_segments<-as.data.frame(cbind(rep(0,dim(Environment)[2]-4),rep(0,dim(Environment)[2]-4),
summary(pigment_rda)$biplot[,1],summary(pigment_rda)$biplot[,2]))/max(abs(summary(pigment_rda)$biplot))

pigment_segments<-pigment_segments/max(abs(pigment_segments))

pigment_var<-c(1,2,3,5,6)

pigment_picks<-c(25,8,13,24,1,7,6)
pigment_names<-c('chl-a','fuco','19-hex','chl-b','chl-c3','19-but','perid')

dpigments<-as.data.frame(cbind(rep(0,length(pigment_picks)),rep(0,length(pigment_picks)),
scores(pigment_rda,display='species')[pigment_picks,1],scores(pigment_rda,display='species')[pigment_picks,2])/ (max(abs(scores(pigment_rda,display='species')[pigment_picks,]))*1.2) )

dpigments<-dpigments/max(abs(dpigments))

pigment_rda_plot<-ggplot(data=d2,aes(x=RDA1,y=RDA2,col=Watermass))+
scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.15, 0.15))+
geom_segment(data=pigment_segments[pigment_var,],aes(x=V1,y=V2,xend=V3,yend=V4),size=2,
alpha=1,color='gray50')+
geom_segment(data=dpigments,aes(x=V1,y=V2,xend=V3,yend=V4),size=2,
alpha=1,color='grey')+
geom_point(aes(shape=Watermass,col=Watermass),size=4,alpha=1)+
scale_colour_manual("Watermass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Watermass", values = c("Arctic" = 15, "Atlantic" = 16, "Barents" = 17, "Coastal" = 18, "Melt" = 19))+
geom_text(label=c(1:40),col='black',size=2)+
geom_label(data=dpigments,aes(x=V3,y=V4),label=pigment_names,color='black',
size=2,fontface='bold',alpha=.5)+
geom_label(data=pigment_segments[pigment_var,],aes(x=V3,y=V4),label=c(rownames(pigment_segments)[pigment_var]),color='black',
size=2,fontface='bold',alpha=.5,)+
theme(legend.position='none',
axis.line.x=element_line(),
axis.line.y=element_line(),
panel.background=element_blank())

# Now we perform the same as above, but for absorption spectra

d3<-as.data.frame(cbind(scores(optics_rda,display='sites'),
hydrog_vector))
d3[,1:2]<-as.matrix(sapply(d3[,1:2], as.character ))  
d3[,1:2]<-as.matrix(sapply(d3[,1:2], as.numeric ))  
d3[,1]<-d3[,1]/max(abs(d3[,1]))
d3[,2]<-d3[,2]/max(abs(d3[,2]))
colnames(d3)[3]<-'Watermass'

optics_segments<-as.data.frame(cbind(rep(0,dim(Environment)[2]-4),rep(0,dim(Environment)[2]-4),
summary(optics_rda)$biplot[,1],summary(optics_rda)$biplot[,2]))/max(abs(summary(optics_rda)$biplot))

optics_segments<-optics_segments/max(abs(optics_segments))

optics_var<-c(1,2,3,5,6)

# We want to obtain some representative absorption spectra to illustrate the plot

mean_spectrum<-colSums(Absorption_data)/(dim(Absorption_data)[1])
# the mean spectrum

x<-d3$RDA1
y<-d3$RDA2
dex<-as.data.frame(cbind(d3,atan2(y,x)*(180/pi) +180 ))
colnames(dex)[4]<-'colly'
# Define angle from origin

NE<-which(dex$colly>180 & dex$colly <300)
E<-which(dex$colly>150 & dex$colly <180)
S<-which(dex$colly >40 & dex$colly <150)
W<-which(dex$colly <40 | dex$colly >300)
# Define different angular regions of the plot-space

d4<- as.data.frame( cbind(c(400:750),
mean_spectrum/mean_spectrum[277],
colSums(Absorption_data[NE,])/colSums(Absorption_data[NE,])[277],
colSums(Absorption_data[E,])/colSums(Absorption_data[E,])[277],
colSums(Absorption_data[S,])/colSums(Absorption_data[S,])[277],
colSums(Absorption_data[W,])/colSums(Absorption_data[W,])[277]
))
# Define a dataset of representative absorption spectra


optics_rda_plot<-
ggplot()+
scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.15, 0.15))+
geom_segment(data=optics_segments[optics_var,],aes(x=V1,y=V2,xend=V3,yend=V4),size=2,
alpha=1,color='gray50')+
labs(x='RDA1',y='RDA2')+
geom_point(data=d3,aes(x=RDA1,y=RDA2,shape=Watermass,col=Watermass),size=4,alpha=1)+
scale_colour_manual("Watermass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Watermass", values = c("Arctic" = 15, "Atlantic" = 16, "Barents" = 17, "Coastal" = 18, "Melt" = 19))+
geom_text(data=d3, aes(x=RDA1,y=RDA2), label=c(1:40),col='black',size=2)+
geom_label(data=optics_segments[optics_var,],aes(x=V3,y=V4),label=c(rownames(optics_segments)[optics_var]),color='black',
size=2,fontface='bold',alpha=.5,)+
geom_line(data=d4,aes(x=(V1/700)+0.01,y=(V3/6)+0.25))+
geom_line(data=d4,aes(x=(V1/700)+0.01,y=(V2/6)+0.25),linetype='dashed')+ # This is the mean spectrum 
geom_line(data=d4,aes(x=(V1/700)-0.0,y=(V4/6)-0.65))+
geom_line(data=d4,aes(x=(V1/700)-0.0,y=(V2/6)-0.65),linetype='dashed')+
geom_line(data=d4,aes(x=(V1/700)-1.0,y=(V5/6)-.65))+
geom_line(data=d4,aes(x=(V1/700)-1.0,y=(V2/6)-.65),linetype='dashed')+
geom_line(data=d4,aes(x=(V1/700)-1.2,y=(V6/6)+.4))+
geom_line(data=d4,aes(x=(V1/700)-1.2,y=(V2/6)+.4),linetype='dashed')+
# These above lines are plotting the representative absorption spectra between 400 and 750 nm
theme(legend.position='none',
axis.line.x=element_line(),
axis.line.y=element_line(),
panel.background=element_blank())

ggarrange(ncol=2,nrow=2,cell_rda_plot,pigment_rda_plot,optics_rda_plot,
legend,labels=c('(A)','(B)','(C)'),font.label=c(size=12),hjust=-1)

# Arrange the final plot
# We can now visualise the way that environmental variation structures explained variance in the cell counts,
# pigment and absorption spectral datasets. 



