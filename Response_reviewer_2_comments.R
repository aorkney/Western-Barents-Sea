# This script conducts statistical analyses in response to the 
# first round of reviewer comments from Reviewer 2
# for the manuscript ID: 860773

setwd('D:/Documents/Final_nutrients_sian')
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

Biomasses<-read.csv('Biomasses_40.csv')
Environment<-read.csv('Environment_40.csv')
Absorption_spectra<-read.csv('Absorption_spectra_40.csv')
Cell_counts<-read.csv('Cell_counts_40.csv')
Pigments<-read.csv('Pigments_40.csv')
# Load the requisite datasets

Absorption_data<-Absorption_spectra[,5:355]/rowSums(Absorption_spectra[,5:355])
# Normalise to integral
Cell_data<- cbind(log10(Cell_counts[,5:110]+1),log10((Cell_counts[,111:114]*1000)+1))
Cell_data<-Cell_data/rowSums(Cell_data)
# Log10 + 1 transform, then normalise to integral
Pigment_data<-Pigments[,5:32]/rowSums(Pigments[,5:32])
# Normalise to integral

explan<-prcomp(Environment[,-c(11:14)],scale=T) 
# PCA decomposition

# Point 1
# What fraction of total variance do Redundancy Analyses explain?

cell_rda<-rda( Cell_data ~ . ,as.data.frame(explan$x))
optics_rda<-rda( Absorption_data ~ . ,as.data.frame(explan$x))
pigment_rda<-rda( Pigment_data ~ . ,as.data.frame(explan$x))
# Run analyses

cell_rda ; optics_rda ; pigment_rda 
# The total constrained variance is reported in each summmary

# Point 2
# What proportion of environmental variance is explained in each PC-axes
# in the PCA decomposition?

summary(explan) 
# Axial variances are reported in this summary

explan$rotation
# The PC-loadings can be viewed to reveal the meaning of individual axes

# Point 3
# Do the results for the analyses derived from pigment ratios change 
# if the pigments are normalised by their variances?

Scaled_Pigment_data<-Pigment_data
for(i in 1:(dim(Pigment_data)[2]) ){
Scaled_Pigment_data[,i]<-(Pigment_data[,i]/sd(Pigment_data[,i]))-mean(Pigment_data[,i])
}
# Perform variance normalisation

Scaled_Pigment_data<-Scaled_Pigment_data[,-5]
# Remove pigments that are invariant ahead of analysis


pigment_rda_scaled<-rda( Scaled_Pigment_data[,-5] ~ . ,as.data.frame(explan$x))
# Perform analysis

pigment_rda_scaled 
# Inspect summary of results

# Let us produce a plot of the RDA space

hydrog_vector<-rep('Barents',dim(Environment)[1] )
hydrog_vector[which(Environment$temp <0 & Environment$salinity <34.8)]<- 'Arctic' #Arw
hydrog_vector[which(Environment$temp >0 & Environment$temp <3 & Environment$salinity <34.4)]<- 'Melt' #Mw
hydrog_vector[which(Environment$temp >3 & Environment$salinity <34.8)]<- 'Coastal' #Cw
hydrog_vector[which(Environment$temp >3 & Environment$salinity >34.8)]<- 'Atlantic' #AtW
# Define a vector of the Watermass classifications of each sample

d2<-as.data.frame(cbind(scores(pigment_rda_scaled,display='sites'),
hydrog_vector))
d2[,1:2]<-as.matrix(sapply(d2[,1:2], as.character ))  
d2[,1:2]<-as.matrix(sapply(d2[,1:2], as.numeric ))  
d2[,1]<-d2[,1]/max(abs(d2[,1]))
d2[,2]<-d2[,2]/max(abs(d2[,2]))
colnames(d2)[3]<-'Watermass'

pigment_segments<-as.data.frame(cbind(rep(0,dim(Environment)[2]-4),rep(0,dim(Environment)[2]-4),
summary(pigment_rda_scaled)$biplot[,1],summary(pigment_rda_scaled)$biplot[,2]))/max(abs(summary(pigment_rda_scaled)$biplot))

pigment_segments<-pigment_segments/max(abs(pigment_segments))

pigment_var<-c(1,2,3,5,6)

#pigment_picks<-c(25,8,13,24,1,7,6)
pigment_picks<-c(24,7,11,23,1,6,5)
pigment_names<-c('chl-a','fuco','19-hex','chl-b','chl-c3','19-but','perid')

dpigments<-as.data.frame(cbind(rep(0,length(pigment_picks)),rep(0,length(pigment_picks)),
scores(pigment_rda_scaled,display='species')[pigment_picks,1],scores(pigment_rda_scaled,display='species')[pigment_picks,2])/ (max(abs(scores(pigment_rda_scaled,display='species')[pigment_picks,]))*1.2) )

dpigments<-dpigments/max(abs(dpigments))

library(ggrepel)
ggplot(data=d2,aes(x=RDA1,y=RDA2,col=Watermass))+
scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.15, 0.15))+
geom_segment(data=pigment_segments[pigment_var,]*2,aes(x=V1,y=V2,xend=V3,yend=V4),size=2,
alpha=1,color='gray50')+
geom_segment(data=dpigments,aes(x=V1,y=V2,xend=V3,yend=V4),size=2,
alpha=1,color='grey')+
geom_point(aes(shape=Watermass,col=Watermass,fill=Watermass),size=5,alpha=1)+
#scale_colour_manual("Watermass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
#scale_shape_manual("Watermass", values = c("Arctic" = 15, "Atlantic" = 16, "Barents" = 17, "Coastal" = 18, "Melt" = 19))+
scale_fill_manual("Water mass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_colour_manual("Water mass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Water mass", values = c("Arctic" = 22, "Atlantic" = 21, "Barents" = 24, "Coastal" = 23, "Melt" = 25))+
geom_text(label=c(1:40),col='black',size=3)+
geom_label_repel(seed=1234,data=dpigments,aes(x=V3,y=V4),label=pigment_names,color='black',
size=3,fontface='bold',alpha=.5)+
geom_label_repel(seed=1234,data=dpigments,aes(x=V3,y=V4),label=pigment_names,color='black',
size=3,fontface='bold',alpha=1,fill=NA)+
geom_label(data=pigment_segments[pigment_var,]*2,aes(x=V3,y=V4),label=c(rownames(pigment_segments)[pigment_var]),color='black',
size=3,fontface='bold',alpha=.5,)+
theme(legend.position='right',
axis.line.x=element_line(),
axis.text.x=element_text(size=12),
axis.title.x=element_text(size=16),
axis.text.y=element_text(size=12),
axis.title.y=element_text(size=16),
legend.title=element_text(size=16),
legend.text=element_text(size=12),
axis.line.y=element_line(),
panel.background=element_blank())

Test_axes<- colnames(explan$rotation)
Env<-explan$x

for(i in 1: (dim(Environment[,-c(11:14)])[2]) ){
	temporary_rda<-anova(rda( Scaled_Pigment_data ~Env[,i],Z=Env[,-i]))
		print( c( 'proportion of variance explained', temporary_rda$V[1]/sum(temporary_rda$V), Test_axes[i] ) )
		print( c( 'p-value', temporary_rda$"Pr(>F)" [1], Test_axes[i] ) )
}

# Pure effects indicate PC1, PC2, PC5, PC6, PC10 are significant

Pigment_term_effects<-anova(pigment_rda_scaled, by="term")
# Term-wise effects indicate PC1, PC2, PC5, PC6, PC10 are significant

Pigment_margin_effects<-anova(pigment_rda_scaled, by="margin")
# Marginal effects indicate PC1, PC2, PC5, PC6, PC10 are significant




# Point 4
# Do Redundancy Analyses predominantly explain variance in phytoplankton community composition
# that could be attributed to variations in overall biomass or in photophysiological adaptation? 

# We will first test whether pigments are a faithful representation of variation in community structure
# (as opposed to forms of photophysiological adaptation orthogonal to community structure) 
# by performing a redundancy analysis to explore whether cell count structure is explained by pigment ratios

cells_pigments<-rda(Cell_data ~., Pigment_data)
# Some 81% of variance in the cell count data is constrained by variation in pigment ratios

anova(cells_pigments,by='term')
# A term-wise effect statistic indicates that a large number of pigment species underly this explained
# structure

PPC<- Pigment_data$a_carot + Pigment_data$b_carot + Pigment_data$zeaxanthin +
Pigment_data$alloxanthin + Pigment_data$diadinoxanthin

summary(lm(PPC~Environment$PAR))
# A significant p-value is reported (0.005) and a relatively low adjusted R^ (0.17)
# Indicating variation in photosynthetically active radiation significantly explains
# at least some variance in photo protective pigments (PPC)

Pigment_data_without_PPC<- Pigments[,c(5:17,19:21,23,25:30)]/rowSums(Pigments[,c(5:17,19:21,23,25:30)])

pigment_rda_without_PPC<-rda( Pigment_data_without_PPC ~ . ,as.data.frame(explan$x))
# about 39% of variance in these pigments is constrained by the environmental data
# this is broadly similar to the 42-46% constrained by the total suite of pigments

anova(pigment_rda_without_PPC,by='term')
# Fewer axes of environmental variance are identifies as explaining variance in the pigment dataset
# when PPC is removed

cells_pigment_rda_without_PPC<-rda( Cell_data ~ . ,Pigment_data_without_PPC)
# when photoprotective pigments are removed only 69% of variance in the cell count data is explained

# What about when PAR is removed from environmental variability?

explan_without_PAR<-prcomp(Environment[,-c(9,11:14)],scale=T) 
# PCA decomposition without PAR

pigment_rda_without_PAR<-rda( Pigment_data ~ . ,as.data.frame(explan_without_PAR$x))
# 41% of variance is constrained


# Let us be thorough and also undertake variance partitioning analyses


PPC<- Pigment_data$a_carot + Pigment_data$b_carot + Pigment_data$zeaxanthin +
Pigment_data$alloxanthin + Pigment_data$diadinoxanthin

PPC_data<-Pigment_data[,c(27,28,20,18,14)]

mod <- varpart (Cell_data, PPC_data, Pigment_data_without_PPC[-which(sapply(Pigment_data_without_PPC,sd)==0)])
mod <- varpart (Cell_data, PPC_data, Pigment_data_without_PPC[-c(5,22)]) # A test eliminating co-linear pigments

mod <- varpart (Absorption_data, PPC_data, Pigment_data_without_PPC[-which(sapply(Pigment_data_without_PPC,sd)==0)])
mod <- varpart (Absorption_data, PPC_data, Pigment_data_without_PPC[-c(5,22)]) # A test eliminating co-linear pigments

# Let us explore pigment packaging responds to variation in PAR, and what the cause of that might be

log_chla<-log10(Pigments$chl_a)
log_spec_abs<-log10(Absorption_data$X676/Pigments$chl_a)
# These two features should have a negative linear relationship 
# if pigment packaging evolves with (chl-a)

mod<-lm(log_spec_abs~log_chla)
# p value <2e-16, and adjusted R^2 of 0.98
# Very strong negative relationship


hydrog_vector<-rep('Barents',dim(Environment)[1] )
hydrog_vector[which(Environment$temp <0 & Environment$salinity <34.8)]<- 'Arctic' #Arw
hydrog_vector[which(Environment$temp >0 & Environment$temp <3 & Environment$salinity <34.4)]<- 'Melt' #Mw
hydrog_vector[which(Environment$temp >3 & Environment$salinity <34.8)]<- 'Coastal' #Cw
hydrog_vector[which(Environment$temp >3 & Environment$salinity >34.8)]<- 'Atlantic' #AtW
# Define a vector of the Watermass classifications of each sample

d<-as.data.frame(cbind(log_chla,log_spec_abs,hydrog_vector))
d[,1:2]<-as.matrix(sapply(d[,1:2], as.character ))  
d[,1:2]<-as.matrix(sapply(d[,1:2], as.numeric )) 
colnames(d)[3]<-'Watermass'
# Produce a dataframe

ggplot(data=d,aes(x=log_chla,y=log_spec_abs))+
#geom_point(size=5,alpha=1)+
geom_point(aes(shape=Watermass,col=Watermass,fill=Watermass),size=5,alpha=1)+ # Plot sample points 
scale_fill_manual("Water mass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_colour_manual("Water mass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Water mass", values = c("Arctic" = 22, "Atlantic" = 21, "Barents" = 24, "Coastal" = 23, "Melt" = 25))+
geom_text(label=c(1:40),col='black',size=3)+ # Label points with unique numbers 
geom_abline(intercept=mod$coef[1],slope=mod$coef[2],size=2,linetype='dashed')+
ylab( expression(paste('log'[10],'(',italic('a'[ph]),italic(''^'*'),'(675))',' (',m^2,'mg'^-1,')')) )+
xlab( expression(paste('log'[10],'(',italic('<Chl>'),')',' (mg m'^-3,')' )))+
theme(legend.position='right',
axis.line.x=element_line(),
axis.text.x=element_text(size=12),
axis.title.x=element_text(size=16),
axis.text.y=element_text(size=12),
axis.title.y=element_text(size=16),
legend.title=element_text(size=16),
legend.text=element_text(size=12),
axis.line.y=element_line(),
panel.background=element_blank())
# Produce a plot

# Take log_spec_abs as an index of packaging
Packaging_index<-1/log_spec_abs


diatom_masses<-rowSums(Biomasses[,5:55])/rowSums(Biomasses[,5:110])
dinoflagellate_masses<-(rowSums(Biomasses[,56:98])/(rowSums(Biomasses[,5:110])))
ciliate_masses<-(rowSums(Biomasses[,99:110])/(rowSums(Biomasses[,5:110])))
# Relative contribution of diatoms/different groups to microphytoplankton


plot(Environment$PAR[which(diatom_masses>.5)],Packaging_index[which(diatom_masses>.5)])
summary(lm(Packaging_index[which(diatom_masses>.5)]~Environment$PAR[which(diatom_masses>.5)]))

plot(Environment$PAR[which(Environment$Depth>20)],Packaging_index[which(Environment$Depth>20)])
summary(lm(Packaging_index[which(Environment$Depth>20)]~Environment$PAR[which(Environment$Depth>20)]))

plot(Environment$PAR[which(Environment$Depth<20)],Packaging_index[which(Environment$Depth<20)])
summary(lm(Packaging_index[which(Environment$Depth<20)]~Environment$PAR[which(Environment$Depth<20)]))

summary(lm(Packaging_index[which(Environment$Depth<20)]~Environment$PAR[which(Environment$Depth<20)]
+diatom_masses[which(Environment$Depth<20)]))


# Diatom-rich and sub-surface subsets do not reveal interesting dependencies of packaing index upon PAR

# Point 5
# Are the results in the presented manuscript mostly driven by diatoms?

# This can be tested by excluding diatom taxa from the cell counts data and re-performing analyses, 
# to see whether they differ substantialy from the original analyses based on cell counts.

Cell_data<- cbind(log10(Cell_counts[,5:110]+1),log10((Cell_counts[,111:114]*1000)+1))
Cell_data_without_diatoms<-Cell_data[,-c(1:51)]/rowSums(Cell_data[,-c(1:51)])

cell_rda_without_diatoms<-rda( Cell_data_without_diatoms ~ . , as.data.frame(explan$x))
# 36% of variance in non-diatom cell count structure is constrained by the environmental variables

for(i in 1: (dim(Environment[,-c(11:14)])[2]) ){ # for each environmental axis
	temporary_rda<-anova(rda( Cell_data_without_diatoms ~Env[,i],Z=Env[,-i])) # produce a Redundancy Analysis which partials-out all other axes
		print( c( 'proportion of variance explained', temporary_rda$V[1]/sum(temporary_rda$V), Test_axes[i] ) )
		print( c( 'p-value', temporary_rda$"Pr(>F)" [1], Test_axes[i] ) )
		# Print the explained variance and p-value derived from permutation 
}
# Pure effect statistics show that PC1, PC2, PC3, PC5, are significant explainers of cell count structure without
# diatoms

anova(cell_rda_without_diatoms, by="term")
# Term-wise effect statistics agree
anova(cell_rda_without_diatoms, by="margin")
# Marginal effect statistics agree

# Now we will produce a triplot for visual inspection


hydrog_vector<-rep('Barents',dim(Environment)[1] )
hydrog_vector[which(Environment$temp <0 & Environment$salinity <34.8)]<- 'Arctic' #Arw
hydrog_vector[which(Environment$temp >0 & Environment$temp <3 & Environment$salinity <34.4)]<- 'Melt' 
hydrog_vector[which(Environment$temp >3 & Environment$salinity <34.8)]<- 'Coastal' #Cw
hydrog_vector[which(Environment$temp >3 & Environment$salinity >34.8)]<- 'Atlantic' #AtW
# Define a vector of the Watermass classifications of each sample

d<-as.data.frame(cbind(scores(cell_rda_without_diatoms,display='sites'),
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



ggplot(data=d,aes(x=RDA1,y=RDA2,col=Watermass))+ # Plot with cell RDA data and watermass properties
ylim(c(-1,1))+
xlim(c(-1,1))+ # Set axis limits
scale_x_continuous(expand = c(0.1, 0.1)) + scale_y_continuous(expand = c(0.15, 0.15))+ 
geom_segment(data=cbind(cell_segments[cell_var,]),aes(x=V1,y=V2,xend=V3,yend=V4),size=2, 
alpha=1,color='gray50')+ # Plot vectors for environmental axes
geom_point(aes(shape=Watermass,col=Watermass,fill=Watermass),size=5,alpha=1)+ # Plot sample points 
scale_fill_manual("Water mass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_colour_manual("Water mass", values = c("Arctic" = "blue", "Atlantic" = "red", "Barents" = 'purple', "Coastal" = 'darkorange', "Melt" = 'cyan'))+
scale_shape_manual("Water mass", values = c("Arctic" = 22, "Atlantic" = 21, "Barents" = 24, "Coastal" = 23, "Melt" = 25))+
geom_text(label=c(1:40),col='black',size=3)+ # Label points with unique numbers 
geom_text(data=cell_segments[cell_var,]*1,aes(x=V3,y=V4),label=c(colnames(explan$x)[cell_var]),color='black',
size=5)+ # Add environmental axis labels 
theme(legend.position='right',
axis.line.x=element_line(),
axis.text.x=element_text(size=12),
axis.title.x=element_text(size=16),
axis.text.y=element_text(size=12),
axis.title.y=element_text(size=16),
legend.title=element_text(size=16),
legend.text=element_text(size=12),
axis.line.y=element_line(),
panel.background=element_blank())

# Which species are controlling variation along RDA-1 in the analysis without diatoms?
cell_rda_without_diatoms$CCA$v[,1][order(cell_rda_without_diatoms$CCA$v[,1])]
# The most extreme values will tell you which species are most distinct along this axis


# Point 6
# Are the results of the redundancy analyses in the original manuscript significantly explained
# by variation in overall biomass?

Biomasses<-read.csv('Biomasses_40.csv')
Overall_biomass<-rowSums(Biomasses[,-c(1:4)])

Biomasses<-read.csv('Biomasses_40.csv')
Environment<-read.csv('Environment_40.csv')
Absorption_spectra<-read.csv('Absorption_spectra_40.csv')
Cell_counts<-read.csv('Cell_counts_40.csv')
Pigments<-read.csv('Pigments_40.csv')
# Load the requisite datasets

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
# Run analyses


summary(lm(summary(cell_rda)$sites[,1]~Overall_biomass))
# significant but weak relationship of RDA-1 cells to overall biomass
summary(lm(summary(cell_rda)$sites[,2]~Overall_biomass))
# no clear dependency of RDA-2 cells on overall biomass

summary(lm(summary(pigment_rda)$sites[,1]~Overall_biomass))
# significant but weak relationship of RDA-1 pigments to overall biomass
summary(lm(summary(pigment_rda)$sites[,2]~Overall_biomass))
# no clear dependency of RDA-2 pigments on overall biomass

summary(lm(summary(optics_rda)$sites[,1]~Overall_biomass))
# significant but weak relationship of RDA-1 absorption spectra to overall biomass
summary(lm(summary(optics_rda)$sites[,2]~Overall_biomass))
# no clear dependency of RDA-2 absorption spectra on overall biomass

# Point 7
# Point 7 did not require new code to answer

# Point 8 
# does scaling the environmental input variables before analysis in RDA change the results?
# (It shouldn't because correlation structure is scaled in vegan's rda function)

cell_rda<-rda( Cell_data ~ . ,Environment[,-c(11:14)])

Scaled_Environment<-Environment[,-c(11:14)]
for(i in 1:(dim(Scaled_Environment)[2]) ){
Scaled_Environment[,i]<-(Scaled_Environment[,i]/sd(Scaled_Environment[,i]))-mean(Scaled_Environment[,i])
}
# scale environmental variables by their variances

cell_rda_scaled_environment<-rda(Cell_data ~ .,Scaled_Environment)
# reperform analysis

plot(cell_rda) ; dev.new() ; plot(cell_rda_scaled_environment)
# You can satisfy yourself the results are arithemtically near-identical

# Point 9
# Can we be certain that seasonal succession and Atlantic influence drive significant structure in the
# RDA analyses?

Julian<-Biomasses$Julian
# Julian day of the year that biomass samples were taken
Order<-order(Julian) # you can choose to order by diatom_masses instead if you wish
J<-(Julian/365)*12 # Quick conversion of Julian day to months (rough)

summary(lm(J ~ summary(cell_rda)$sites[,1] + summary(cell_rda)$sites[,2]))
# The month the samples were collected in is a strong function of both 
# leading RDA axes in the cell count data

summary(lm(J ~ summary(pigment_rda)$sites[,1] + summary(pigment_rda)$sites[,2]))
# The month the samples were collected in is a moderate function of RDA-1, but is 
# not clearly related to RDA-2

summary(lm(J ~ summary(optics_rda)$sites[,1] + summary(optics_rda)$sites[,2]))
# The month the samples were collected in is a moderate function of RDA-1, but is 
# not clearly related to RDA-2

# How can we tell whether Atlantic identify of samples is a significant variable structuring the
# distribution of site scores in the redundancy analyses? 
# We can use multiple anova, to see whether there is a separation of means by water mass group

Atlantic_vector<-hydrog_vector
Atlantic_vector[which(Atlantic_vector != 'Atlantic')]<-'Not_Atlantic'
# Atlantic_vector[which(hydrog_vector == 'Coastal')]<-'Coastal'
# Uncomment this line if you want to consider coastal water mass samples independently
# Atlantic_vector[which(hydrog_vector == 'Coastal')]<-'Atlantic'
# Uncomment this line if you wish to groap coastal water mass samples with Atlantic samples

cell_manova_data<-as.data.frame(cbind(summary(cell_rda)$sites,Atlantic_vector))
colnames(cell_manova_data)<-c('RDA1','RDA2','RDA3','RDA4','RDA5','RDA6','WaterMass')
cell_manova_data$WaterMass<-as.character(cell_manova_data$WaterMass)
for(i in 1:6){
	cell_manova_data[,i]<-as.numeric(as.character(cell_manova_data[,i]))
}
rownames(cell_manova_data)<-NULL
cell_manova<-manova( cbind(RDA1,RDA2,RDA3,RDA4,RDA5,RDA6) ~ WaterMass, data=cell_manova_data)
summary.aov(cell_manova)
# This result suggests that Atlantic status significantly determines point dispersion along 
# axes RDA1 and RDA2
# The p values are 0.03 and 0.004


pigment_manova_data<-as.data.frame(cbind(summary(pigment_rda)$sites,Atlantic_vector))
colnames(pigment_manova_data)<-c('RDA1','RDA2','RDA3','RDA4','RDA5','RDA6','WaterMass')
pigment_manova_data$WaterMass<-as.character(pigment_manova_data$WaterMass)
for(i in 1:6){
	pigment_manova_data[,i]<-as.numeric(as.character(pigment_manova_data[,i]))
}
rownames(pigment_manova_data)<-NULL
pigment_manova<-manova( cbind(RDA1,RDA2,RDA3,RDA4,RDA5,RDA6) ~ WaterMass, data=pigment_manova_data)
summary.aov(pigment_manova)
# This result suggests that Atlantic status significantly structures variance over
# RDA1 but not RDA2 (p = 0.003 and p = 0.07)

optics_manova_data<-as.data.frame(cbind(summary(optics_rda)$sites,Atlantic_vector))
colnames(optics_manova_data)<-c('RDA1','RDA2','RDA3','RDA4','RDA5','RDA6','WaterMass')
optics_manova_data$WaterMass<-as.character(optics_manova_data$WaterMass)
for(i in 1:6){
	optics_manova_data[,i]<-as.numeric(as.character(optics_manova_data[,i]))
}
rownames(optics_manova_data)<-NULL
optics_manova<-manova( cbind(RDA1,RDA2,RDA3,RDA4,RDA5,RDA6) ~ WaterMass, data=optics_manova_data)
summary.aov(optics_manova)
# This result shows that Atlantic status is not significantly related to variance over
# RDA1 and RDA2, (p = 0.13 and p = 0.6)
# However, we note on inspecting original manuscript Figure 6C that there is an over-lap of Atlantic
# samples with Coastal samples in RDA1-2 space.
# If Coastal samples are treated as their own class, then water mass variation is significantly 
# related to RDA1 (p = 0.00014)
# and if Coastal samples are grouped with Atlantic samples, then water mass is very
# clearly related to RDA1 (p = 9.4e-05). 


# What about if we try to remove the effects of diachronous sampling? 
# For RDA axes that are significantly related to the time-of-sampling we need to 
# perform a linear model of axis score against time-of-sampling, and then extract the 
# residual.
# The MANOVA analysis can then be re-performed to see whether signal remains. 

RDA1r<-lm(summary(cell_rda)$sites[,1]~J)$residuals
RDA2r<-lm(summary(cell_rda)$sites[,2]~J)$residuals

cell_manova_data<-as.data.frame(cbind(RDA1r,RDA2r,Atlantic_vector))
colnames(cell_manova_data)<-c('RDA1','RDA2','WaterMass')
cell_manova_data$WaterMass<-as.character(cell_manova_data$WaterMass)
for(i in 1:2){
	cell_manova_data[,i]<-as.numeric(as.character(cell_manova_data[,i]))
}
rownames(cell_manova_data)<-NULL
cell_manova<-manova( cbind(RDA1,RDA2) ~ WaterMass, data=cell_manova_data)
summary.aov(cell_manova)


RDA1r<-lm(summary(pigment_rda)$sites[,1]~J)$residuals
RDA2r<-summary(pigment_rda)$sites[,2]

pigment_manova_data<-as.data.frame(cbind(RDA1r,RDA2r,Atlantic_vector))
colnames(pigment_manova_data)<-c('RDA1','RDA2','WaterMass')
pigment_manova_data$WaterMass<-as.character(pigment_manova_data$WaterMass)
for(i in 1:2){
	pigment_manova_data[,i]<-as.numeric(as.character(pigment_manova_data[,i]))
}
rownames(pigment_manova_data)<-NULL
pigment_manova<-manova( cbind(RDA1,RDA2) ~ WaterMass, data=pigment_manova_data)
summary.aov(pigment_manova)


RDA1r<-lm(summary(optics_rda)$sites[,1]~J)$residuals
RDA2r<-summary(optics_rda)$sites[,2]

optics_manova_data<-as.data.frame(cbind(RDA1r,RDA2r,Atlantic_vector))
colnames(optics_manova_data)<-c('RDA1','RDA2','WaterMass')
optics_manova_data$WaterMass<-as.character(optics_manova_data$WaterMass)
for(i in 1:2){
	optics_manova_data[,i]<-as.numeric(as.character(optics_manova_data[,i]))
}
rownames(optics_manova_data)<-NULL
optics_manova<-manova( cbind(RDA1,RDA2) ~ WaterMass, data=optics_manova_data)
summary.aov(optics_manova)


# Now we move on to the minor comments section


# Let us produce a table of the most abdundant microphytoplankton taxa in each water mass

setwd('D:/Documents/Final_nutrients_sian')
Cell_counts<-read.csv('Cell_counts_40.csv')
Biomasses<-read.csv('Biomasses_40.csv')
Environment<-read.csv('Environment_40.csv')

hydrog_vector<-rep('Barents',dim(Environment)[1] )
hydrog_vector[which(Environment$temp <0 & Environment$salinity <34.8)]<- 'Arctic' #Arw
hydrog_vector[which(Environment$temp >0 & Environment$temp <3 & Environment$salinity <34.4)]<- 'Melt' #Mw
hydrog_vector[which(Environment$temp >3 & Environment$salinity <34.8)]<- 'Coastal' #Cw
hydrog_vector[which(Environment$temp >3 & Environment$salinity >34.8)]<- 'Atlantic' #AtW
# Define a vector of the Watermass classifications of each sample


colSums(Cell_counts[which(hydrog_vector=='Arctic'),5:110])[
order(colSums(Cell_counts[which(hydrog_vector=='Arctic'),5:110]))]

colSums(Cell_counts[which(hydrog_vector=='Atlantic'),5:110])[
order(colSums(Cell_counts[which(hydrog_vector=='Atlantic'),5:110]))]

colSums(Cell_counts[which(hydrog_vector=='Barents'),5:110])[
order(colSums(Cell_counts[which(hydrog_vector=='Barents'),5:110]))]

colSums(Cell_counts[which(hydrog_vector=='Coastal'),5:110])[
order(colSums(Cell_counts[which(hydrog_vector=='Coastal'),5:110]))]

colSums(Cell_counts[which(hydrog_vector=='Melt'),5:110])[
order(colSums(Cell_counts[which(hydrog_vector=='Melt'),5:110]))]

colSums(Biomasses[which(hydrog_vector=='Arctic'),5:110])[
order(colSums(Biomasses[which(hydrog_vector=='Arctic'),5:110]))]

colSums(Biomasses[which(hydrog_vector=='Atlantic'),5:110])[
order(colSums(Biomasses[which(hydrog_vector=='Atlantic'),5:110]))]

colSums(Biomasses[which(hydrog_vector=='Barents'),5:110])[
order(colSums(Biomasses[which(hydrog_vector=='Barents'),5:110]))]

colSums(Biomasses[which(hydrog_vector=='Coastal'),5:110])[
order(colSums(Biomasses[which(hydrog_vector=='Coastal'),5:110]))]

colSums(Biomasses[which(hydrog_vector=='Melt'),5:110])[
order(colSums(Biomasses[which(hydrog_vector=='Melt'),5:110]))]

# What is the minimum number of cells recorded in a sample?

range(rowSums(Cell_counts[,5:110]))


