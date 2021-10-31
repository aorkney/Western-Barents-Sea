# This script will perform Redundancy Analyses, 
# relating phytoplankton community composition and environmental conditions,
# including an additional analysis for which environmental conditions are
# decomposed by PCA to eliminate co-linearity. 
# Pure-effect, Term-effect and Marginal-effect significance testing
# will then be applied to determine which environmental variables
# (or which axes of variance in the environmental data) explain
# structure in the phytoplankton community composition datasets. 

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

cell_rda<-rda( Cell_data ~ . ,Environment[,-c(11:14)])
optics_rda<-rda( Absorption_data ~ . ,Environment[,-c(11:14)])
pigment_rda<-rda( Pigment_data ~ . ,Environment[,-c(11:14)])
# Perform Redundancy Analyses 
Test_axes<- colnames(Environment[-c(11:14)])
Env<-Environment[,-c(11:14)]

# Uncomment the following lines if you want to perform 
# PCA decomposition on the Environmental data

# !!!
# explan<-prcomp(Environment[,-c(11:14)],scale=T) 
# PCA decomposition
# cell_rda<-rda( Cell_data ~ . ,as.data.frame(explan$x))
# optics_rda<-rda( Absorption_data ~ . ,as.data.frame(explan$x))
# pigment_rda<-rda( Pigment_data ~ . ,as.data.frame(explan$x))
# Run analyses
# Test_axes<- colnames(explan$rotation)
# Env<-explan$x
# !!!

# We will begin by assessing the 'Pure-effects' associated with each environmental variable 
# The significance is assessed by permutation, and your results may differ by a small amount from those reported
# in any publication or manuscript associated with this code. 

for(i in 1: (dim(Environment[,-c(11:14)])[2]) ){ # for each environmental axis
	temporary_rda<-anova(rda( Cell_data ~Env[,i],Z=Env[,-i])) # produce a Redundancy Analysis which partials-out all other axes
		print( c( 'proportion of variance explained', temporary_rda$V[1]/sum(temporary_rda$V), Test_axes[i] ) )
		print( c( 'p-value', temporary_rda$"Pr(>F)" [1], Test_axes[i] ) )
		# Print the explained variance and p-value derived from permutation 
}

for(i in 1: (dim(Environment[,-c(11:14)])[2]) ){
	temporary_rda<-anova(rda( Pigment_data ~Env[,i],Z=Env[,-i]))
		print( c( 'proportion of variance explained', temporary_rda$V[1]/sum(temporary_rda$V), Test_axes[i] ) )
		print( c( 'p-value', temporary_rda$"Pr(>F)" [1], Test_axes[i] ) )
}

for(i in 1: (dim(Environment[,-c(11:14)])[2]) ){
	temporary_rda<-anova(rda( Absorption_data ~Env[,i],Z=Env[,-i]))
		print( c( 'proportion of variance explained', temporary_rda$V[1]/sum(temporary_rda$V), Test_axes[i] ) )
		print( c( 'p-value', temporary_rda$"Pr(>F)" [1], Test_axes[i] ) )
}



# Now we will calculate 'Term-wise' Effects, for which the significance of each term/environmental axis is calculated
# in series.

Cell_term_effects<-anova(cell_rda, by="term")
Optics_term_effects<-anova(optics_rda, by="term")
Pigment_term_effects<-anova(pigment_rda, by="term")

print(colnames(Env))
Cell_term_effects$Variance[1:10]/sum(Cell_term_effects$Variance) # Variance
Cell_term_effects$"Pr(>F)"[1:10] # p-values
Optics_term_effects$Variance[1:10]/sum(Optics_term_effects$Variance)
Optics_term_effects$"Pr(>F)"[1:10]
Pigment_term_effects$Variance[1:10]/sum(Pigment_term_effects$Variance)
Pigment_term_effects$"Pr(>F)"[1:10]

# As these statistics are derived from permutation your results may differ from those reported in any manuscripts
# that use this code.


# Now we will calculate 'Marginal' Effects, for which the significance of each term/environmental axis is calculated
# in a multivariate context. 

Cell_margin_effects<-anova(cell_rda, by="margin")
Optics_margin_effects<-anova(optics_rda, by="margin")
Pigment_margin_effects<-anova(pigment_rda, by="margin")

print(colnames(Env))
Cell_margin_effects$Variance[1:10]/sum(Cell_margin_effects$Variance) # Variance
Cell_margin_effects$"Pr(>F)"[1:10] #p-value
Optics_margin_effects$Variance[1:10]/sum(Optics_margin_effects$Variance) 
Optics_margin_effects$"Pr(>F)"[1:10] 
Pigment_margin_effects$Variance[1:10]/sum(Pigment_margin_effects$Variance) 
Pigment_margin_effects$"Pr(>F)"[1:10] 

# As these statistics are derived from permutation your results may differ from those reported in any manuscripts
# that use this code.


