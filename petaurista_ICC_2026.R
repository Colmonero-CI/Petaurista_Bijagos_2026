##Author: Ivo Colmonero-Costeira
##Date: 2026-09
##Discription: This script contains the code used to analyze the str and mitochondrial DNA data 
##from Cercopithecus petaurista of the Bijagós Archpelago, Guinea-Bissau. doi:10.1038/s41598-026-69922-4

#Example files are deposited in the GitHub page

##### Packages #####

library(hierfstat)
library(PopGenReport)
library(adegenet)
library(ape)
library(vegan)

#load the data
data.t <-read.genetable("petaurista_str_ICC_2026.csv", pop=1, ind=2, long=3, lat=4, other.min = 5, other.max = 6,
                        oneColPerAll=FALSE, ploidy=2, ncode =3, NA.char = "0") #Change NA.char according to your database!

#check the data
is.genind(data.t) 
summary(data.t)

####################
##Allelic Richness##
####################

library(hierfstat)

allele_rich<-allelic.richness(data.t) #With rearefaction
allele_rich

####################
### Pairwise FST ###
####################

library(hierfstat)

fst<-pairwise.neifst(data.t)
fst

###########
### PCA ###
###########

#retrieve the data in a tabular format, replace NA with the mean of the locus
data.pca <- tab(data.t, freq=FALSE, NA.method="mean") 

#compute the PCA
pca1 <-rda(data.pca, center=TRUE, scale=FALSE)

#check the ordination plot
summary(pca1)

###########
### RDA ###
###########

#remove samples with missing coordinates
data.pca_w_coords <- data.pca[!is.na(data.t$other$latlong[,1]),]
summary(data.pca_w_coords)

#extract the coordinates in x,y format
space.xy=data.t$other$latlong[!is.na(data.t$other$latlong[,1]),]

#Next we make the polinomials to fit the data
#start by centering the coordinates
x = scale(space.xy[,1], scale=FALSE)
y = scale(space.xy[,2], scale=FALSE)

#make the polinomials, in this case we are using a 3rd degree polynomial
space = poly(cbind(x,y), degree=3, raw=FALSE)

#rename the columns for better understanding
colnames(space) <- c("X","X2","X3","Y","XY","X2Y","Y2","XY2","Y3")
plot(space[,4] ~space[,1]) #check if projection is ok


#compute the RDA
ord.spa<-rda(data.pca_w_coords ~., data.frame(space), add=TRUE)

#check the ordination plot
plot(ord.spa)

#summary of the ordination plot
summary(ord.spa)

#adjusted R2
(R2.all.spa = RsquareAdj(ord.spa)$adj.r.squared) 

#now we will select the variables that are significant using a forward selection procedure with p critical value of 0.01 and R2 stoping criterion
stp.spa = ordistep(rda(data.pca_w_coords ~ 1, data.frame(space)), scope = formula(ord.spa), scale= FALSE, 
                  direction="forward", pstep = 9999, Pin = 0.01, R2scope = TRUE)

#retain the selected variables
(selected.spa = attributes(stp.spa$terms)$term.labels)

#change the space matrix to include only the selected variables
space = space[,selected.spa]

#re-do the db-RDA with the selected variables
ord.spa<-rda(data.pca_w_coords ~., data.frame(space), add=TRUE)

#now let's check the collinearity of the retained variables

vif.cca(ord.spa) ##remove variables with VIF >5
space=space[,c("X","X2", "Y", "XY")] ##change according to vif results

#re-do the db-RDA with the selected non-colinear variables
ord.spa<-rda(data.pca_w_coords ~., data.frame(space), add=TRUE)

#check the ordination plot
plot(ord.spa)

#summary of the ordination plot
summary(ord.spa)

#adjusted R2
(R2.all.spa = RsquareAdj(ord.spa)$adj.r.squared) 

#now let's check the significance of the db-RDA
a<-anova.cca(ord.spa, permutations = how(nperm=9999))
a

#check the significance of each axis
axis.test = anova.cca(ord.spa, by="axis", step=9999)

#percentage of variance explained by each axis
tot_inert = ord.spa$tot.chi
orig_eigs = ord.spa$CCA$eig
orig_inert = sum(orig_eigs)
all_perc = (orig_eigs/orig_inert)*100 

#################
### mtDNA MDS ###
#################

#Load the data
#I find it easier to just "copy and paste" data directly from the excel file to R.

camp<-as.matrix(read.table(pipe("pbpaste"), sep="\t", row.names = 1)) #add the DNA alignment
#No header!
#Include sample name!

#transform the data to compute genetic distances
camp<- as.alignment(camp)
camp<-as.DNAbin(camp)
Dgen <- dist.dna(camp, model="TN93", gamma=4, pairwise.deletion = TRUE) #TN93 model for C. petaurista

#compute the MDS
mds<-cmdscale(Dgen, eig = TRUE)
