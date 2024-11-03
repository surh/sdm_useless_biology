# Original Copyright (c) 2021 Santiago Ramírez Barahona under MIT License
# Modifications Copyright (c) 2024 Sur Herrera Paredes under GPL3 License




#### Step 1. Initial configuration ####
#BiocManager::install("LEA")
#install.packages("gradientForest", repos="http://R-Forge.R- project.org",dependendices=T)
#devtools::install_github("AndiKur4/MaizePal")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)
library(LEA);
library(adegenet);
library(maps); 
library(dismo); 
library(gplots); 
library(raster); 
library(gradientForest); 
library(gdistance); 
library(geosphere);library(MaizePal) 
library(tidyverse);
library(ggplot2);
library(igraph);
library(ggridges);
library(UpSetR);
library(here);
library(rasterVis)
library(ggridges)




library(tidyverse)



# setwd(here::here())
# create a raster mask with values set to zero, 
#but we the same extent and resolution as
#the bioclimatic variables uses.

mask <- raster::raster("../Climate-Change-Genomics/datasets/input/present/bio_1.asc") %>% replace(.,.,0)
# source functions that will be used
source("datasets/code/climate_change_functions.R") # create all directories that will be used 
if(!dir.exists("datasets/output")){
dir.create("datasets/output") }
if(!dir.exists("datasets/conservation")){ dir.create("datasets/conservation")
} if(!dir.exists("datasets/output/admixture")){
  dir.create("datasets/output/admixture")}
if(!dir.exists("datasets/maxent")){ dir.create("datasets/maxent")}
#### Step 1. Run admixture analyses #### 
# Load genomic dataset and obtain genotypes ($tab slot).
load("datasets/input/species_input.R")
# Since adegenet shows genotypes as pairs of alleles, we need to select columns containing the frequencies of the p alleles (SNP1.G,SNP1.C; SNP2.T; SNP3.C,SNP3.A). # The LEA package can generate the input data from different datasets (plink, lfmm, data.frame), but it is important to remove fixed and missing SNPs. Check the LEA tutorial to understand how to run the analyses and see example for input formats.
geno <- species_input$genind$tab
geno <- geno[,seq(1,ncol(geno)-1,2)]
# Obtain the coordinates and names of populations
coords <- species_input$genind$other
pops <- rownames(coords)
# In the next sections we create datasets for two clusters, but the user should be able to set the number of clusters accordingly (see Functions).
# remove fixed SNPs and set NA to 9 (needed by snmf)
fixed_SNP <- apply(geno,2,sd,na.rm=T) 
geno <- geno[,-which(fixed_SNP==0)] 
geno[which(is.na(geno))] <- 9
# write geno object for admixture analysis
write.geno(geno,output.file = "datasets/output/admixture/species.geno")
# create project to run the admixture analysis and run the analysis for a K number of groups (for more details on the method see the vignette of the LEA package)
project = NULL
if(FALSE){
  project = snmf("datasets/output/admixture/species.geno",
                 K = 1:10, #set number of K
                 entropy = TRUE,
                 repetitions = 1, #set number of repetitions
                 project = "new")
  # The LEA package creates a file project that contains the information of the analysis and can be loaded for future analyses
}
warm_col=MaizePal::maize_pal("JimmyRed",4)[2] 
cold_col=MaizePal::maize_pal("MaizAzul",6)[6]
project <- load.snmfProject("datasets/output/admixture/species.snmfProject")
# the Q function generates a matrix showing the ancestry proportions of each individual. The number of columns corresponds to K, and each column contains the proportions of an individual to each population. The highest proportion indicates the membership to a group. i=2
ques <- Q(object = project,K = i)
#define individual that have a higher contribution to a given genetic group (50% split) and get the name of each individual
g1 <- which(ques[,1]>=0.5)
g2 <- which(ques[,2]>=0.5)
g1 <- rownames(geno)[ g1]
g2 <- rownames(geno)[ g2] # if you have three clusters add a new group
# we know that “Ixtlan” grows in the North (make sure that g1 belongs to G1-N).
if(length(grep("Ixtlan", g1)) ==0){ 
temp <- g1
g1 <- g2
g2 <- temp
ques <- ques[,c(2,1)]
rm(temp)
}
# The next code will extract the information (lon, lat, bio_1) of each individual from the population data and add it to the data frame order inds.
order_inds <- data.frame(rownames(geno),ques, pop=NA,lon=NA,lat=NA,bio_1=NA) 
names(order_inds)[1:3] <- c("inds","g1","g2")
#coords and pops were defined at the beginning.
coords <- coords[,c("longitude","latitude","bio_1")]
# For each individual, set the name of the pop, and the longitude, latitude and bio_1, where it grows. The next loop takes the info of each population, and searches the individuals that belong to the population (grep), finally in the rows where those individuals are, it adds the information
for(i in 1:length(pops)){
  tpop <- pops[i] #one pop at the time
  temp <- grep(paste(tpop,"_",sep = ""),order_inds$inds) #find all inds that belong to pop i. 
  order_inds[temp,"pop"] <- rownames(coords[tpop,]) # add name of pop (temp has the index of all individuals belonging to the population 
  order_inds[temp,"lon"] <- coords[tpop,"longitude"] #add longitude 
  order_inds[temp,"lat"] <- coords[tpop,"latitude"] 
  order_inds[temp,"bio_1"] <- coords[tpop,"bio_1"]
}
#finally we plot the ancestry proportions based on the annual mean temperature (BIO1) at which populations grow.
my.colors <- c(warm_col,cold_col)
bio_1 <- order_inds[order(order_inds$bio_1),]
bio_1 <- as.matrix(bio_1[,c("g1","g2")])
par(mai=c(0.5,1,0.6,1)) 
barplot(t(bio_1),col=my.colors,border=NA,main="Individuals (ordered by mean
 temperature)",ylab="Ancestry proportions (k=2)",xlab="",xaxt="n")
