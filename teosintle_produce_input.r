# BiocManager::install("LEA")
# install.packages("gradientForest", repos="http://R-Forge.R-project.org",dependendices=T)
# devtools::install_github("AndiKur4/MaizePal")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)
library(LEA);
library(adegenet);
library(maps);
library(dismo);
library(gplots);
library(raster);
library(gradientForest);
library(gdistance);
library(geosphere);
library(MaizePal)
library(tidyverse);
library(igraph);
library(ggridges);
library(UpSetR);
library(here);
library(rasterVis);

source("teosintle/datasets/code/climate_change_functions.R")
# create all directories that will be used



#' Load datasets
load("teosintle/datasets/input/species_input.R")
# Since adegenet shows genotypes as pairs of alleles, we need to select columns
# containing the frequencies of the p alleles (SNP1.G,SNP1.C; SNP2.T; SNP3.C,SNP3.A).
# The LEA package can generate the input data from different datasets (plink, lfmm,
# data.frame), but it is important to remove fixed and missing SNPs. Check the LEA tutorial
# to understand how to run the analyses and see example for input formats.
geno <- species_input$genind$tab
geno <- geno[,seq(1,ncol(geno)-1,2)]
# Obtain the coordinates and names of populations
coords <- species_input$genind$other
pops <- rownames(coords)
# In the next sections we create datasets for two clusters, but the user should be able to
# set the number of clusters accordingly (see Functions).
# remove fixed SNPs and set NA to 9 (needed by snmf)
fixed_SNP <- apply(geno,2,sd,na.rm=T)
geno <- geno[,-which(fixed_SNP==0)]
geno[which(is.na(geno))] <- 9




project <- load.snmfProject("data/teosintle_species.snmfProject")
i = 2
ques <- Q(object = project,K = i)
#define individual that have a higher contribution to a given genetic group (50% split) and
# get the name of each individual
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
# The next code will extract the information (lon, lat, bio_1) of each individual from the
# population data and add it to the data frame order inds.
order_inds <- data.frame(rownames(geno),ques, pop=NA,lon=NA,lat=NA,bio_1=NA)
names(order_inds)[1:3] <- c("inds","g1","g2")
#coords and pops were defined at the beginning.
coords <- coords[,c("longitude","latitude","bio_1")]
# For each individual, set the name of the pop, and the longitude, latitude and bio_1,
# where it grows. The next loop takes the info of each population, and searches the
# individuals that belong to the p



for(i in 1:length(pops)){
  tpop <- pops[i] #one pop at the time
  temp <- grep(paste(tpop,"_",sep = ""),order_inds$inds) #find all inds that belong to pop i.
  order_inds[temp,"pop"] <- rownames(coords[tpop,]) #add name of pop (temp has the
  #index of all individuals belonging to the population
  order_inds[temp,"lon"] <- coords[tpop,"longitude"]#add longitude
  order_inds[temp,"lat"] <- coords[tpop,"latitude"]
  order_inds[temp,"bio_1"] <- coords[tpop,"bio_1"]
}
#finally we plot the ancestry proportions based on the annual mean temperature (BIO1) at
#which populations grow.
warm_col <- "red"
cold_col <- "blue"
my.colors <- c(warm_col,cold_col)
bio_1 <- order_inds[order(order_inds$bio_1),]
bio_1 <- as.matrix(bio_1[,c("g1","g2")])
par(mai=c(0.5,1,0.6,1))
jpeg("admixture.jpeg", width = 1000, height = 600)
barplot(t(bio_1),
        col=my.colors,
        border=NA, 
        width = 1,
        space = 0,
        main="Individuals (orderedtemperature)",
        ylab="Ancestry proportions (k=2)",xlab="",xaxt="n")
dev.off()





# Run SDM per population

# get all individuals belonging a to genetic group (%) and select the column pop, finally,
# get only the unique data
# first identify the groups with function which (if you increase the ancestry threshold above
# 0.5 you can remove potential admixed populations, and if you decrease the threshold you
# can include admixed populations.
g1 <- order_inds[which(order_inds$g1 >= 0.5),]
g2 <- order_inds[which(order_inds$g2 >= 0.5),]
# get the names of populations
coord_g1 <- unique(g1[,"pop"])
coord_g2 <- unique(g2[,"pop"])
# repeat if you have more clusters
coord_g1 <- coords[coord_g1,c("longitude","latitude")]
coord_g2 <- coords[coord_g2,c("longitude","latitude")]
# change the names of the columns to ID, lon and lat that will be used for the input of
# maxent
coord_g1 <- data.frame(ID="g1",lon= coord_g1 $longitude,lat= coord_g1$latitude)
coord_g2 <- data.frame(ID="g2",lon= coord_g2 $longitude,lat= coord_g2$latitude)
# one of the g2 populations grow within the opposite genetic group. They might be outlier
# populations, so we remove them to avoid an issue with the convex hull model.
coord_g2 <- coord_g2[-which(coord_g2$lat>19.75),]
apply_buffer = FALSE # We can apply a 'buffer' around the occurrence localities; if FALSE a convexHull is estimated.
# get all known teosintes coordinates from CONABIO and select the Mexicana group
mexicana <- read.table("teosintle/datasets/input/coordinates.csv",header = T,sep = ",")
mexicana <- mexicana[grep("mexicana",mexicana$Taxa),]
names(mexicana) <- c("ID","lon","lat")
# we use convex hull models to identify populations that are putatively within the
# distribution of each cluster to increase the number of populations & accuracy SDMs
# obtain the convex hull
size_buffer = 3 # Size of buffer
if(apply_buffer) ch_g1 <- buffer(SpatialPoints(coord_g1[,c("lon","lat")]),width=size_buffer)
if(!apply_buffer) ch_g1 <- convHull(coord_g1[,c("lon","lat")])
if(apply_buffer) ch_g2 <- buffer(SpatialPoints(coord_g2[,c("lon","lat")]),size_buffer)
if(!apply_buffer) ch_g2 <- convHull(coord_g2[,c("lon","lat")])
# predict populations belonging to the models, this function searches all population from
# the CONABIO dataset that grow inside the mexicana genetic groups.
if(apply_buffer) inside <- sp::over(SpatialPoints(mexicana[,2:3]),ch_g1)
if(!apply_buffer) inside <- predict(ch_g1,mexicana[,2:3])
g1 <- mexicana[which(inside==1),2:3]
g1 <- data.frame(ID="g1", g1)
#same for g2
if(apply_buffer) inside <- sp::over(SpatialPoints(mexicana[,2:3]),ch_g2)
if(!apply_buffer) inside <- predict(ch_g2,mexicana[,2:3])
g2 <- mexicana[which(inside==1),2:3]
g2 <- data.frame(ID="g2", g2)
#get coordinates for each genetic cluster (g1, g2, all)
g1 <- rbind(g1,coord_g1) #this function combines the coordinates from the conabio and geno populations
g2 <- rbind(g2,coord_g2)
maxt <- rbind(g1,g2) # we combine all data (all the popualtions)

# create the maxent input containing each of the genetic groups and the the entire
# coordinates
all <- maxt
all$ID <- "all" #we change the name from g1 and g2 to all obtain the SDM of the combined dataset
maxt <- rbind(maxt,all) #we combine all data (g1,g2, and all)
write.table(maxt,file = "data/teosintle_maxent_input.csv",row.names = F,col.names = T,sep = ",")











