#######################################################################################################
#                                                                                                     #
#                   ANALYSIS OF TAXONOMIC AND FUNCTIONALDIVERSITY                     
#                              IN FKNMS REEF FISHES                                      
#                                                                                                     #
#######################################################################################################

# Megan Hepner
# mhepner@mail.usf.edu 

#######################################################################################################
#                                        TABLE OF CONTENTS                                            #
#   Line 21:  Required libraries                                                                      #
#   Line 52:  Importing and formatting the data                                                       #                                                #
#   Line 67:  Calculate Species Richness 
#   Line 112: Calculate Effective Simpson Diversity 
#   Line 156: Calculate Effective Shannon Diversity 
#   Line 107: Functional dendrogram  
#   Line 234: Calculate Gowers distance                                                                  #                                                            #
#   Line 376: Characterizing bivariate relationships (SPLOMs, Mantel tests)                           #                                                                #
#                                                                                                     #
#######################################################################################################

library(ade4) #Calls: mantel.rtest, dudi.mix, disc 
library(FD) #Calls: gowdis -
library(vegan) #Calls: simpson, species richness, shannon, PERMANOVA
library(magrittr) #pipping 
library(secr) 
library(mgcv) #GAM 
library(base) #hclust - builds dendrograms 
library(clue) #cl_consensus - consensus dendogram is built 
library(tidyverse) 
#cl_ultrametric - dendrogram was scaled to 1 
#cl_dissimilarity - extract pairwaise distances using a matrix norm, the 2-norm

#######################################################################################################
#                                IMPORTING AND FORMATTING THE DATA                                    
#######################################################################################################

rm(list=ls())

#Import abundance from file
relative_abund=read.csv("~/Dropbox/R_Thesis_Project/data/relative_abundance.csv")

#Import species list from file
species_list=read.csv("~/Dropbox/R_Thesis_Project/data/species_code_species_common_name.csv")

#Import functional trait matrix from file
traits=read.csv("~/Dropbox/R_Thesis_Project/data/species_trait_matrix.csv")

#######################################################################################################
#                                    TAXANOMIC DIVERSITY                                           #
#######################################################################################################

#Calculate Species Richness 
Data2Plot <- array(dim=c(4,16))
for (whatarea in 1:4){ #subregions upper, middle, lower, DT
  for (whatyear in 1999:2014){
    
    thisdata1 = relative_abund[relative_abund$Year == whatyear,]
    thisdata2 = thisdata1[thisdata1$Subregion == whatarea,]
    thisdata3 = thisdata2$Relative.Abundance
    
    species_richness = specnumber(thisdata3) #specnumber from vegan 
    
    yearindex = whatyear - 1998
    Data2Plot[whatarea,yearindex] = species_richness
  }
}

Data2Plot=Data2Plot[,c(1:14,16)]
YearsToPlot = c(seq(from=1999,to=2012),2014)
PlotTitle = c("Upper Keys","Middle Keys","Lower Keys","Dry Tortugas")

par(new=FALSE)
tiff(filename = "plot1.tif",width=4, height=8, units="in", res=600, compression = "lzw")

par(cex.main = 1.4) #size of the main title
par(omi=c(.6,.7,.5,.2)) # outer margin areas
par(mar=c(2.2,2,2.8,1)) # inner margin areas  #bottom, left, top, right
par(mfrow=c(4,1)) # number of rows and columns

for (whatarea in 1:4){
  
  plot(x=YearsToPlot,y=Data2Plot[whatarea,],ylim=c(0,250),xlim=c(1998,2016),main=PlotTitle[whatarea],xlab="",ylab="",axes=FALSE)
  box()
  axis(side=1,at=c(1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2014,2016))
  axis(2)
  
}
mtext("Year",side=1,line=1,outer=T,adj=.5,cex=1.2)
mtext("Number of Species",outer=T,side=2,line=2.3,cex=1.2)
mtext("Species Richness", side = 3, outer = T, cex = 1.2)

dev.off()

write.csv(Data2Plot, file = "~/Dropbox/Thesis/Functional_Diversity/Results with four Subregions/species_richness.csv") #records species richness 

#Calculate Effective Simpson Diversity 
DataToPlot <- array(dim=c(4,16))
for (whatarea in 1:4){ #subregions upper, middle, lower, DT
  for (whatyear in 1999:2014){
    
    thisdata1 = relative_abund[relative_abund$Year == whatyear,]
    thisdata2 = thisdata1[thisdata1$Subregion == whatarea,]
    thisdata3 = thisdata2$Relative.Abundance
    
    effective_simpson = 1/(1-diversity(thisdata3, index = "simpson"))
    
    yearindex = whatyear - 1998
    DataToPlot[whatarea,yearindex] = effective_simpson
  }
}

DataToPlot=DataToPlot[,c(1:14,16)]
YearsToPlot = c(seq(from=1999,to=2012),2014)
PlotTitle = c("Upper Keys","Middle Keys","Lower Keys","Dry Tortugas")

par(new=FALSE)
tiff(filename = "plot2.tif",width=4, height=8, units="in", res=600, compression = "lzw")

par(cex.main = 1.4) #size of the main title
par(omi=c(.6,.7,.5,.2)) # outer margin areas
par(mar=c(2.2,2,2.8,1)) # inner margin areas  #bottom, left, top, right
par(mfrow=c(4,1)) # number of rows and columns

for (whatarea in 1:4){
  
  plot(x=YearsToPlot,y=DataToPlot[whatarea,],ylim=c(0,30),xlim=c(1998,2016),main=PlotTitle[whatarea],xlab="",ylab="",axes=FALSE)
  box()
  axis(side=1,at=c(1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2014,2016))
  axis(2)
  
}
mtext("Year",side=1,line=1,outer=T,adj=.5,cex=1.2)
mtext("Effective Number of Species",outer=T,side=2,line=2.3,cex=1.2)
mtext("Simpson Diversity", side = 3, outer = T, cex = 1.2)

dev.off()

write.csv(DataToPlot, file = "~/Dropbox/Thesis/Functional_Diversity/Results with four Subregions/effective_simpson.csv") 

#Calculate Shannon Effective Diversity 
DataPlot <- array(dim=c(4,16))
for (whatarea in 1:4){ #subregions upper, middle, lower, DT
  for (whatyear in 1999:2014){
    
effective_shannon = exp(diversity(thisdata3, index = "shannon"))

    yearindex = whatyear - 1998
    DataPlot[whatarea,yearindex] = effective_shannon
    
  }
}

DataPlot=DataPlot[,c(1:14,16)]
YearsToPlot = c(seq(from=1999,to=2012),2014)
PlotTitle = c("Upper Keys","Middle Keys","Lower Keys","Dry Tortugas")

par(new=FALSE)
tiff(filename = "plot3.tif",width=4, height=8, units="in", res=600, compression = "lzw")

par(cex.main = 1.4) #size of the main title
par(omi=c(.6,.7,.5,.2)) # outer margin areas
par(mar=c(2.2,2,2.8,1)) # inner margin areas  #bottom, left, top, right
par(mfrow=c(4,1)) # number of rows and columns

for (whatarea in 1:4){
  
  plot(x=YearsToPlot,y=DataPlot[whatarea,],ylim=c(5,15),xlim=c(1998,2016),main=PlotTitle[whatarea],xlab="",ylab="",axes=FALSE)
  box()
  axis(side=1,at=c(1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2014,2016))
  axis(2)
  
}
mtext("Year",side=1,line=1,outer=T,adj=.5,cex=1.2)
mtext("Effective Number of Species",outer=T,side=2,line=2.3,cex=1.2)
mtext("Shannon Diversity", side = 3, outer = T, cex = 1.2)

dev.off()

write.csv(DataPlot, file = "~/Dropbox/Thesis/Functional_Diversity/Results with four Subregions/effective_shannon.csv") 


#######################################################################################################
#                                     FUNCTIONAL DENDROGRAM                                           #
#######################################################################################################

#Replace species codes with species names in the relative abundance dataset
sp_sr_ra = relative_abund %>%
  left_join(
    species_list %>%
      select(species_code, spec_common_name),
    by=c('Species_CD'='species_code')) %>%
  as_tibble() %>%
  mutate(
    spec_common_name = tolower(as.character(spec_common_name))) %>%
  arrange(spec_common_name) %>% #alphabeticaly
  select(spec_common_name, Year, Subregion, Relative.Abundance) %>%
  nest(-Year) %>%
  mutate(
    data_wide        = map(
      data, ~ spread(data=.x, spec_common_name, Relative.Abundance) %>%
        # conform to: subregions x species
        select(-Subregion) %>% 
        as.data.frame())) #relative abundance matrix for species*communities

# readr::read_csv() # keeps strings as character, vs evil factors. yay! from readr
# utils::read.csv() # evil! base read of csv
sp_traits = traits %>%
  as_tibble() %>%
  select(-c(1:3)) %>%      # remove extra columns
  arrange(Common_name) %>% # arrange alphabetically by species
  filter(!Common_name == 'Na#') %>%
  mutate( 
    # factors are EVIL!
    Common_name = tolower(as.character(Common_name)),
    # change misclassified categorical traits to numeric 
    Maxlength     = as.numeric(as.character(Maxlength)),
    Trophic_level = as.numeric(as.character(Trophic_level)),
    # ordinal traits
    Water_column   = ordered(
      Water_column, levels=c("Benthic","Demersal","Pelagic site attached","Pelagic non-site attached")),
    Complexity     = ordered(Complexity, levels=c("Low","Medium","High")),
    Gregariousness = ordered(Gregariousness, levels=c("1","2","3"))) %>%
  as.data.frame()

# conform to: species x traits
rownames(sp_traits) = sp_traits$Common_name
sp_traits = sp_traits %>% select(-Common_name)
head(sp_traits)

sp_sr_ra$data_wide[[1]]

#Calculate Gower distances
traits.dist = gowdis(traits, ord="podani")

#Functional Dispersion

#Functional Diversity
sp_sr_ra %>%
  mutate(
    fd = map(data_wide, dbFD)
  )

#functional_div = dbFD(traits.dist, sp_sr_ra$data_wide[[1]])
functional_div = dbFD(sp_traits, sp_sr_ra$data_wide[[1]])

wtf = tibble(
  sp_traits = rownames(sp_traits),
  sp_sr_ra  = colnames(sp_sr_ra$data_wide[[1]]),
  match     = sp_traits == sp_sr_ra) %>%
  filter(!match)




#######################################################################################################
#                               CHARACTERIZING BIVARIATE RELATIONSHIPS                                #
#######################################################################################################

#Create scatterplot matrix of different diversity measures against one another (9" x 9")
lapply(seq_along(alphadiv.list),function(i) {
  pairs.panels(alphadiv.list[[i]][,68:73],#main=names(alphadiv.list)[i],
               labels=c("Richness","Evenness","Simpson","Shannon", "Functional"),
               hist.col="grey80",rug=T,smooth=F,ellipses=F,lm=F,method="spearman") } )

#Test significance of Spearman rank correlations (H0 that rho = 0)
do.call(rbind,lapply(list("Richness","Evenness","species.div","func.div"),function(j) {
  do.call(rbind,lapply(rev(list("richness","evenness","species.div","func.div")),function(k) {
    data.frame(title=paste(j,k,sep="~"),
               correlation=cor.test(alphadiv.list[[2]][,j],alphadiv.list[[2]][,k])$estimate,
               p.value=cor.test(alphadiv.list[[2]][,j],alphadiv.list[[2]][,k])$p.value) } ) ) } ) ) 

#Examine correlations between matrices using Mantel's test
mantel.rtest(taxo.dist,func.dist) #R^2 = 0.74