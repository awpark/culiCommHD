## ----setup, include=FALSE----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----loadLibs----------------------------------------------------------------------------------------
library(dplyr)
library(maptools)
library(raster)
library(metacom)


## ----loadData----------------------------------------------------------------------------------------
load("get_db.Rda")

#check for UTM typos
Vial.ID.issues <- db$Vial.ID[which(db$UTM.1==db$UTM.2)]

#fix a few database typos confirmed by email from Stacey Vigil to Andrew Park (11/11/2022)
#AP's copy has UTM.2=UTM.1, Stacey looked up correct UTM.2 in original data sheets
db$UTM.2[which(db$Vial.ID=="FL0276")] <- 2957038
db$UTM.2[which(db$Vial.ID=="FL0315")] <- 2957057
db$UTM.2[which(db$Vial.ID=="FL0567")] <- 3001295
db$UTM.2[which(db$Vial.ID=="FL0568")] <- 3002091
db$UTM.2[which(db$Vial.ID=="FL0816")] <- 3034831
db$UTM.2[which(db$Vial.ID=="FL0817")] <- 3034831
db$UTM.2[which(db$Vial.ID=="FL0818")] <- 3034831


#here I [AP] just convert UTM lat/long using http://rcn.montana.edu/resources/Converter.aspx
#FL is un UTM zone 17 (northern hemisphere)
#vial FL0276
db$Lat[which(db$Vial.ID=="FL0276")] <- 26.73170
db$Long[which(db$Vial.ID=="FL0276")] <- -81.90239
#vial FL0315
db$Lat[which(db$Vial.ID=="FL0315")] <- 26.73190
db$Long[which(db$Vial.ID=="FL0315")] <- -81.89777
#vial FL0567
db$Lat[which(db$Vial.ID=="FL0567")] <- 27.12999
db$Long[which(db$Vial.ID=="FL0567")] <- -82.08142
#vial FL0568
db$Lat[which(db$Vial.ID=="FL0568")] <- 27.13710
db$Long[which(db$Vial.ID=="FL0568")] <- -82.09118
#vial FL0816
db$Lat[which(db$Vial.ID=="FL0816")] <- 27.43069
db$Long[which(db$Vial.ID=="FL0816")] <- -82.31852
#vial FL0817
db$Lat[which(db$Vial.ID=="FL0817")] <- 27.43069
db$Long[which(db$Vial.ID=="FL0817")] <- -82.31852
#vial FL0818
db$Lat[which(db$Vial.ID=="FL0818")] <- 27.43069
db$Long[which(db$Vial.ID=="FL0818")] <- -82.31852

db %<>% dplyr::filter(Lat<35) # four sites are up in Michigan and spatially disconnected from this metacommunity analysis



## ----wmaLatLong--------------------------------------------------------------------------------------

centroids <- db |> 
  group_by(Site) |> 
  summarize(cLat=mean(Lat),cLon=mean(Long))


save(centroids,file="get_siteCentroids.Rda")



## ----------------------------------------------------------------------------------------------------

sites <- unique(db$Site)
sub.db <- db[,11:62]

#for loop that outputs interaction matrix combining all years
ret <- matrix(0, ncol= 52, nrow=length(sites)) 
latitude=vector(); longitude=vector()
for(i in 1:length(sites)){
  temp=sub.db[which(db$Site==sites[i]),]	
  ret[i,] = colSums(temp)
  latitude[i] = mean(db[which(db$Site == sites[i]),'Lat'])
  longitude[i] = mean(db[which(db$Site == sites[i]),'Long'])
}
rownames(ret)=sites; colnames(ret)=colnames(sub.db)

#ret.lat=ret[order(latitude),]
ret.nozero <- ret[-which(rowSums(ret) == 0),]
ret.nozero <- ret.nozero > 0 
sites.nozero <- sites[-which(rowSums(ret)==0)]
latitude <- latitude[-which(rowSums(ret)==0)]
longitude <- longitude[-which(rowSums(ret)==0)]



## ----------------------------------------------------------------------------------------------------
#pull out site scores from reciprocal averaging
scores <- OrderMatrix(ret.nozero, scores=1, outputScores=TRUE)

#Coherence 
coh.ret.nozero <- Coherence(ret.nozero > 0, order=TRUE, orderNulls=TRUE)
#Species Turnover
turn.ret.nozero <- Turnover(ret.nozero > 0, order=TRUE, orderNulls=TRUE)
#Boundary clumping
clump.ret.nozero <- BoundaryClump(ret.nozero > 0, order=TRUE)



## ----------------------------------------------------------------------------------------------------
ret.spatial <- ret.nozero[order(latitude),order(scores$speciesscores)]

#Coherence 
coh.ret.spatial <- Coherence(ret.spatial > 0, order=FALSE, sims=100, orderNulls=TRUE)
#Species Turnover
turn.ret.spatial <- Turnover(ret.spatial > 0, order=FALSE, orderNulls=TRUE)
#Boundary clumping
clump.ret.spatial <- BoundaryClump(ret.spatial > 0, order=FALSE)


