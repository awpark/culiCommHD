## ----setup, include=FALSE------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ------------------------------------------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(raster)


## ----loadData------------------------------------------------------------------------------------------
#morbidity mortality data for HD in deer from SCWDS by county 2008-2012 (coinciding light trap data)
load("get_hdData.Rda")


## ----smooth--------------------------------------------------------------------------------------------
library(MASS)
dens<-kde2d(hdData$X_centroid,hdData$Y_centroid,h=c(5.0,5.0),n=200,lims=c(min(hdData$X_centroid)-5.0,max(hdData$X_centroid)+5.0,min(hdData$Y_centroid)-5.0,max(hdData$Y_centroid)+5.0))
dens$z<-dens$z[,dim(dens$z)[2]:1]
projCRS=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
r=raster::raster(t(dens$z),xmn=min(dens$x),xmx=max(dens$x),ymn=min(dens$y),ymx=max(dens$y),crs=projCRS)

