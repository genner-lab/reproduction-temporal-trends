#!/usr/bin/env Rscript

# good post on maps
#https://hansenjohnson.org/post/bathymetric-maps-in-r/

# load libs and data
library("marmap")
library("oce")
library("ocedata")
data("coastlineWorldFine")

# set box
left <- -20
right <- 10
bottom <- 40
top <- 60

# get bathy data
base <- getNOAA.bathy(lon1=left,lon2=right,lat1=bottom,lat2=top,resolution=1)

# convert bathymetry
bathyLon <- as.numeric(rownames(base))
bathyLat <- as.numeric(colnames(base))
bathyZ <- as.numeric(base)
dim(bathyZ) <- dim(base)


# plot UK and Europe map
pdf(file="temp/results/figures/uk-map-temp.pdf",width=8,height=8,useDingbats=FALSE)
    plot(coastlineWorldFine,clon=-5,clat=50,span=2000,projection="+proj=merc",col="lightgrey",lonlabels=FALSE,latlabels=FALSE)
    mapPolygon(longitude=c(-6.5,-6.5,-3.5,-3.5),latitude=c(49,51,51,49),lty=1,lwd=2,border="#CC3311")
dev.off()


# plot study site
pdf(file="temp/results/figures/bathymap-temp.pdf",width=8,height=8,useDingbats=FALSE)
    plot(coastlineWorldFine,clon=-5,clat=50,span=200,projection="+proj=merc",col="lightgrey",grid=TRUE)
    # plot bathymetry
    mapContour(bathyLon, bathyLat, bathyZ,
        levels=c(-25,-50,-75,-100),
        lwd=c(1,1,1,2),
        lty=c(3,2,1,1),
        col='darkgray')
    # add depth legend
    legend("bottomright",seg.len=3,cex=0.7,
        lwd=c(1,1,1,2),
        lty=c(3,2,1,1),
        legend=c("25","50","75","100"),
        col='darkgray',title="Depth [m]",bg="white")
    # add scalebar
    mapScalebar(x="topright",length=20,cex=0.7)
    # add map data
    mapPoints(longitude=-4.136271,latitude=50.367561,pch=15,col="black",cex=2) # Plymouth
    mapPoints(longitude=-5.533078,latitude=50.116452,pch=15,col="black",cex=1.5) # Penzance
    mapPoints(longitude=-4.22,latitude=50.25,pch=17,col="#CC3311",cex=1.5) # L4
    mapPoints(longitude=-4.3,latitude=50.18,pch=17,col="#CC3311",cex=1.5) # L5
    mapPoints(longitude=-4.37,latitude=50.03,pch=17,col="#CC3311",cex=1.5) # E1
dev.off()
