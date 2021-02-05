
setwd("C:/Users/etelford.IC.003/Dropbox/Phd/Corona_chapter_1X/R/SDM_workshop")


#1. Use the getData function to download the Worldclim climate data set at 2.5 arc-minute resolution.
#You want the bioclimatic variables (use var="bio", see ?getData in the dismo package).

library(dismo)
library(raster)
bioRasts <- getData(name="worldclim", #other options available 
                       res = 2.5, # resolution
                       var = "bio") # which variable(s)?

 

#2. Make a raster stack of the bio10, bio11, bio18, and bio19 bioclimatic variables and clip this raster
#stack to the outline of Australia (not the extent) using the shapefile provided (~data/oz_nz_aea.shp).
#(Note that the shapefile also contains New Zealand, so you will have to do something about that before
#you perform the clipping operation, among other things. . . ).

bioRasts <- bioRasts[[c(10, 11,18,19)]]

aus_nz <- shapefile("C:/Users/etelford.IC.003/Dropbox/Phd/Corona_chapter_1X/R/SDM_workshop/oz_nz_aea.shp")
aus <- subset(aus_nz, NAME=="Australia")
plot(aus)
projection(aus)
projection(bioRasts)
#they are not the same projection

# project the vector data
ausProj <-spTransform(aus, CRSobj =projection(bioRasts))
# create a mask raster to use for clipping
ausRast <-raster(ausProj, res=res(bioRasts))
ausRast <-setValues(ausRast, 1)
ausRast <-mask(ausRast, ausProj)

origin(ausRast)
origin(bioRasts)
origin(ausRast) <-origin(bioRasts)

# now we clip by multiplying the mask and the bioRasts together
bioRasts <- bioRasts*ausRast
names(bioRasts) <-c("bio10", "bio11", "bio18", "bio19")
plot(bioRasts)


#3. Use the gbif function to download records for the 
#Austral grass tree(Xanthorrhoea australis).Clean up the resulting
#data frame by removing records without geographic coordinates 
#those that fall outside the Australian mainland.
#Convert the data to aSpatialPointsDataFramewith  the  correct  CRS
#and  containing  only these  attributes:acceptedScientificName,
#institutionCode,lon,lat, andyear. Save yourSpatialPointsDataFrameas 
#a shapefile.

grassTree <-gbif(genus="Xanthorrhoea", species = "australis", geo=T)
# select columns
grassTree <- grassTree[,c("acceptedScientificName", "institutionCode", "lon","lat",  "year")]
# remove all NA's
grassTree <-na.omit(grassTree)

# convert to a spatial object
coordinates(grassTree) <-c("lon", "lat")
projection(grassTree) <-projection(ausProj)

ausGrassTree <-over(grassTree, ausProj)# overlay
ausGrassTree <-which(is.na(ausGrassTree))# get records that are NA
grassTree <- grassTree[-ausGrassTree,] # remove NA records

#get the cell number of each record
cells <-cellFromXY(bioRasts$bio10, grassTree)
# use the duplicated function to find duplicate cell numbers
dups <-duplicated(cells)
# returns'TRUE'if the cell number is duplicated
grassTreeND <- grassTree[!dups,]# the ! = not, so here we keep dups == FALSE
shapefile(grassTreeND, "grassTree.shp", overwrite=T)# save as shapefil

#4. Make a simple map of the cleaned species occurrence records from GBIF, using a color ramp or
#symbolization scheme to indicate the year the record was collected. Make sure to include the polygon
#of Australia and bio10 as the background. NOTE: All data in this map should be in the original
#projection of the Australia & New Zealand shapefile (i.e., not WGS84). Save the transformed bio10
#raster as a GeoTiff.

# transform the grass tree data from WGS84
grassTreeProj <-spTransform(grassTreeND, CRSobj =projection(aus))
bio10 <- bioRasts$bio10
bio10Proj <-projectRaster(bio10, crs=projection(aus)) # transform

writeRaster(bio10Proj, "bio10Proj.tif", overwrite=T)# save raster as GeoTif

library(colorRamps)
plot(bio10Proj/10, col=rgb.tables(1000), alpha=0.5)
# alpha  sets transparencyplot(aus, add=T)
# add the polygon
# here's where I setup the color scheme# this part: [nrow(grassTreeProj):1] reverses the order of the color ramp

yearCol <-topo.colors(n=nrow(grassTreeProj))[nrow(grassTreeProj):1]# add the points with the color order set to the year columnpoints(grassTreeProj, pch=21, lwd=1, cex=0.8,bg=yearCol[order(grassTreeProj$year)])

points(grassTreeProj, pch=21, lwd=1, cex=0.8,bg=yearCol[order(grassTreeProj$year)])

#5. Use the point data to extract the bioclimatic variables from the raster stack and compare the climate
#conditions where this species has been observed to the broader climate of Australia. A few hints: Have
#a look at the sampleRandom function. To perform the comparison between climates where the species
#is present and Australia more broadly, you have a number of options. You might try scatter plots, box
#plots or histograms, but you do not need to do any statistical analyses (in other words, see what you
#can learn from simple plots alone). Answer the question: How does the climate where X. australis
#has been observed differ from that of Australian climates more generally?

sppDat <-data.frame(extract(bioRasts, grassTreeND))
sppDat <-na.omit(sppDat)

bgDat <-data.frame(sampleRandom(bioRasts, size=10000))

par(mfrow=c(2,2))

# bio10
plot(density(bgDat$bio10/10), col="black", ylim=c(0,0.24),
     main="bio10: Summer Temperature",
     xlab="Summer Temperature (C)"),
lines(density(sppDat$bio10/10), col="red")

# bio11
plot(density(bgDat$bio11/10), col="black", ylim=c(0,0.24),
     main="bio11: Winter Temperature",
     xlab="Winter Temperature (C)")
lines(density(sppDat$bio11/10), col="red")
legend(x=13, y=0.2, legend=c("Grass tree", "Random"),col=c("red", "black"), lwd=2, cex=0.8, bty="n")

# bio18
plot(density(bgDat$bio18), col="black", ylim=c(0,0.01),
     main="bio18: Summer Precipitation",
     xlab="Summer Precipitation (mm)")
lines(density(sppDat$bio18), col="blue")

# bio19
plot(density(bgDat$bio19), col="black", ylim=c(0,0.015),
     main="bio19: Winter Precipitation",
     xlab="Winter Precipitation (mm)")
lines(density(sppDat$bio19), col="blue")
legend(x=500, y=0.01, legend=c("Grass tree", "Random"),
       col=c("blue", "black"), lwd=2, cex=0.8, bty="n")

#6. Create a raster of the number of species observations in each grid cell.

ozRast <- bioRasts$bio10
ozRast[which(!is.na(ozRast[]))] <- 0
# now we get the cell numbers for each downloaded record
cells <-cellFromXY(ozRast, grassTree)

# count the number of times each cell number occurs using the'table'function
counts <-data.frame(table(cells))

# need to some work here to get the cell numbers to be numeric because they
# appear as names in the counts object and so are character strings
cellNum <-as.numeric(as.character(counts$cells))

# we index the raster by cellNum and assign the counts (counts$Freq) as follows:
ozRast[cellNum] <- counts$Freq

# the cell sizes are small, so here I plot a small area to show the pattern
plot(ozRast, xlim=c(145,152), ylim=c(-37,-34), 
     col=rgb.tables(1000),maxpixels=ncell(ozRast)
