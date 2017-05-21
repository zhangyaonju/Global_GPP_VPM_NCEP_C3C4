#######
##this code is used to estimate the C4 fraction for crop
##before using this code, a global cropland fraction 
##at the target 0.0833333 degree is needed.

##first read in the 0.5 deg
pct2008cro<-raster("V:/users/yzhang/DATA/mcd12c1/MCD12C1.A2008001.051.2012283145256.IGBP.pct.tif",band=13)
pct2008cnv<-raster("V:/users/yzhang/DATA/mcd12c1/MCD12C1.A2008001.051.2012283145256.IGBP.pct.tif",band=15)

cro<-(pct2008cro+pct2008cnv)/100
writeRaster(cro,"V:/VPM/driving_data/C3C4/crop/crop.tif")
#####resample this layer to 0.083333
c4crop<-raster("V:/VPM/driving_data/C3C4/crop/c4fraction.tif")
rescrop<-raster("V:/VPM/driving_data/C3C4/crop/resampledcrop.tif")
cropc4<-c4crop/rescrop
cropc4[cropc4>1]<-1
cropc4[rescrop==0]<-0
writeRaster(cropc4,"V:/VPM/driving_data/C3C4/cropc4_pct.tif")
######resample nature vegetation to 0.083333
nat<-raster("V:/VPM/driving_data/C3C4/nature_vegetation/naturec4.tif")
out<-nat/100
out[out<0]=0
writeRaster(out,"V:/VPM/driving_data/C3C4/natc4_pct.tif")

###run this for the entire 296 tiles

library(raster)
library(rgdal)
lc<-list.files("/data/ifs/VPM/driving_data/MCD12Q1/2004/",pattern="*.tif",full.names=T)
tile<-substr(basename(lc),18,23)
for (i in 1:length(lc)){
  lcdata<-raster(lc[i])
  ext<-extent(lcdata)
  system(paste("gdalwarp -t_srs '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs' -tr 463.3127 463.3127 -te ",
               ext[1]," ",ext[3]," ",ext[2]," ",ext[4]," /data/ifs/VPM/driving_data/C3C4/prepare/cropc4_pct.tif /data/ifs/VPM/driving_data/C3C4/crop/",
               tile[i],".crop.tif",sep=""))
  system(paste("gdalwarp -t_srs '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs' -tr 463.3127 463.3127 -te ",
               ext[1]," ",ext[3]," ",ext[2]," ",ext[4]," /data/ifs/VPM/driving_data/C3C4/prepare/natc4_pct.tif /data/ifs/VPM/driving_data/C3C4/nature/",
               tile[i],".nature.tif",sep=""))
  
}
