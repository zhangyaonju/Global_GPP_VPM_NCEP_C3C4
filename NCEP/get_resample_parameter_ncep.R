library(ncdf4)
library(rgdal)
library(raster)
get_resample_parameter_ncep<-function(hv){
  #   hv_file <- read.csv('/data/ifs/ClimateData/NCEP2/resample_para/tile2.csv')
  #   tile_num <- dim(hv_file)[1]
  #   for (k in 92:92){
  #     hv = hv_file[k,2]
  #     hv_sub = substr(hv,2,7)
  #h/v number for each tile
  modisdata.dir<-paste('/data/ifs/modis/products/mod09a1/geotiff/evi/2000/',hv,'/',sep='')
  #modisdata.dir<-paste('/data/ifs/modis/products/mod09a1/geotiff/evi/2000/',noquote(hv_sub),'/',sep='')
  #Gaussian grid for reanalysis data
  ncfile<-nc_open(paste('/data/ifs/ClimateData/NCEP2/Origin/Daily/Ta/air.2m.gauss.2000.nc',sep=''))
  lat<-ncfile$var[[2]]$dim[[2]]$vals
  delta_lat<-array(0,dim=c(93))
  for(i in 1:93){
    delta_lat[i]=lat[i+1]-lat[i]
  }
  
  
  sample.file<-list.files(modisdata.dir)[1]
  sample.data<-readGDAL(paste(modisdata.dir,sample.file,sep=''))
  
  coor.orig<-summary(sample.data)[[5]]
  proj<-array(,dim=c(2400*2400,2))
  
  ####################The centre of the pixel
  for (i in 1:2400){
    for (j in 1:2400){
      proj[(i-1)*2400+j,1]<-coor.orig[1,1]+(j-1)*463.3127
      proj[(i-1)*2400+j,2]<-coor.orig[2,1]+(2400-i)*463.3127 
    }
  }
  ##############obtain the geographical projection
  reproject<-project(proj,summary(sample.data)$proj4string,inv=TRUE)
  ###############################Also need the centre coordinate 
  #################Some problems for the longitude###############################  }
  reproject[,1]<-reproject[,1]/1.875
  j=1
  while(lat[j+1]>reproject[1,2]){
    j=j+1
  }
  for (i in 1:2400){
    if(lat[j+1]>reproject[2400*(i-1)+1,2]){
      j=j+1
    }      
    reproject[2400*(i-1)+1:2400,2]=j+(reproject[2400*(i-1)+1,2]-lat[j])/delta_lat[j]
  }
  
  resample.parameter<-array(,c(5760000,5))
  resa.temp.parameter<-array(,c(4))
  for (i in 1:2400){
    for (j in 1:2400){
      resa.temp.parameter[1]<-cos(pi/2*sqrt((reproject[(i-1)*2400+j,1]-floor(reproject[(i-1)*2400+j,1]))^2+
                                              (reproject[(i-1)*2400+j,2]-floor(reproject[(i-1)*2400+j,2]))^2)/1.414214)^4
      resa.temp.parameter[2]<-cos(pi/2*sqrt((reproject[(i-1)*2400+j,1]-ceiling(reproject[(i-1)*2400+j,1]))^2+
                                              (reproject[(i-1)*2400+j,2]-floor(reproject[(i-1)*2400+j,2]))^2)/1.414214)^4
      resa.temp.parameter[3]<-cos(pi/2*sqrt((reproject[(i-1)*2400+j,1]-ceiling(reproject[(i-1)*2400+j,1]))^2+
                                              (reproject[(i-1)*2400+j,2]-ceiling(reproject[(i-1)*2400+j,2]))^2)/1.414214)^4
      resa.temp.parameter[4]<-cos(pi/2*sqrt((reproject[(i-1)*2400+j,1]-floor(reproject[(i-1)*2400+j,1]))^2+
                                              (reproject[(i-1)*2400+j,2]-ceiling(reproject[(i-1)*2400+j,2]))^2)/1.414214)^4
      resample.parameter[(i-1)*2400+j,1]<-resa.temp.parameter[1]/sum(resa.temp.parameter)
      resample.parameter[(i-1)*2400+j,2]<-resa.temp.parameter[2]/sum(resa.temp.parameter)
      resample.parameter[(i-1)*2400+j,3]<-resa.temp.parameter[3]/sum(resa.temp.parameter)
      resample.parameter[(i-1)*2400+j,4]<-resa.temp.parameter[4]/sum(resa.temp.parameter)
      resample.parameter[(i-1)*2400+j,5]<-(floor(reproject[(i-1)*2400+j,1])+96)+192*(floor(reproject[(i-1)*2400+j,2]))
    }
  }
  
  parameter.out.file<-paste('/data/ifs/ClimateData/NCEP2/resample_para/',noquote(hv),'.tif',sep='')
  
  sample.data@data[[1]]<-as.vector(resample.parameter[,1])
  sample.data@data[[2]]<-as.vector(resample.parameter[,2])
  sample.data@data[[3]]<-as.vector(resample.parameter[,3])
  sample.data@data[[4]]<-as.vector(resample.parameter[,4])
  sample.data@data[[5]]<-as.vector(resample.parameter[,5])
  
  writeGDAL(sample.data,parameter.out.file,drivername='GTiff',type='Float32')
}
