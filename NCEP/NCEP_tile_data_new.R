#VPM global simulation 
#1st stage: data preparation
#4th chapter: 
#developed by Yao Zhang
#Apr 28, 2015
library(raster)
library(rgdal)
library(ncdf4)

NCEP_tile_data_new<-function(year,hv){
  #####=====================================================================
  ## we need to read in the parameter for resampling and apply the algorithm
  #####
  #we first search whether it is western hemisphere or eastern hemisphere
  #  0~180 -180~0
  # ew: 0 for western hemisphere, 1 for eastern hesmisphere
  ew<-as.numeric(substr(hv,2,3))>17
  
  parameter.file<-paste('/data/ifs/ClimateData/NCEP2/resample_para/',hv,'.tif',sep='')
  
  input.8day.file<-c(paste('/data/ifs/ClimateData/NCEP2/Origin/Daily/TaDay/ta.daytime.2m.8day.',year,'.nc',sep=''),
                     paste('/data/ifs/ClimateData/NCEP2/Origin/Daily/DSSW/rad.8day.',year,'.nc',sep=''))
  
  input.temp.data<-nc_open(input.8day.file[1])
  temp.8day.data<-ncvar_get(input.temp.data)
  input.radi.data<-nc_open(input.8day.file[2])
  radiation.8day.data<-ncvar_get(input.radi.data)
  
  
  resample.parameter.data<-readGDAL(parameter.file)
  coor.orig<-summary(resample.parameter.data)[[2]]
  resample.parameter<-resample.parameter.data@data
  coor.x<-array(,dim=c(2400))
  coor.y<-array(,dim=c(2400))
  coor.x<-coor.orig[1,1]+0:2399*463.3127
  coor.y<-coor.orig[2,2]-0:2399*463.3127 
  
  ###recalculate the parameter because of the difference between the sinuside projection and gaussian grid
  if (ew > 0){
    resample.parameter[,5]<-resample.parameter[,5] - 96
  }else{
    resample.parameter[,5]<-resample.parameter[,5] + 96
  }
  ###add by Yao Zhang on 20170420
  resample.parameter[,5][resample.parameter[,5]<0]=0
  
  dir.create(paste('/data/ifs/VPM/driving_data/NCEP2/Tile/DSSW/',year,'/',sep=""),recursive = T)  
  dir.create(paste('/data/ifs/VPM/driving_data/NCEP2/Tile/TaDay/',year,'/',sep=""),recursive = T)  
  radiation.fine.output<-paste('/data/ifs/VPM/driving_data/NCEP2/Tile/DSSW/',year,'/rad.',hv,'.',year,'.nc',sep='')
  temp.fine.output<-paste('/data/ifs/VPM/driving_data/NCEP2/Tile/TaDay/',year,'/temp.',hv,'.',year,'.nc',sep='')
  
  fine.air.temp.data<-array(,dim=c(2400,2400,46))
  fine.radiation.data<-array(,dim=c(2400,2400,46))
  
  for (scene.num in 1:46){
    fine.air.temp.data[,,scene.num]<-round((temp.8day.data[,,scene.num][resample.parameter[,5]+1]*resample.parameter[,1]+
                                              temp.8day.data[,,scene.num][resample.parameter[,5]+2]*resample.parameter[,2]+
                                              temp.8day.data[,,scene.num][resample.parameter[,5]+194]*resample.parameter[,3]+
                                              temp.8day.data[,,scene.num][resample.parameter[,5]+193]*resample.parameter[,4])*100)
    fine.radiation.data[,,scene.num]<-round((radiation.8day.data[,,scene.num][resample.parameter[,5]+1]*resample.parameter[,1]+
                                              radiation.8day.data[,,scene.num][resample.parameter[,5]+2]*resample.parameter[,2]+
                                              radiation.8day.data[,,scene.num][resample.parameter[,5]+194]*resample.parameter[,3]+
                                              radiation.8day.data[,,scene.num][resample.parameter[,5]+193]*resample.parameter[,4])*100)
  }
  fine.radiation.data[is.na(fine.radiation.data)]<- -32767
  fine.air.temp.data[is.na(fine.air.temp.data)]<- -32767
  
  #=================================create the ncdf file
  # define the dimension
  x<-ncdim_def('sin_x','meter',coor.x)
  y<-ncdim_def('sin_y','meter',coor.y)
  t<-ncdim_def('Time',paste('days since ',year,'-01-01',sep=''),(0:45)*8, unlim=TRUE)
  
  #define variables
  fillvalue<- -32767
  dlname<-'Photosynthetical active radiation'
  PAR<-ncvar_def('PAR','0.01 PPFD',list(x,y,t),fillvalue,dlname,prec='integer')
  
  #create netCDF file and put arrays
  ncout<-nc_create(radiation.fine.output,PAR)
  
  #put variables
  ncvar_put(ncout,PAR,fine.radiation.data)
  nc_close(ncout)
  
  ####=================================================================
  
  #define variables
  fillvalue<- -32767
  dlname<-'daytime mean temperature'
  temperature<-ncvar_def('Temp','0.01 degreeC',list(x,y,t),fillvalue,dlname,prec='integer')
  
  #create netCDF file and put arrays
  ncout<-nc_create(temp.fine.output,temperature)
  
  #put variables
  ncvar_put(ncout,temperature,fine.air.temp.data)
  nc_close(ncout) 
}

