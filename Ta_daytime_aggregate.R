library(ncdf4)
library(rgdal)
library(raster)
Ta_daytime_aggregate<-function(year){
  ## read data
  ## air.temp.file<-nc_open(paste('D:/research/NCEP2/Daily/Ta/air.2m.gauss.',year,'.nc',sep=''))
#  air.temp.file<-nc_open(paste('/data/ifs/ClimateData/NCEP2/Origin/Daily/Ta/air.2m.gauss.',year,'.nc',sep=''))
  tmax.temp.file<-nc_open(paste('/data/ifs/ClimateData/NCEP2/Origin/Daily/Tamax/tmax.2m.gauss.',year,'.nc',sep=''))
  tmin.temp.file<-nc_open(paste('/data/ifs/ClimateData/NCEP2/Origin/Daily/Tamin/tmin.2m.gauss.',year,'.nc',sep=''))
#  radiation.file<-nc_open(paste('/data/ifs/ClimateData/NCEP2/Origin/Daily/DSSW/dswrf.sfc.gauss.',year,'.nc',sep=''))
  outfile<-c(paste('/data/ifs/ClimateData/NCEP2/Origin/Daily/TaDay/ta.daytime.2m.8day.',year,'.nc',sep=''))
             ###you need to know the variables and dimension of the nc file in PanoplyWin,wxc
  #  air.temp.raw.data<-ncvar_get(air.temp.file,'air')
 # air.temp.raw.data<-ncvar_get(air.temp.file,'air')
  tmax.temp.raw.data<-ncvar_get(tmax.temp.file,'tmax')
  tmin.temp.raw.data<-ncvar_get(tmin.temp.file,'tmin')
  # air.temp.day.data<-array(0,dim=c(192,94,8))    ###save every 8-day data
  # air.temp.8day.data<-array(0,dim=c(192,94,46))  ###all the 8-day data in a year,wxc
 # air.temp.8day.data<-array(0,dim=c(192,94,46))
  ta_daytime.temp.8day.data<-array(0,dim=c(192,94,46))
  #####timezone.offset is to calculate the daily average from 6-hourly data, I delete it because we use daily data now.wxc 
  # wxc timezone.offset<--round((air.temp.file$var[[1]]$dim[[1]]$vals*(air.temp.file$var[[1]]$dim[[1]]$vals<180)-
  # wxc                            (360-air.temp.file$var[[1]]$dim[[1]]$vals)*(air.temp.file$var[[1]]$dim[[1]]$vals>179))/15/6)
  
  #############################################Average temperature in daytime####################################
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  #########wxc,num of 8-day###############
  obs.num<-dim(tmax.temp.raw.data)[3]
  eightday.num<-floor(obs.num/8);
  # delete by wxc  
  #  for(x in 1:192){
  #    for (y in 1:94){
  #      for (i in 2:45){
  #        for (j in 1:8){
  #          air.temp.day.data[x,y,j]<-min(air.temp.raw.data[x,y,(timezone.offset[x]+(i-1)*8*4+(j-1)*4+1):(timezone.offset[x]+(i-1)*8*4+4*j)])*0.25+
  #            max(air.temp.raw.data[x,y,(timezone.offset[x]+(i-1)*8*4+(j-1)*4+1):(timezone.offset[x]+(i-1)*8*4+4*j)])*0.75
  #        }
  #        air.temp.8day.data[x,y,i]<-mean(air.temp.day.data[x,y,1:min(8,obs.num/4-i*8)])
  #      }
  #    } 
  #  }
  ############################################Modified for daily data,wxc####################################
  for (x in 1:192){
    for (y in 1:94){
      for (i in 1:eightday.num){
        #        for (j in 1:8){
        #          air.temp.day.data[x,y,j]<-air.temp.raw.data[x,y,(i-1)*8+j]
        #        }
        ta_daytime.temp.8day.data[x,y,i]<-mean(tmax.temp.raw.data[x,y,((i-1)*8+1):(i*8)]*0.75
                                              +tmin.temp.raw.data[x,y,((i-1)*8+1):(i*8)]*0.25)
      }
    }
  }
  
  #  for(x in 1:192){
  #    for (y in 1:94){
  #      for(j in 1:8){
  #        air.temp.day.data[x,y,j]<-min(air.temp.raw.data[x,y,max(1,(timezone.offset[x]+(j-1)*4+1)):(timezone.offset[x]+j*4)])*0.25+
  #          max(air.temp.raw.data[x,y,max(1,(timezone.offset[x]+(j-1)*4+1)):(timezone.offset[x]+j*4)])*0.75
  #      }
  #      air.temp.8day.data[x,y,1]<-mean(air.temp.day.data[x,y,1:8])
  #    }
  #  }
  
  #  for(x in 1:192){
  #    for (y in 1:94){
  #      for(j in 1:(obs.num-360)){
  #        air.temp.day.data[x,y,j]<-min(air.temp.raw.data[x,y,(timezone.offset[x]+1440+(j-1)*4+1):min(obs.num,(timezone.offset[x]+1440+j*4))])*0.25+
  #          max(air.temp.raw.data[x,y,(timezone.offset[x]+1440+(j-1)*4+1):min(obs.num,(timezone.offset[x]+1440+j*4))])*0.75
  #      }
  #      air.temp.8day.data[x,y,46]<-mean(air.temp.day.data[x,y,1:(obs.num/4-360)])
  #    }
  #  }
  ##########################################Modified for daily data,wxc####################################
  ########################if else sentence is for incomplete data in some years############################
  if (eightday.num==45){
    for (x in 1:192){
      for (y in 1:94){
        #          for (j in 1:(obs.num-360)){
        #            air.temp.day.data[x,y,j]<-air.temp.raw.data[x,y,1:(obs.num-360)]
        #          }
        ta_daytime.temp.8day.data[x,y,46]<-mean(tmax.temp.raw.data[x,y,361:obs.num]*0.75
                                               +tmin.temp.raw.data[x,y,361:obs.num]*0.25)
      }
    }
  }else{
    for (x in 1:192){
      for (y in 1:94){
        for (i in eightday.num+1:45){
          #          for (j in 1:8){
          #            air.temp.day.data[x,y,j]<-air.temp.raw.data[x,y,(i-1)*8+j]            
          #          }
          ta_daytime.temp.8day.data[x,y,i]<-mean(tmax.temp.raw.data[x,y,((i-1)*8+1):(i*8)]*0.75
                                                +tmin.temp.raw.data[x,y,((i-1)*8+1):(i*8)]*0.25)
        }
      }
    }
    for (x in 1:192){
      for (y in 1:94){
        #        for (j in 1:(obs.num-360)){
        #          air.temp.day.data[x,y,j]<-air.temp.raw.data[x,y,1:(obs.num-360)]
        #        }
        ta_daytime.temp.8day.data[x,y,46]<-mean(tmax.temp.raw.data[x,y,361:obs.num]*0.75
                                               +tmin.temp.raw.data[x,y,361:obs.num]*0.25)
      }
    }
  }
  
  ##################################################################################################################
  
  time.output<-tmax.temp.file$var[[2]]$dim[[4]]$vals[(0:45)*8+1]
  
  x<-ncdim_def('Lon','degreesE',tmax.temp.file$var[[2]]$dim[[1]]$vals)
  y<-ncdim_def('Lat','degreesN',tmax.temp.file$var[[2]]$dim[[2]]$vals)
  t<-ncdim_def("Time",'hours since 1800-01-01 00:00:0.0',time.output,unlim=TRUE)
  
  ncair<-ncvar_def('Ta_daytime',"degK",list(x,y,t),missval=-9.96920996838687e+36,longname='8-day 2m mean temperature at the daytime',prec='float')
  
  ncair.out<-nc_create(outfile[1],ncair)
  ncvar_put(ncair.out,'Ta_daytime',ta_daytime.temp.8day.data)
  
  nc_close(ncair.out)
  
  ##################################################Tmin#########################################################
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  #########wxc,num of 8-day###############
 
  gc()
  
}


