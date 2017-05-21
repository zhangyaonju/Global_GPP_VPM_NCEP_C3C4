climate_data_ncep_aggregate<-function(year){

  air.temp.file<-nc_open(paste('/data/eomf/users/yzhang/NCEP/origin/air.2m.gauss.',year,'.nc',sep=''))
  radiation.file<-nc_open(paste('/data/eomf/users/yzhang/NCEP/origin/dswrf.sfc.gauss.',year,'.nc',sep=''))
  outfile<-c(paste('/data/eomf/users/yzhang/NCEP/temperature/air.2m.8day.',year,'.nc',sep=''),
             paste('/data/eomf/users/yzhang/NCEP/radiation/rad.8day.',year,'.nc',sep=''))
  
  air.temp.raw.data<-ncvar_get(air.temp.file,'air')
  air.temp.day.data<-array(0,dim=c(192,94,8))
  air.temp.8day.data<-array(0,dim=c(192,94,46))
  
  timezone.offset<--round((air.temp.file$var[[1]]$dim[[1]]$vals*(air.temp.file$var[[1]]$dim[[1]]$vals<180)-
                              (360-air.temp.file$var[[1]]$dim[[1]]$vals)*(air.temp.file$var[[1]]$dim[[1]]$vals>179))/15/6)
  
  obs.num<-dim(air.temp.raw.data)[3]
  
  
  for(x in 1:192){
    for (y in 1:94){
      for (i in 2:45){
        for (j in 1:8){
          air.temp.day.data[x,y,j]<-min(air.temp.raw.data[x,y,(timezone.offset[x]+(i-1)*8*4+(j-1)*4+1):(timezone.offset[x]+(i-1)*8*4+4*j)])*0.25+
            max(air.temp.raw.data[x,y,(timezone.offset[x]+(i-1)*8*4+(j-1)*4+1):(timezone.offset[x]+(i-1)*8*4+4*j)])*0.75
        }
        air.temp.8day.data[x,y,i]<-mean(air.temp.day.data[x,y,1:min(8,obs.num/4-i*8)])
      }
    } 
  }
  
  
  for(x in 1:192){
    for (y in 1:94){
      for(j in 1:8){
        air.temp.day.data[x,y,j]<-min(air.temp.raw.data[x,y,max(1,(timezone.offset[x]+(j-1)*4+1)):(timezone.offset[x]+j*4)])*0.25+
          max(air.temp.raw.data[x,y,max(1,(timezone.offset[x]+(j-1)*4+1)):(timezone.offset[x]+j*4)])*0.75
      }
      air.temp.8day.data[x,y,1]<-mean(air.temp.day.data[x,y,1:8])
    }
  }
  
  for(x in 1:192){
    for (y in 1:94){
      for(j in 1:(obs.num/4-360)){
        air.temp.day.data[x,y,j]<-min(air.temp.raw.data[x,y,(timezone.offset[x]+1440+(j-1)*4+1):min(obs.num,(timezone.offset[x]+1440+j*4))])*0.25+
          max(air.temp.raw.data[x,y,(timezone.offset[x]+1440+(j-1)*4+1):min(obs.num,(timezone.offset[x]+1440+j*4))])*0.75
      }
      air.temp.8day.data[x,y,46]<-mean(air.temp.day.data[x,y,1:(obs.num/4-360)])
    }
  }
  
  time.output<-air.temp.file$var[[1]]$dim[[3]]$vals[(0:45)*32+1]
  
  x<-ncdim_def('Lon','degreesE',air.temp.file$var[[1]]$dim[[1]]$vals)
  y<-ncdim_def('Lat','degreesN',air.temp.file$var[[1]]$dim[[2]]$vals)
  t<-ncdim_def("Time",'hours since 1800-01-01 00:00:0.0',time.output,unlim=TRUE)
  
  ncair<-ncvar_def('air',"degK",list(x,y,t),missval=-9.96920996838687e+36,longname='8-day 2m air temperature',prec='float')
  
  ncair.out<-nc_create(outfile[1],ncair)
  ncvar_put(ncair.out,'air',air.temp.8day.data)
  
  nc_close(ncair.out)
  
  
  radiation.raw.data<-ncvar_get(radiation.file,'dswrf')
  radiation.8day.data<-array(0,dim=c(192,94,46))

  
  obs.num<-dim(radiation.raw.data)[3]
  
  
  for(x in 1:192){
    for (y in 1:94){
      for (i in 1:45){
        radiation.8day.data[x,y,i]<-mean(radiation.raw.data[x,y,(8*i-7):(i*8)])
      }
    } 
  }


  
  for(x in 1:192){
    for (y in 1:94){
      radiation.8day.data[x,y,46]<-mean(radiation.raw.data[x,y,361:obs.num])
    }
  }
  
  time.output<-radiation.file$var[[1]]$dim[[3]]$vals[(0:45)*8+1]
  
  x<-ncdim_def('Lon','degreesE',radiation.file$var[[1]]$dim[[1]]$vals)
  y<-ncdim_def('Lat','degreesN',radiation.file$var[[1]]$dim[[2]]$vals)
  t<-ncdim_def("Time",'hours since 1800-01-01 00:00:0.0',time.output,unlim=TRUE)
  
  ncrad<-ncvar_def('dswrf',"W/m^2",list(x,y,t),missval=-9.96920996838687e+36,longname='8-day Downward Solar Radiation Flux at surface',prec='float')
  
  ncrad.out<-nc_create(outfile[2],ncrad)
  ncvar_put(ncrad.out,'dswrf',radiation.8day.data)
  
  nc_close(ncrad.out)
  
  gc()
}


