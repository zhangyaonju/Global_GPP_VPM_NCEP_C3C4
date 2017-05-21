#-------------------------------------------------------------------------------
# Name:        Global GPP simulation using the VPM model
# Inputs: LUE, Topt, Tmin, Tmax, LSWImax, 8-day multi-EVI, 8-day multi-LSWI, 8-day multi-NCEP2(Rs and 2-m Temp)
#
# Author:      Yao Zhang
#
# Created:     12/01/2015
# Modified:    05/04/2016
# Modified for C3C4:
#              04/28/2017
# Copyright:   (c) eomf 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import multiprocessing
import os
from os import listdir
from os.path import isfile, join
from osgeo import gdal
from osgeo.gdalconst import *
import numpy
import numpy.ma as ma
from netCDF4 import Dataset
import time
#import h5py

startTime=time.time()

# input the parent directory
root='/data/ifs/'
# input the year
# input the simulation tile
#tile="h10v04"

# function to read Lc and Lswi max
def VPMprmt (directory,tile):
    tempfile=None
    for path, subdirs, files in os.walk(directory):
        for name in files:
            if (tile in name) and (".tif" == name[-4:]) : tempfile=os.path.join(path,name)
    print tempfile
    inprmtdata=gdal.Open(tempfile)
    prmtdata=inprmtdata.ReadAsArray()
    return prmtdata


#Function to build vrt
def buildVrtFile (root,year,tile,product):
    fileList=[]
    for path, subdirs, files in os.walk(root):
        for name in files:
            if ".tif" == name[-4:]: fileList.append([os.path.join(path,name)])
    fileList.sort()
    print len(fileList),'files were built into a vrt file'
    filename=[]
    filename.append(os.path.join('/data/ifs/users/yzhang/TEMP/VRT',year+tile+product+'_list.txt'))
    outFilelist=open(filename[0],'w')
    for file in fileList:
        outFilelist.write(file[0]+'\r\n')
        filename.append(os.path.basename(file[0])[13:16])
    outFilelist.close()
    return filename


# Function to wirte array to tiff file
def write_file(output_name,output_array,GeoT,xsize,ysize,proJ,driverName='GTiff'):
    print "creating", output_name
    dr=gdal.GetDriverByName(driverName)
    dr.Register()
    do=dr.Create(output_name,xsize,ysize,1,gdal.GDT_UInt16,options = [ 'COMPRESS=LZW' ])

    do.SetGeoTransform(GeoT)
    do.SetProjection(proJ)
    do.GetRasterBand(1).WriteArray(output_array)
    do.GetRasterBand(1).SetNoDataValue(65535)
    do=None


# Function to wirte multi-dimesion array to tiff file
def export_array(Rasters,directory,type,tile,year,index):
    fileNum=Rasters.shape[0]
    for i in range(fileNum):
        fileName=os.path.join(directory,type+'.'+year+str(1000+index[i])[1:]+'.'+tile+'.tif')
        write_file(fileName,Rasters[i,:,:],geoTran,cols,rows,geoProj,"GTiff" )


def GPPprocess(tile):
    #input dirs 
    #VI
    #for VI using the SG_filled dataset
    dirEVI=root+'/VPM/driving_data/MOD09A1_006_BISE_SG/'+year+'/'+tile
    #dirEVI=root+'/modis/products/mod09a1/geotiff/evi/'+year+'/'+tile
	#########this is for C6 version but using C5 LSWI for 2000.
    #########use lswi_2000_2001 for missing tiles
    dirLSWI=root+'/modis/products_006/mod09a1/geotiff/lswi/'+year+'/'+tile
    if int(year) < 2002:
        dirLSWI=root+'/VPM/driving_data/lswi_2000_2001/'+year+'/'+tile
    #climate data R and T
    fileRAD = root+'/VPM/driving_data/NCEP2/DSSW/'+year+'/rad.'+tile+'.'+year+'.nc'
    fileTEMP = root+'VPM/driving_data/NCEP2/TaDay/'+year+'/temp.'+tile+'.'+year+'.nc'

    # The directory for LSWImax and landcover
    dirmaxLSWI=root+'/VPM/driving_data/LSWImax_MA_06/'+year
    dirLC = root+'/VPM/driving_data/MCD12Q1/'+year
    
    if int(year) < 2003:
        dirmaxLSWI=root+'/VPM/driving_data/LSWImax_MA_06/2003'
    if int(year) > 2014:
        dirmaxLSWI=root+'/VPM/driving_data/LSWImax_MA_06/2014'
    
    if int(year) < 2001:
        dirLC = root+'/VPM/driving_data/MCD12Q1/2001'
    if int(year) > 2013:
        dirLC = root+'/VPM/driving_data/MCD12Q1/2013'
    
    # read Lc and lswimax
    LC=VPMprmt(dirLC,tile)
    LSWImax=VPMprmt(dirmaxLSWI,tile)
    
    # read the C3/C4 percentage;;;; added by Yao on 04/28/2017
    C4grassdir = root+'/VPM/driving_data/C3C4/nature/'
    C4cropdir = root+'/VPM/driving_data/C3C4/crop/'
    C4grass = VPMprmt(C4grassdir,tile)
    C4crop = VPMprmt(C4cropdir,tile)
    C4grass[C4grass<0] = 0
    C4crop[C4crop<0] = 0

    if (LC!=None and LSWImax!=None):
        #output dir for GPP apar and epsilon g
        dirGPPTILE = root+'/VPM/product/v20/gpp/'+year+'/'+tile
        dirGPPpotTILE = root+'/VPM/product/v20/gpppot/'+year+'/'+tile
        dirstressTILE = root+'/VPM/product/v20/stress/'+year+'/'+tile
        
        if not os.path.exists(dirGPPTILE):
            os.makedirs(dirGPPTILE)
        if not os.path.exists(dirGPPpotTILE):
            os.makedirs(dirGPPpotTILE)
        if not os.path.exists(dirstressTILE):
            os.makedirs(dirstressTILE)
    
        #using LUT to get key parameters
        #replace all the fill value(255) with 0 (water)
        LC[LC==255]=0    
    
        LUTminT=numpy.array([0,-1,2,-1,-1,-1,-1,1,-1,1,0,0,-1,-1,0,0,0])
        LUToptT=numpy.array([30,20,28,20,20,19,25,31,24,30,27,20,30,27,27,20,30])
        LUTmaxT=numpy.array([48,40,48,40,40,48,48,48,48,48,48,48,48,48,48,40,48])
        LUTepsilon=numpy.array([float('nan'),0.078,0.078,0.078,0.078,0.078,0.078,0.078,0.078,0.078,0.078,0.078,0.078,0.078,0.078,0.078,0.078])
        LUTepsiloncrop=numpy.array([float('nan'),0,0,0,0,0,0,0,0,0,0,0,0.078,0,0.078,0,0])
        LUTepsilongrass=numpy.array([float('nan'),0,0,0,0,0,0,0,0.078,0.078,0.078,0.078,0,0,0,0,0])
        
        Tmin=LUTminT[LC]*100
        Tmax=LUTmaxT[LC]*100
        Topt=LUToptT[LC]*100
        LUE=LUTepsilon[LC]+LUTepsiloncrop[LC]*C4crop*0.5+LUTepsilongrass[LC]*C4grass*0.5
        
        # build EVI vrt file and read as an array
        EVIfile=buildVrtFile(dirEVI,year,tile,'evi')
        vrtEVI=os.path.join(os.path.dirname(EVIfile[0]),year+tile+'EVI_vrt.vrt')
        print "Building the vrt file: ", vrtEVI
        os.system('gdalbuildvrt -separate -input_file_list '+EVIfile[0]+' '+vrtEVI)
    
        # build LSWI vrt file and read as an array
        LSWIfile=buildVrtFile(dirLSWI,year,tile,'lswi')
        vrtLSWI=os.path.join(os.path.dirname(LSWIfile[0]),year+tile+'LSWI_vrt.vrt')
        print "Building the vrt file: ", vrtLSWI
        os.system('gdalbuildvrt -separate -input_file_list '+LSWIfile[0]+' '+vrtLSWI)
    
        #get the DOYs from the EVI file    
        EVIdoy=map(int,EVIfile[1:])
        LSWIdoy=map(int,LSWIfile[1:])
        
        # Define the OUTPUT Georeference information
    
        global rows, cols, geoProj,geoTran
    
        inEVI=gdal.Open(vrtEVI)
        print "reading the multi-EVI..."
        EVI=inEVI.ReadAsArray()
        EVI=1.25*(EVI*0.0001-0.1)
        EVI[EVI<0]=0
        EVI[EVI>1]=1
    
        rows = 2400
        cols = 2400
        geoTran=inEVI.GetGeoTransform()
        geoProj=inEVI.GetProjection()
    
    
        inLSWI=gdal.Open(vrtLSWI)
        print "reading the multi-LSWI file..."
        LSWI=inLSWI.ReadAsArray()
        
        print "reading the multi-NARR temperature file..."
        netCDFTEMP= Dataset(fileTEMP, 'r')
        TEMP=netCDFTEMP.variables['Temp'][:]
        TEMP=TEMP-27315 # from Kevin  to Celsius
    
        print "reading the multi-NARR radiation file..."
        netCDFRAD= Dataset(fileRAD, 'r')
        RAD=netCDFRAD.variables['PAR'][:]
        RAD=RAD/100 #convert the Radiation to w/m-2
        
    
        allexist=list(set(LSWIdoy).intersection(set(EVIdoy)))
        allexist.sort()
        if not (len(LSWIdoy)==len(EVIdoy)==46):    
    
            exist=(numpy.asarray(allexist)-1)/8
            EVIindex=[]
            LSWIindex=[]
            for doy in allexist:
                EVIindex.append(EVIdoy.index(doy))
                LSWIindex.append(LSWIdoy.index(doy))
            EVI=EVI[numpy.asarray(EVIindex),:,:]
            LSWI=LSWI[numpy.asarray(LSWIindex),:,:]
            TEMP=TEMP[exist,:,:]
            RAD=RAD[exist,:,:]
        
        
        Tscalar=((TEMP-Tmax)*(TEMP-Tmin))*1.000/(((TEMP-Tmax)*(TEMP-Tmin))-(TEMP-Topt)*(TEMP-Topt))
        Tscalar=numpy.where(Tscalar>1,1,Tscalar)
        Tscalar=numpy.where(Tscalar<0,0,Tscalar)
        
        Wscalar =(10000+LSWI)*1.000/(10000+LSWImax)    # if lswi > 0, use the original function
        #WscalarCROP=numpy.where(LSWI<=0,LSWI+LSWImax,Wscalar) # if lswi<=0, adjust the wscalar using lswi+lswimax
        Wscalar=numpy.where(Wscalar>1,1,Wscalar)
        Wscalar=numpy.where(Wscalar<0,0,Wscalar)
        #WscalarCROP=numpy.where(WscalarCROP<0,0,WscalarCROP)
        #WscalarCROP=numpy.where(WscalarCROP>1,1,WscalarCROP)
        #WscalarCROP=numpy.where((WscalarCROP>1)|(WscalarCROP<0),1,WscalarCROP)
    
        #scalar for GPPpot is 0.001 g C m-2 day-1
        GPPpot=LUE*EVI*RAD*1000
        #scalar for stress is 0.0001 
        stress=Tscalar*Wscalar*10000
        #scalar for GPP is 0.001 g C m-2 day-1
        GPP=GPPpot*stress/10000
        export_array (Rasters=numpy.uint(GPPpot),directory=dirGPPpotTILE,type='GPPpot',year=year,tile=tile,index=allexist)
        GPPpot=None
        export_array (Rasters=numpy.uint(stress),directory=dirstressTILE,type='stress',year=year,tile=tile,index=allexist)
        stress=None
        export_array (Rasters=numpy.uint(GPP),directory=dirGPPTILE,type='GPP',year=year,tile=tile,index=allexist)
        GPP=None
    

def process_list(tile= None, mp = True, count = 1):    
    if mp:
        #count = multiprocessing.cpu_count()-save_cpus
        pool = multiprocessing.Pool(processes=count)
        pool.map(GPPprocess, tile)


#
tile = ['h16v00','h17v00','h18v00','h19v00','h12v01','h13v01','h14v01','h15v01','h16v01','h17v01',\
    'h18v01','h19v01','h20v01','h21v01','h22v01','h23v01','h09v02','h10v02','h11v02','h12v02','h13v02',\
    'h14v02','h15v02','h16v02','h17v02',\
    'h18v02','h19v02','h20v02','h21v02','h22v02','h23v02','h24v02','h25v02','h26v02','h06v03','h07v03',\
    'h08v03','h09v03','h10v03','h11v03','h12v03','h13v03','h14v03','h15v03','h17v03','h18v03','h19v03',\
    'h20v03','h21v03','h22v03','h23v03','h24v03','h25v03','h26v03','h27v03','h28v03','h29v03','h08v04',\
    'h09v04','h10v04','h11v04','h12v04','h13v04','h14v04','h17v04','h18v04','h19v04','h20v04','h21v04',\
    'h22v04','h23v04','h24v04','h25v04','h26v04','h27v04','h28v04','h07v05','h08v05','h09v05','h10v05',\
    'h11v05','h12v05','h15v05','h16v05','h17v05','h18v05','h19v05','h20v05','h21v05','h22v05','h23v05',\
    'h24v05','h25v05','h26v05','h27v05','h28v05','h29v05','h30v05','h02v06','h03v06','h07v06','h08v06',\
    'h09v06','h10v06','h11v06','h16v06','h17v06','h18v06','h19v06','h20v06','h21v06','h22v06','h23v06',\
    'h24v06','h25v06','h26v06','h27v06','h28v06','h29v06','h30v06','h31v06','h01v07','h03v07','h07v07',\
    'h08v07','h09v07','h10v07','h11v07','h12v07','h15v07','h16v07','h17v07','h18v07','h19v07','h20v07',\
    'h21v07','h22v07','h23v07','h24v07','h25v07','h26v07','h27v07','h28v07','h29v07','h30v07','h31v07',\
    'h32v07','h33v07','h34v07','h00v08','h01v08','h02v08','h08v08','h09v08','h10v08','h11v08','h12v08',\
    'h13v08','h16v08','h17v08','h18v08','h19v08','h20v08','h21v08','h22v08','h23v08','h25v08','h26v08',\
    'h27v08','h28v08','h29v08','h30v08','h31v08','h32v08','h33v08','h34v08','h35v08','h00v09','h01v09',\
    'h02v09','h03v09','h04v09','h08v09','h09v09','h10v09','h11v09','h12v09','h13v09','h14v09','h16v09',\
    'h18v09','h19v09','h20v09','h21v09','h22v09','h23v09','h25v09','h27v09','h28v09','h29v09','h30v09',\
    'h31v09','h32v09','h33v09','h34v09','h35v09',\
    #southhemisphere
    'h00v10','h01v10','h02v10','h03v10','h04v10','h05v10','h10v10','h11v10','h12v10','h13v10','h14v10',\
    'h17v10','h19v10','h20v10','h21v10','h22v10','h23v10','h27v10','h28v10','h29v10','h30v10','h31v10',\
    'h32v10','h33v10','h34v10','h35v10','h01v11','h02v11','h03v11','h04v11','h05v11','h06v11','h08v11',\
    'h10v11','h11v11','h12v11','h13v11','h14v11','h15v11','h19v11','h20v11','h21v11','h22v11','h23v11',\
    'h27v11','h28v11','h29v11','h30v11','h31v11','h32v11','h33v11','h11v12','h12v12','h13v12','h16v12',\
    'h17v12','h19v12','h20v12','h24v12','h27v12','h28v12','h29v12','h30v12','h31v12','h32v12','h05v13',\
    'h12v13','h13v13','h17v13','h20v13','h21v13','h22v13','h28v13','h29v13','h30v13','h31v13','h13v14',\
    'h14v14','h15v14','h16v14','h18v14','h22v14','h27v14','h28v14']
#for year in ["2001","2002"]:
'''missing tiles for LC
'h11v01','h14v00','h15v00','h20v00','h21v00','h24v01',
'''

for year in ["2012","2013","2014"]:
	process_list(tile=tile,mp = True, count = 5) 