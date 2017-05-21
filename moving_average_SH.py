
#-------------------------------------------------------------------------------
# Name:        calculate the moving average of LSWI
# Inputs: LSWImax for each tile for each year (2001~2014)
#
# Author:      Yao Zhang
#
# Created:     12/01/2015
# Modified:    05/04/2016
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



# input the parent directory
root='/data/ifs/VPM/'
# input the yearwww.iplaysoft.com

#Function to build vrt
def buildVrtFile (year,tile):
    fileList=[]
    for yr in range(int(year)-2,int(year)+3):
        fileList.append([os.path.join('/data/ifs/users/xcwu/VPM_GPP/LSWImax/LSWI/',str(yr),tile+'.'+str(yr)+'.maxLSWI_5d_10d.tif')])
    #print len(fileList),'files were built into a vrt file'
    filename=os.path.join(dirLSWIMA,year+tile+'_list.txt')
    outFilelist=open(filename,'w')
    for file in fileList:
        outFilelist.write(file[0]+'\r\n')
    outFilelist.close()
    return filename

# Function to write array to tiff file
def write_file(output_name,output_array,GeoT,xsize,ysize,proJ,driverName='GTiff'):
    print "creating", output_name
    dr=gdal.GetDriverByName(driverName)
    dr.Register()
    do=dr.Create(output_name,xsize,ysize,1,gdal.GDT_Float32)
    do.SetGeoTransform(GeoT)
    do.SetProjection(proJ)
    do.GetRasterBand(1).WriteArray(output_array)
    do.GetRasterBand(1).SetNoDataValue(32767)
    do=None

# Function to calculate the moving average
def movingaverage(tile):
    # Output directories for moving average LSWImax
 
    # if the output directories don't exist, create the new directories
    if not os.path.exists(dirLSWIMA):
        os.makedirs(dirLSWIMA)
    # build LSWImax vrt file and read as an array
    file=buildVrtFile(year,tile)
    vrtLSWImax=os.path.join(os.path.dirname(file),year+tile+'LSWImax_vrt.vrt')
    #print "Building the vrt file: ", year+tile+vrtLSWImax
    os.system('gdalbuildvrt -separate -input_file_list '+file+' '+vrtLSWImax)
    
    global rows, cols, geoProj,geoTran
    inLSWImax=gdal.Open(vrtLSWImax)
    #print "reading the multi-LSWI..."
    LSWImax=inLSWImax.ReadAsArray()

    rows = 2400
    cols = 2400
    geoTran=inLSWImax.GetGeoTransform()
    geoProj=inLSWImax.GetProjection()
    
    #find the second largest LSWImax
    secLSWImax = numpy.sort(LSWImax, axis=0, kind='quicksort', order=None)[3,:,:]
    
    write_file(dirLSWIMA+'/'+tile+'.'+year+'.maxLSWI_MA.tif',secLSWImax,geoTran,rows,cols,geoProj,driverName='GTiff')
    
    secLSWImax=None
    LSWImax=None
 
def process_list(tiles = None, mp = True, save_cpus = 0):    
    if mp:
        count = multiprocessing.cpu_count()-save_cpus
        #manager = multiprocessing.Manager()
        #lock = manager.Lock()
        #map(lambda f: f.append(lock), tiles)
        pool = multiprocessing.Pool(processes=count)
        pool.map(movingaverage, tiles)

tile=['h00v09','h01v09','h02v09','h03v09','h04v09','h08v09','h09v09','h10v09','h11v09','h12v09','h13v09',\
    'h14v09','h16v09','h18v09','h19v09','h20v09','h21v09','h22v09','h23v09','h25v09','h27v09','h28v09',\
    'h29v09','h30v09','h31v09','h32v09','h33v09','h34v09','h35v09','h00v10','h01v10','h02v10','h03v10',\
    'h04v10','h05v10','h10v10','h11v10','h12v10','h13v10','h14v10','h17v10','h19v10','h20v10','h21v10',\
    'h22v10','h23v10','h27v10','h28v10','h29v10','h30v10','h31v10','h32v10','h33v10','h34v10','h35v10',\
    'h01v11','h02v11','h03v11','h04v11','h05v11','h06v11','h08v11','h10v11','h11v11','h12v11','h13v11',\
    'h14v11','h15v11','h19v11','h20v11','h21v11','h22v11','h23v11','h27v11','h28v11','h29v11','h30v11',\
    'h31v11','h32v11','h33v11','h11v12','h12v12','h13v12','h16v12','h17v12','h19v12','h20v12','h24v12',\
    'h27v12','h28v12','h29v12','h30v12','h31v12','h32v12','h05v13','h12v13','h13v13','h17v13','h20v13',\
    'h21v13','h22v13','h28v13','h29v13','h30v13','h31v13','h13v14','h14v14','h15v14','h16v14','h18v14',\
    'h22v14','h27v14','h28v14']
#for year in ['2003','2004','2005','2006','2007','2008','2009','2010','2011','2012']:
for year in ['2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']:
	dirLSWIMA=root+'/driving_data/LSWImax_MA_06/'+year
	process_list(tiles=tile, mp = True, save_cpus = 0)
