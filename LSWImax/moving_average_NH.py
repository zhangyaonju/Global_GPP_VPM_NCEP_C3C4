
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
# input the year

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

		
''''h09v02','h10v02','h11v02','h12v02','h13v02','h14v02','h15v02','h16v02','h17v02',\
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
    'h27v08','h28v08','h29v08','h30v08','h31v08','h32v08','h33v08','h34v08','h35v08'		
'''		


tile=['h11v01','h12v01','h13v01','h14v01','h15v01','h16v01','h17v01','h18v01','h19v01','h20v01','h21v01',\
    'h22v01','h23v01','h24v01','h14v00','h15v00','h16v00','h17v00','h18v00','h19v00','h20v00','h21v00',\
    'h09v02','h10v02','h11v02','h12v02','h13v02','h14v02','h15v02','h16v02','h17v02',\
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
    'h27v08','h28v08','h29v08','h30v08','h31v08','h32v08','h33v08','h34v08','h35v08']

for year in ['2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']:
	dirLSWIMA=root+'/driving_data/LSWImax_MA_06/'+year
	process_list(tiles=tile, mp = True, save_cpus = 0)
