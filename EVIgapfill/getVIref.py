#-------------------------------------------------------------------------------
# Name:        Preprocessing for the EVI reference 
# Inputs:      EVI for each 8-day from all tiles and quality layers
#
# Author:      Yao Zhang
#
# Created:     3/29/2017
# Modified:    
# Copyright:   (c) eomf 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import multiprocessing
import os
from os import listdir
from os.path import isfile, join
from osgeo import gdal
from osgeo.gdalconst import *
from scipy.signal import savgol_filter
from ctypes import *
import numpy as np
import numpy.ma as ma
#from netCDF4 import Dataset
import time
import pandas as pd

startTime = time.time()

root = '/data/ifs/modis/products_006/mod09a1/geotiff/'
dirref = '/data/ifs/VPM/driving_data/EVI_ref/'

'''
def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]
'''

def VIsmooth(ndvi):
    rdays = c_int(len(ndvi))
    fun = cdll.LoadLibrary(os.getcwd() + '/bise.so')
    outndvi = (c_double * len(ndvi))()
    slidingperiod = c_int(np.sum(ndvi == 0)/(len(ndvi)/20))
    #apply the bise algorithm
    fun.bise(byref(rdays), (c_double * len(ndvi))(*ndvi), byref(slidingperiod), outndvi)
    bisendvi = np.frombuffer(outndvi)
    #print bisendvi
    bisendvi[bisendvi == -1] = np.nan
    peaks = []
    threshold = 1.5
    check = np.argwhere(np.isnan(bisendvi))
    #print check
    if len(check) < 3:
        return ndvi
    else:
        for i in range(0, len(check)):
            if i == 0:
                if bisendvi[check[i]] > (threshold * np.mean(bisendvi[np.array([(check[len(check)-1], check[i+1])])])):
                    if bisendvi[check[i]] > 3000: peaks.append(check[i])
            else:
                if i == (len(check)-1):
                    if bisendvi[check[i]] > (threshold * np.mean(bisendvi[np.array([(check[i-1], check[1])])])):
                        if bisendvi[check[i]] > 3000: peaks.append(check[i])
                else:
                    if bisendvi[check[i]] > (threshold * np.mean(bisendvi[check[np.array([i-1, i+1])]])):
                        if bisendvi[check[i]] > 3000: peaks.append(check[i])
        bisendvi[peaks] = np.nan
    return bisendvi

    #


def buildVrtFile(root, doy, tile, product):
    fileList = []
    for year in range(2000, 2017):
        tiledir = os.path.join(root, product, str(year), tile)
        for path, subdirs, files in os.walk(tiledir):
            for name in files:
                if (str(1000+doy)[1:] == name[13:16]) and (".tif" == name[-4:]): fileList.append([os.path.join(path, name)])
    fileList.sort()
    print len(fileList), 'files were built into a vrt file'
    if len(fileList) == 0: return 0
    filename = os.path.join('/data/ifs/users/yzhang/TEMP/VRT', str(1000+doy)[1:]+tile+product+'_list.txt')
    outFilelist = open(filename, 'w')
    for file in fileList:
        outFilelist.write(file[0]+'\r\n')
    outFilelist.close()
    return filename

def write_file(output_name, output_array, GeoT, xsize, ysize, proJ, driverName='GTiff'):
    print "creating", output_name
    dr = gdal.GetDriverByName(driverName)
    dr.Register()
    do = dr.Create(output_name, xsize, ysize, 1, gdal.GDT_UInt16, options=['COMPRESS=LZW'])

    do.SetGeoTransform(GeoT)
    do.SetProjection(proJ)
    do.GetRasterBand(1).WriteArray(output_array)
    do.GetRasterBand(1).SetNoDataValue(32767)
    do = None

def export_array (Rasters, directory, prod, tile, index):
    fileNum = Rasters.shape[0]
    for i in range(fileNum):
        fileName=os.path.join(directory, 'MOD09A1.'+str(1000+index[i])[1:]+'.'+tile+'.'+prod+'.tif')
        write_file(fileName, Rasters[i, :, :], geoTran, cols, rows, geoProj, "GTiff")

def parallelize_dataframe(df):
    df_split = np.array_split(df, 5, axis=1)
    pool = multiprocessing.Pool(5)
    df = np.concatenate(pool.map(dataframeapply, df_split), axis=1)
    pool.close()
    pool.join()
    return df

def dataframeapply(df):
    df = pd.DataFrame(np.concatenate([df[23:46, :], df, df[:23, :]]))
    df_smoothed = df.apply(VIsmooth)
    df_smoothed = df_smoothed.interpolate(axis=0)
    #make a SG filter
    df_select = df_smoothed.as_matrix()[23:69, :]
    df_select[np.isnan(df_select)] = 0
    bisendviSG = savgol_filter(df_select, window_length=5, polyorder=3)
    #bisendvi = None
    bisendviSG[bisendviSG < 0] = 0
    return bisendviSG

def import_all_year_data(tile):
    temp = np.zeros([46, 2400*2400], np.dtype(float))
    if int(tile[5:6])<2:
        temp[:]=np.nan
    for doy in range(1, 369, 8):
        evifile = buildVrtFile(root, doy, tile, 'evi')
        cloudfile = buildVrtFile(root, doy, tile, 'cloudmask')
        aerosolfile = buildVrtFile(root, doy, tile, 'aerosolmask')
        #if no file found for this DOY
        if evifile == 0: continue
        #doyList.append(doy)
        #build vrt for EVI
        vrtEVI = os.path.join(os.path.dirname(evifile), str(1000+doy)[1:]+tile+'EVI_vrt.vrt')
        print "Building the vrt file: ", evifile
        os.system('gdalbuildvrt -separate -input_file_list '+evifile+' '+vrtEVI)
        inEVI = gdal.Open(vrtEVI)
        EVI = inEVI.ReadAsArray()
        #build vrt for cloudmask
        vrtcloud = os.path.join(os.path.dirname(cloudfile), str(1000+doy)[1:]+tile+'cloud_vrt.vrt')
        print "Building the vrt file: ", cloudfile
        os.system('gdalbuildvrt -separate -input_file_list '+cloudfile+' '+vrtcloud)
        incloud = gdal.Open(vrtcloud)
        cloud = incloud.ReadAsArray()
        #build vrt for aerosol
        vrtaerosol = os.path.join(os.path.dirname(aerosolfile), str(1000+doy)[1:]+tile+'aerosol_vrt.vrt')
        print "Building the vrt file: ", aerosolfile
        os.system('gdalbuildvrt -separate -input_file_list '+aerosolfile+' '+vrtaerosol)
        inaerosol = gdal.Open(vrtaerosol)
        aerosol = inaerosol.ReadAsArray()
        global rows, cols, geoProj, geoTran
        rows = 2400
        cols = 2400
        geoTran = inEVI.GetGeoTransform()
        geoProj = inEVI.GetProjection()
		#mask for bad quality
        EVIgood = ma.masked_where((cloud != 1)|(aerosol == 0)|(EVI < 0)|(EVI > 10000), EVI)
        EVIgood = EVIgood.reshape(EVIgood.size/2400/2400, 2400*2400)
        medianEVI = np.nanmedian(EVIgood, axis=0)
        EVI = None
        aerosol = None
        cloud = None
        EVIgood = None
	    #assign to the 46 layer of matrix
        temp[(doy-1)/8, :] = medianEVI
        meanEVI = None
    return temp



def smooth(tile):
    #first use this function to get mean and save it in an array
    temp = import_all_year_data(tile)
    ####after get the mean value for all doy, I will run a bise gapfill first
    print temp.size
    ##when using the single processing
    #inputVI = pd.DataFrame(temp)
    #VIsmoothed = inputVI.apply(VIsmooth, axis=0)
    #VIsmoothed = VIsmoothed.as_matrix()
    #VIsmoothed = parallelize_dataframe(temp)
    ##when using the multiprocessing
    VIsmoothed = dataframeapply(temp)
    VIsmoothed = VIsmoothed.reshape(VIsmoothed.size/2400/2400, 2400, 2400)
    TILEdir = os.path.join(dirref, tile)
    if not os.path.exists(TILEdir):
        os.makedirs(TILEdir)
    export_array (Rasters=np.int16(VIsmoothed), directory=TILEdir, \
        prod='EVI.BISE.SG', tile=tile, index=range(1, 369, 8))
    temp = None
    inputVI = None
    VIsmoothed = None


def process_list(tile=None, mp=True, count=1):
    if mp:
        #count = multiprocessing.cpu_count()-save_cpus
        pool = multiprocessing.Pool(processes=count)
        pool.map(smooth, tile)

#

'''
tile = ['h17v00','h12v01','h13v01','h14v01','h15v01','h16v01','h17v01','h18v01','h19v01','h20v01','h21v01',\
    'h22v01','h23v01','h09v02','h10v02','h11v02','h12v02','h13v02','h14v02','h15v02','h16v02','h17v02',\
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
    '''
'''
tile = ["h17v00","h13v01","h10v02","h21v02","h22v02","h20v04","h21v04","h23v04",\
    "h24v04","h27v04","h08v05","h10v05","h11v05","h17v05","h19v05","h20v05","h21v05",\
    "h22v05","h23v05","h24v05","h25v05","h26v05","h27v05","h28v05","h29v05","h30v05",\
    "h02v06","h03v06","h07v06","h08v06","h09v06","h11v06","h16v06","h17v06","h18v06",\
    "h19v06","h20v06","h21v06","h22v06","h23v06","h24v06","h25v06","h26v06","h27v06",\
    "h28v06","h29v06","h30v06","h31v06","h01v07","h03v07","h08v07","h12v07","h24v07"]

tile = ["h00v10","h01v10","h02v06","h02v10","h03v10","h04v10","h07v05","h08v06","h09v05",\
    "h09v06","h10v05","h10v09","h10v10","h11v01","h11v05","h11v10","h12v09","h12v10",\
    "h13v10","h14v00","h14v04","h14v10","h15v00","h16v00","h17v04","h17v05","h18v00",\
    "h18v06","h19v00","h19v04","h19v05","h19v06","h20v00","h20v06","h21v00","h21v05",\
    "h21v06","h21v10","h22v04","h22v05","h22v06","h23v04","h23v05","h23v06","h23v09",\
    "h24v01","h24v05","h25v04","h25v05","h25v09","h26v04","h27v04","h27v05","h27v10",\
    "h28v04","h28v09","h28v10","h29v09","h29v10","h30v05","h30v09","h30v10","h31v10",\
    "h32v10","h35v09"]
'''
tile = ["h11v01","h12v01","h13v01","h14v00","h14v01","h15v00","h15v01","h16v00","h16v01",\
    "h17v00","h17v01","h18v00","h18v01","h19v00","h19v01","h20v00","h20v01","h21v00",\
    "h21v01","h22v01","h23v01","h24v01"]

#i = np.arange(0,5)
#segtiles = tile[0:60] #lotus
#segtiles = tile[60:120]   #for peony
#segtiles = tile[120:180] #for cattle
#segtiles = tile[180:240]  # crane
#segtiles = tile[240:287] #lily
#segtiles = tile[0:12] #for cattle
#segtiles = tile[12:24] #for lily
#segtiles = tile[24:36] #for crane
#segtiles = tile[36:48] #for lotus
#segtiles = tile[48:65] #for poeny
#smooth(segtiles)
#process_list(segtiles, mp=True, count=6)
#segtiles = tile[0:5] #for cattle
#segtiles = tile[5:10] #for lily
#segtiles = tile[10:15] #for crane
#segtiles = tile[15:20] #for lotus
segtiles = tile[20:22] #for poeny
process_list(segtiles, mp=True, count=5)
