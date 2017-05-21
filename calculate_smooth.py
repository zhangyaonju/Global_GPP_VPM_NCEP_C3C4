#-------------------------------------------------------------------------------
# Name:        comparing the yearly and the reference to get the smoothed VI
# Inputs:      1. EVI for each 8-day from all tiles and quality layers
#              2. EVI reference calculated from getVIref.py
# Author:      Yao Zhang
#
# Created:     4/01/2017
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
import gc
#from getVIref import write_file

#startTime = time.time()

global root, dirref, smootheddir
root = '/data/ifs/modis/products_006/mod09a1/geotiff/'
dirref = '/data/ifs/VPM/driving_data/EVI_ref/'
smootheddir = '/data/ifs/VPM/driving_data/MOD09A1_006_BISE_SG'

def nan_helper(y):
    return np.isnan(y), lambda z: z.nonzero()[0]

def write_file(output_name, output_array, GeoT, xsize, ysize, proJ, driverName='GTiff'):
    print "creating", output_name
    dr = gdal.GetDriverByName(driverName)
    dr.Register()
    do = dr.Create(output_name, xsize, ysize, 1, gdal.GDT_UInt16, options=['COMPRESS=LZW'])

    do.SetGeoTransform(GeoT)
    do.SetProjection(proJ)
    do.GetRasterBand(1).WriteArray(output_array)
    do.GetRasterBand(1).SetNoDataValue(65535)
    do = None

def export_array(Rasters, directory, prod, year, tile, index):
    fileNum = Rasters.shape[0]
    for i in range(fileNum):
        fileName = os.path.join(directory, 'MOD09A1.A'+str(year)+str(1000+index[i])[1:]+'.'+\
            tile+'.'+prod+'.tif')
        write_file(fileName, Rasters[i, :, :], geoTran, cols, rows, geoProj, "GTiff")


def buildVrtFile(dir, year, tile, product):
    fileList = []
    if year == 2000:
        #pre year for 2000, we use 2000
        dirprod = os.path.join(dir, product, str(year), tile)
        for path, subdirs, files in os.walk(dirprod):
            for name in files:
                if (".tif" == name[-4:]) and (int(name[13:16]) > 177):
                    fileList.append([os.path.join(path, name)])
        fileList.sort()
        #add current year 2000
        newlist = []
        for name in files:
            if ".tif" == name[-4:]:
                newlist.append([os.path.join(path, name)])
        newlist.sort()
        fileList = fileList + newlist
        #add next year year+1
        newlist = []
        dirprod = os.path.join(dir, product, str(year+1), tile)
        for path, subdirs, files in os.walk(dirprod):
            for name in files:
                if (".tif" == name[-4:]) and (int(name[13:16]) < 185):
                    newlist.append([os.path.join(path, name)])
        newlist.sort()
        fileList = fileList + newlist
    if year == 2016:
        #pre year for 2016, we use year-1
        dirprod = os.path.join(dir, product, str(year-1), tile)
        for path, subdirs, files in os.walk(dirprod):
            for name in files:
                if (".tif" == name[-4:]) and (int(name[13:16]) > 177):
                    fileList.append([os.path.join(path, name)])
        #add current year 2016
        dirprod = os.path.join(dir, product, str(year), tile)
        for path, subdirs, files in os.walk(dirprod):
            for name in files:
                if ".tif" == name[-4:]:
                    fileList.append([os.path.join(path, name)])
        fileList.sort()
        #add next year, still use 2016
        newlist = []
        for name in files:
            if (".tif" == name[-4:]) and (int(name[13:16]) < 185):
                newlist.append([os.path.join(path, name)])
        newlist.sort()
        fileList = fileList + newlist
    if (year != 2000) and (year != 2016):
        dirprod = os.path.join(dir, product, str(year-1), tile)
        for path, subdirs, files in os.walk(dirprod):
            for name in files:
                if (".tif" == name[-4:]) and (int(name[13:16]) > 177):
                    fileList.append([os.path.join(path, name)])
        dirprod = os.path.join(dir, product, str(year), tile)
        for path, subdirs, files in os.walk(dirprod):
            for name in files:
                if ".tif" == name[-4:]:
                    fileList.append([os.path.join(path, name)])
        dirprod = os.path.join(dir, product, str(year+1), tile)
        for path, subdirs, files in os.walk(dirprod):
            for name in files:
                if (".tif" == name[-4:]) and (int(name[13:16]) < 185):
                    fileList.append([os.path.join(path, name)])
        fileList.sort()
    print len(fileList), 'files were built into a vrt file'
    filename = []
    filename.append(os.path.join('/data/ifs/users/yzhang/TEMP/VRT',\
        str(year)+tile+product+'_list.txt'))
    outFilelist = open(filename[0], 'w')
    for file in fileList:
        outFilelist.write(file[0]+'\r\n')
        filename.append(os.path.basename(file[0])[13:16])
    outFilelist.close()
    return filename



def buildrefVrtFile(tile):
    #get ref VI
    fileList = []
    for path, subdirs, files in os.walk(os.path.join(dirref, tile)):
        for name in files:
            if ".tif" == name[-4:]:
                fileList.append([os.path.join(path, name)])
    fileList.sort()
    filename = os.path.join('/data/ifs/users/yzhang/TEMP/VRT', 'ref.'+tile+'_list.txt')
    outFilelist = open(filename, 'w')
    for file in fileList:
        outFilelist.write(file[0]+'\r\n')
    outFilelist.close()
    return filename


def VIsmooth_ref(x):
    #the size of EVIgood is 92*5760000, the size of the reference data is 46*5760000
    x[x == -9999] = np.nan
    EVIgood = x[0:92]
    reference = np.concatenate([x[115:], x[92:], x[92:115]])
    if np.sum(np.isnan(EVIgood)) == 92:
        return np.concatenate([x[92:], x[23:69], x[92:]])
    ############################
    #here require complicated algorithm
    #first get the difference between these two
    diff = EVIgood - reference
    #fun = cdll.LoadLibrary(os.getcwd() + '/bise.so')
    #outdiff = (c_double * len(EVIgood))()
    #nans, y = nan_helper(diff)
    #diff[nans] = np.interp(y(nans), y(~nans), diff[~nans])
    diff[reference == 0] = 0
    diff = pd.Series(diff)
    reconstructVI = reference+diff.interpolate()
    SGVI = savgol_filter(np.array(reconstructVI[23:69]), window_length=5, polyorder=3)
    SGVI[SGVI < 0] = 0
    return np.concatenate([SGVI, x[23:69], x[92:]])

def reject_outliers(data, m=3):
    data[abs(data - np.nanmean(data, axis=0)) > m * np.nanstd(data, axis=0)] = np.nan
    return data

def sep_sg(df):
    #df = df.as_matrix()
    df[np.isnan(df)] = 0
    df = savgol_filter(df, axis=0, window_length=5, polyorder=3)
    return df
'''
def parallelize_dataframe(df, func):
    df_split = np.array_split(df, 5, axis=1)
    pool = multiprocessing.Pool(5)
    df = np.concatenate(pool.map(func, df_split), axis=1)
    pool.close()
    pool.join()
    return df
'''
def gapfill_VI(tile):
    #read in the reference data
    refVIfile = buildrefVrtFile(tile)
    vrtrefVI = os.path.join(os.path.dirname(refVIfile),\
        tile+'vrtrefVI_vrt.vrt')
    print "Building the vrt file: ", refVIfile
    os.system('gdalbuildvrt -separate -input_file_list '+refVIfile+' '+vrtrefVI)
    inrefVI = gdal.Open(vrtrefVI)
    refVI = inrefVI.ReadAsArray()
    refVI = refVI.reshape(46, 5760000)
    #read the VI for each year, do the gap filling, and then
    for year in range(2000, 2017):
        #first use this function to get mean and save it in an array
        temp = np.ones([92, 2400*2400], np.dtype(float))*(-9999)
        evifile = buildVrtFile(root, year, tile, 'evi')
        cloudfile = buildVrtFile(root, year, tile, 'cloudmask')
        aerosolfile = buildVrtFile(root, year, tile, 'aerosolmask')
        #if no file found for this year
        if evifile == 0: continue
        #doyList.append(doy)
        #build vrt for EVI
        vrtEVI = os.path.join(os.path.dirname(evifile[0]),\
            str(year)+tile+'EVI_vrt.vrt')
        print "Building the vrt file: ", vrtEVI
        os.system('gdalbuildvrt -separate -input_file_list '+evifile[0]+' '+vrtEVI)
        inEVI = gdal.Open(vrtEVI)
        EVI = inEVI.ReadAsArray()
        EVIdoy = np.array(map(int, evifile[1:]))
        #build vrt for cloudmask
        vrtcloud = os.path.join(os.path.dirname(cloudfile[0]),\
            str(year)+tile+'cloud_vrt.vrt')
        print "Building the vrt file: ", vrtcloud
        os.system('gdalbuildvrt -separate -input_file_list '+cloudfile[0]+' '+vrtcloud)
        incloud = gdal.Open(vrtcloud)
        cloud = incloud.ReadAsArray()
        #build vrt for aerosol
        vrtaerosol = os.path.join(os.path.dirname(aerosolfile[0]),\
            str(year)+tile+'aerosol_vrt.vrt')
        print "Building the vrt file: ", vrtaerosol
        os.system('gdalbuildvrt -separate -input_file_list '+aerosolfile[0]+' '+vrtaerosol)
        inaerosol = gdal.Open(vrtaerosol)
        aerosol = inaerosol.ReadAsArray()
        global rows, cols, geoProj, geoTran
        rows = 2400
        cols = 2400
        geoTran = inEVI.GetGeoTransform()
        geoProj = inEVI.GetProjection()
        #mask for bad quality
        EVIgood = ma.masked_where((cloud != 1)|(aerosol == 0)|(EVI < 0)|(EVI > 10000), EVI)
        EVIgood.set_fill_value(-9999)
        EVIgood = ma.filled(EVIgood.reshape(EVIgood.size/2400/2400, 2400*2400))
        EVI = None
        aerosol = None
        cloud = None
        #put the EVIgood into a 92 layer matrix
        #temp = np.empty([92, 2400*2400], np.dtype(float))
        EVIdoy = (EVIdoy-1)/8
        doy3year = np.empty(len(EVIdoy), np.dtype(int))
        k = 0
        doy3year[0] = EVIdoy[0]
        for i in range(1, (len(EVIdoy))):
            if EVIdoy[i] < EVIdoy[i-1]:
                k = k+46
            doy3year[i] = k + EVIdoy[i]
        #this is position of the layers that should be inserted into the target
        doy3year = doy3year - 23
        #this is the indicator of the layer of the EVIgood
        lid = np.arange(0, len(doy3year))
        temp[doy3year[(doy3year >= 0)&(doy3year < 92)], :] = \
            EVIgood[lid[(doy3year >= 0)&(doy3year < 92)], :]
        EVIgood = None
        #apply the new gapfill function
        ######################################
        ##this region is for new algorithm
        temp[temp == -9999] = np.nan
        temp[0:23, :] = temp[0:23, :] - refVI[23:46, :]
        temp[23:69, :] = temp[23:69, :] - refVI
        temp[69:, :] = temp[69:, :] - refVI[0:23, :]
        ## we need to gap fill the NA value first spikes!
        #reject outliers > than 3*sigma
        temp = reject_outliers(temp)
        diffVI = pd.DataFrame(temp)
        temp = None
        gc.collect()
        diffVI = diffVI.interpolate()
        diffVI[np.isnan(diffVI)] = 0
        reconstructVI = refVI + diffVI.as_matrix()[23:69, :]
        VIsmoothed = sep_sg(reconstructVI)
        #VIsmoothed = parallelize_dataframe(reconstructVI, sep_sg)
        '''
        pool2 = multiprocessing.Pool(processes=10)
        df_split = np.array_split(reconstructVI, 10, axis=1)
        #SG_VI = savgol_filter(reconstructVI[:,0:57600], axis=0, window_length=5, polyorder=3)
        df = pandas.concat(pool2.map(sep_sg, df_split))
        pool2.close()
        pool2.join()
        diffVI = diffVI.interpolate()
        reconstructVI = refVI + diffVI.as_matrix()[23:69,:]
        reconstructVIgood = ma.masked_where((refVI == 0)|(reconstructVI < 0), reconstructVI)
        reconstructVIgood.set_fill_value(0)
        reconstructVIgood = ma.filled(reconstructVIgood.reshape(46, 2400*2400))
        #####
        inputVI = pd.DataFrame(reconstructVIgood[:,1:57600])
        VIsmoothed = inputVI.apply(savgol_filter, axis=0, args=(5, 3))'''
        #VIsmoothed = VIsmoothed.as_matrix()
        VIsmoothed = VIsmoothed.reshape(VIsmoothed.size/5760000, 2400, 2400)
        TILEdir = os.path.join(smootheddir, str(year), tile)
        if not os.path.exists(TILEdir):
            os.makedirs(TILEdir)
        export_array(Rasters=np.int16(VIsmoothed), directory=TILEdir,\
            prod='EVI.BISE.SG', year=year, tile=tile, index=range(1, 369, 8))
        inputVI = None
        VIsmoothed = None
        temp = None


def process_list(tile=None, mp=True, count=1):
    if mp:
        pool = multiprocessing.Pool(processes=count)
        pool.map(gapfill_VI, tile)


tile = ['h14v00','h15v00','h16v00','h17v00','h18v00','h19v00','h20v00','h21v00','h11v01','h12v01','h13v01',\
    'h14v01','h15v01','h16v01','h17v01','h18v01','h19v01','h20v01','h21v01','h22v01','h23v01','h24v01',\
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
tile = ["h07v06","h09v06","h11v01","h11v06","h14v00","h15v00","h16v00","h17v04","h17v06",\
    "h18v00","h18v04","h18v06","h19v00","h19v06","h20v00","h20v04","h21v00","h21v04","h21v06",\
    "h23v04","h24v01","h24v04","h26v04","h27v04"]

tile = ["h11v01","h12v01","h13v01","h14v00","h14v01","h15v00","h15v01","h16v00","h16v01",\
    "h17v00","h17v01","h18v00","h18v01","h19v00","h19v01","h20v00","h20v01","h21v00",\
    "h21v01","h22v01","h23v01","h24v01"]
'''

#segtiles = tile[0:60] #lotus
#segtiles = tile[60:120]   #for peony
#segtiles = tile[120:180] #for cattle
#segtiles = tile[180:240]  # crane
segtiles = tile[240:287] #lily
#process_list(segtile, mp=False)
#segtiles = tile[0:5] #lotus
#segtiles = tile[5:10] #for lily
#segtiles = tile[10:15] #peony
#segtiles = tile[15:20] #cattle
#segtiles = tile[20:22] #crane
process_list(segtiles, mp=True, count=5)
#gapfill_VI("h11v08")