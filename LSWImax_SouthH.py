
import os
import string
from os import listdir
from os.path import isfile, join
from osgeo import gdal
from osgeo.gdalconst import *
import numpy
import numpy.ma as ma


# input the parent directory
root='/data/ifs/users/xcwu'



#Function to build vrt
def buildVrtFile (root1,root2,product):
    fileList1=[]
    for path1, subdirs1, files1 in os.walk(root1):
        for name1 in files1:
            if ".tif" == name1[-4:]: fileList1.append([os.path.join(path1,name1)])
    fileList1.sort()
    print len(fileList1),'files were built into a vrt file'
    filename=os.path.join('/data/ifs/users/xcwu/VPM_GPP/LSWImax/temp',year+tile+product+'_list.txt')
    outFilelist=open(filename,'w')
    for file1 in fileList1[23:46]:
        outFilelist.write(file1[0]+'\r\n')
    fileList2=[]
    for path2, subdirs2, files2 in os.walk(root2):
        for name2 in files2:
            if ".tif" == name2[-4:]: fileList2.append([os.path.join(path2,name2)])
    fileList2.sort()
    print len(fileList2),'files were built into a vrt file'
    for file2 in fileList2[0:23]:
        outFilelist.write(file2[0]+'\r\n')	
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
''' 
    'h12v01','h13v01','h14v01','h15v01','h16v01','h17v01','h18v01','h19v01','h20v01','h21v01','h22v01',\
    'h23v01','h09v02','h10v02','h11v02','h12v02','h13v02','h14v02','h15v02','h16v02','h17v02',\
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
    'h31v09','h32v09','h33v09','h34v09','h35v09'
'''

#for year in ['2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014']:
for year in ['2003','2004','2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015']:
    for tile in ['h00v10','h01v10','h02v10','h03v10','h04v10','h05v10','h10v10','h11v10',
	'h12v10','h13v10','h14v10','h17v10','h19v10','h20v10','h21v10','h22v10','h23v10','h27v10',
	'h28v10','h29v10','h30v10','h31v10','h32v10','h33v10','h34v10','h35v10','h01v11','h02v11',
	'h03v11','h04v11','h05v11','h06v11','h08v11','h10v11','h11v11','h12v11','h13v11','h14v11',
	'h15v11','h19v11','h20v11','h21v11','h22v11','h23v11','h27v11','h28v11','h29v11','h30v11',
	'h31v11','h32v11','h33v11','h11v12','h12v12','h13v12','h16v12','h17v12','h19v12','h20v12',
	'h24v12','h27v12','h28v12','h29v12','h30v12','h31v12','h32v12','h05v13','h12v13','h13v13',
	'h17v13','h20v13','h21v13','h22v13','h28v13','h29v13','h30v13','h31v13','h13v14','h14v14',
	'h15v14','h16v14','h18v14','h22v14','h27v14','h28v14']:
        # Output directories for LSWImax
        dirLSWImax=root+'/VPM_GPP/LSWImax/LSWI/'+year
        # Output directories for SOS and EOS
        dirDOYSOSEOS=root+'/VPM_GPP/LSWImax/SOSEOS/'+year+'/DOYSOSEOS'
        # if the output directories don't exist, create the new directories
        if not os.path.exists(dirLSWImax):
            os.makedirs(dirLSWImax)
        if not os.path.exists(dirDOYSOSEOS):
            os.makedirs(dirDOYSOSEOS)
        # The directory for the multi-nighttime LST, LSWI, Cloud
        nextyear=string.atoi(year)+1		
        nextyear=str(nextyear)				
        dirLST1='/data/ifs/modis/products_006/myd11a2/'+year+'/'+tile
        dirLST2='/data/ifs/modis/products_006/myd11a2/'+nextyear+'/'+tile
        dirLSWI1='/data/ifs/modis/products_006/mod09a1/geotiff/lswi/'+year+'/'+tile
        dirLSWI2='/data/ifs/modis/products_006/mod09a1/geotiff/lswi/'+nextyear+'/'+tile
        dirNDVI1='/data/ifs/modis/products_006/mod09a1/geotiff/ndvi/'+year+'/'+tile
        dirNDVI2='/data/ifs/modis/products_006/mod09a1/geotiff/ndvi/'+nextyear+'/'+tile
        dirCloud1='/data/ifs/modis/products_006/mod09a1/geotiff/cloudmask/'+year+'/'+tile
        dirCloud2='/data/ifs/modis/products_006/mod09a1/geotiff/cloudmask/'+nextyear+'/'+tile
        # build LSWI vrt file and read as an array
        file=buildVrtFile(dirLSWI1,dirLSWI2,'lswi')
        vrtLSWI=os.path.join(os.path.dirname(file),year+tile+'LSWI_vrt.vrt')
        print "Building the vrt file: ", vrtLSWI
        os.system('gdalbuildvrt -separate -input_file_list '+file+' '+vrtLSWI)

        global rows, cols, geoProj,geoTran
        inLSWI=gdal.Open(vrtLSWI)
        print "reading the multi-LSWI..."
        LSWI=inLSWI.ReadAsArray()
        LSWIorg=LSWI
        rows = 2400
        cols = 2400
        geoTran=inLSWI.GetGeoTransform()
        geoProj=inLSWI.GetProjection()

        # build NDVI vrt file and read as an array
        file=buildVrtFile(dirNDVI1,dirNDVI2,'ndvi')
        vrtNDVI=os.path.join(os.path.dirname(file),year+tile+'NDVI_vrt.vrt')
        print "Building the vrt file: ", vrtNDVI
        os.system('gdalbuildvrt -separate -input_file_list '+file+' '+vrtNDVI)
        inNDVI=gdal.Open(vrtNDVI)
        print "reading the multi-NDVI..."
        NDVI=inNDVI.ReadAsArray()


        # build cloud vrt file and read as an array
        file=buildVrtFile(dirCloud1,dirCloud2,'cloud')
        vrtCloud=os.path.join(os.path.dirname(file),year+tile+'cloud_vrt.vrt')
        print "Building the vrt file: ", vrtCloud
        os.system('gdalbuildvrt -separate -input_file_list '+file+' '+vrtCloud)
        inCloud=gdal.Open(vrtCloud)
        print "reading the multi-Cloud..."
        Cloud=inCloud.ReadAsArray()

        # build nighttime LST vrt file and read as an array
        file=buildVrtFile(dirLST1,dirLST2,'ngtLST')
        vrtLST=os.path.join(os.path.dirname(file),year+tile+'ngtLST_vrt.vrt')
        print "Building the vrt file: ", vrtLST
        os.system('gdalbuildvrt -separate -input_file_list '+file+' '+vrtLST)
        inLST=gdal.Open(vrtLST)
        print "reading the multi-night LST..."
        LST=inLST.ReadAsArray()

        # Convert the MODIS LST using scalar and offset 0.02, 273
        LST=LST*0.02-273


        # calculate the end of growing season
        # The first day when nightLST <10 DEG during the 2nd half year
        LSTEOS=numpy.where((LST==-273)|(LST>=10),0, 1)


        # calculate the start of growing season
        # The first day when there are three points with nightLST >5 DEG
        LST=numpy.where(LST<0,0, LST)
        LST=LST/5
        LST=LST.astype(int)
        LST=numpy.where(LST>0,1, LST)
        iniLST=LST[0,:,:]*LST[1,:,:]*LST[2,:,:]
        for i in range(43):
            LST=numpy.roll(LST,-1,axis=0)
            tempLST=LST[0,:,:]*LST[1,:,:]*LST[2,:,:]
            iniLST=numpy.append(iniLST,tempLST,axis=0)
        iniLST=iniLST.reshape(44,1200,1200)
        # to calculate if all year nightLST < 5 DEG. Yes->-10 (no growing season/LSWI=0.25)
        maskSOS=numpy.sum(iniLST,axis=0)
        # calculate the SOS

        inSOSLST=numpy.argmax(iniLST,axis=0)
        inSOSLST=numpy.where(maskSOS==0,-10,inSOSLST)

        iniLST=None
        maskSOS=None

        # resample 1-km LST to 500 m
        temp=numpy.zeros((2400,2400))
        SOSLST = numpy.array([inSOSLST[x/2,y/2] for x, y in numpy.ndindex(temp.shape)])
        SOSLST=SOSLST.reshape(2400,2400)
        print "saving the SOS"
        write_file(dirDOYSOSEOS+'/'+tile+'.'+year+'.SOS_ngtLST_5d.tif',SOSLST,geoTran,rows,cols,geoProj,driverName='GTiff')
        temp=None

        print "start to calculate the EOS falling below 10 degree"

        iniLSTEOS=LSTEOS[23:,:,:]
        EOSLST = numpy.argmax(iniLSTEOS,axis=0)
        EOSLST=EOSLST+23

        # if iniLSTEOS==0 means that there is no date below <10 DEG, Use the last DOY as EOS
        midEOSLST=numpy.sum(iniLSTEOS,axis=0)
        EOSLST=numpy.where(midEOSLST==0,45,EOSLST )
        # here is no growing season
        EOSLST=numpy.where(inSOSLST==-10,-10,EOSLST )
        midEOSLST=None



        # resample 1-km LST to 500 m
        temp=numpy.zeros((2400,2400))
        EOSLSTOUT = numpy.array([EOSLST[x/2,y/2] for x, y in numpy.ndindex(temp.shape)])
        EOSLSTOUT=EOSLSTOUT.reshape(2400,2400)
        temp=None
        EOSLST=None
        print "saving the EOS below 10 degree"
        write_file(dirDOYSOSEOS+'/'+tile+'.'+year+'.EOS_ngtLST_10d.tif',EOSLSTOUT,geoTran,rows,cols,geoProj,driverName='GTiff')

        # for the v01 tiles, the scene numbers are not 46
        if tile in ['h16v01','h15v01','h14v01','h13v01','h12v01','h17v01','h18v01','h19v01','h20v01','h21v01','h22v01','h23v01']:
            before=numpy.zeros((2400,2400))
            LSWI=numpy.insert(LSWI,0,before,axis=0)
            Cloud=numpy.insert(Cloud,0,before,axis=0)
            LSWIorg=numpy.insert(LSWIorg,0,before,axis=0)
            NDVI=numpy.insert(NDVI,0,before,axis=0)
            before=None
##            after=numpy.arange(23040000)*0
##            after=after.reshape((4,2400,2400))
            after=numpy.arange(23040000)*0
            after=after.reshape((4,2400,2400))
            Cloud=numpy.append(Cloud,after,axis=0)
            LSWI=numpy.append(LSWI,after,axis=0)
            LSWIorg=numpy.append(LSWIorg,after,axis=0)
            NDVI=numpy.append(NDVI,after,axis=0)
            print LSWI.shape
            print Cloud.shape
            print NDVI.shape
            after=None
            
        #if LSWI.shape[1] == 45:
        #    middle=-1*numpy.ones((2400,2400))
        #    LSWI=numpy.insert(LSWI,22,middle,axis=0)
        #    Cloud=numpy.insert(Cloud,22,middle,axis=0)
        #    LSWIorg=numpy.insert(LSWIorg,22,middle,axis=0)
        # exclude the LSWI affected by cloud
        
        LSWI = ma.masked_where(Cloud > 1,LSWI)
        print "calculating the maximum LSWI"

        #Build (46,2400,2400) temporal indices
        temp=numpy.arange(46)
        temp=numpy.repeat(temp,2400*2400)
        #maskout the LSWI beyond the range [SOSLSTO
        temp=temp.reshape(46,2400,2400)

        # calculate the maximum NDVI during June to August
        NDVIgrowingseason=ma.masked_where((temp < 23)|(temp > 30),NDVI)
        NDVImaxday=numpy.argmax(NDVIgrowingseason,axis=0)
        # if all year nighttime LST<5 DEG, get LSWImax on the day of NDVImax.
        SOSLSTOUT=numpy.where(SOSLST==-10,NDVImaxday-1,SOSLST)
        EOSLSTOUT1=EOSLSTOUT
        EOSLSTOUT=numpy.where(EOSLSTOUT==-10,NDVImaxday+1,EOSLSTOUT)
        # if SOS>EOS, SOS=NDVI max day -1, EOS=NDVI max day +1
        SOSLSTOUT=numpy.where(SOSLSTOUT>EOSLSTOUT,NDVImaxday-1,SOSLSTOUT)
        EOSLSTOUT=numpy.where(SOSLSTOUT>EOSLSTOUT,NDVImaxday+1,EOSLSTOUT)

        LSWI=ma.masked_where((temp < SOSLSTOUT)|(temp > EOSLSTOUT), LSWI)
        LSWIorg=ma.masked_where((temp < SOSLSTOUT)|(temp > EOSLSTOUT), LSWIorg)
        #return the LSWI maximum without considering NAA
        maxLSWI=numpy.nanmax(LSWI,axis=0)
        maxLSWIorg=numpy.nanmax(LSWIorg,axis=0)
        temp=None
        SOSLSTOUT=None
        # Set 0.25 to fill the no growing season, LSWI not [1,1]
        maxLSWI=numpy.where((maxLSWI>=10000)|(maxLSWI<=-10000),maxLSWIorg,maxLSWI)
        EOSLSTOUT=None
        EOSLSTOUT1=None
        write_file(dirLSWImax+'/'+tile+'.'+year+'.maxLSWI_5d_10d.tif',maxLSWI,geoTran,rows,cols,geoProj,driverName='GTiff')
        maxLSWI=None
        LSWI=None




