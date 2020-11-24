from osgeo import gdal
import numpy as np

lulc= gdal.Open('LULC.tif')
field= lulc.ReadAsArray()
mapper= {
	11: 0.040,
	21: 0.0404,
	22: 0.0678,
	23: 0.0678,
	24: 0.0404,
	31: 0.0113,
	41: 0.36,
	42: 0.32,
	43: 0.40,
	52: 0.40,
	71: 0.368,
	81: 0.325,
	82: 0.035,
	90: 0.086,
	95: 0.1825
	}
	
manning= np.ones(field.shape)*-9999.
for key in mapper:
	rows,cols= np.where(field==key)
	manning[rows,cols]=mapper[key]

def arr2tif(dst, arr, sample):
	rows,cols= arr.shape
	
	driver = gdal.GetDriverByName("GTiff")
	outdata = driver.Create(dst, cols, rows, 1, gdal.GDT_Float32)
	outdata.SetGeoTransform(sample.GetGeoTransform())##sets same geotransform as input
	outdata.SetProjection(sample.GetProjection())##sets same projection as input
	outdata.GetRasterBand(1).WriteArray(arr)
	outdata.GetRasterBand(1).SetNoDataValue(-9999.)##if you want these values transparent
	outdata.FlushCache() ##saves to disk!!
	outdata = None

arr2tif('manningn.tif', manning, lulc)
