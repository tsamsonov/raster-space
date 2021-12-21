import numpy as np
import rspace
from osgeo import gdal, ogr, osr
from WBT.whitebox_tools import WhiteboxTools

PINF = -3.402823466e+38

wbt = WhiteboxTools()
wbt.set_whitebox_dir('/Users/tsamsonov/GitHub/raster-space/raster_space/WBT/')
wbt.set_verbose_mode(False)

outRasterSRS = osr.SpatialReference()
# srs = buildings.crs()
# wkt = srs.toWkt()
# outRasterSRS.ImportFromWkt(wkt)

drv = gdal.GetDriverByName('GTiff')

bld_ds = gdal.Open("/Volumes/Work/__UCLIM/Kosheleva/mask_blds_5m.tif")
srcband_buildings = bld_ds.GetRasterBand(1)
outblds = '/Users/tsamsonov/GitHub/raster-space/dist_buildings.tif'
bld_dist_ds = drv.Create(outblds,
                        bld_ds.RasterXSize, bld_ds.RasterYSize, 1,
                        gdal.GetDataTypeByName('Float32'))

bld_dist_ds.SetGeoTransform(bld_ds.GetGeoTransform())
bld_dist_ds.SetProjection(bld_ds.GetProjectionRef())
dstband_blds = bld_dist_ds.GetRasterBand(1)
gdal.ComputeProximity(srcband_buildings, dstband_blds, ["DISTUNITS=GEO"])

npdist_buildings = np.array(dstband_blds.ReadAsArray())  # length testing

buildings_alloc = '/Users/tsamsonov/GitHub/raster-space/bld_alloc.tif'
wbt.euclidean_allocation('/Users/tsamsonov/GitHub/raster-space/dist_buildings.tif', buildings_alloc)

hgtfile_blds = gdal.Open(buildings_alloc)
hgtband_blds = hgtfile_blds.GetRasterBand(1)
npheight_blds = np.array(hgtband_blds.ReadAsArray())

StepX = bld_ds.GetGeoTransform()[1]
nodata = -1.0
ndirs = 4
dist = 10000
threads = 4

npres = rspace.estimate_space(npdist_buildings, npheight_blds, npdist_buildings, npheight_blds,
                              StepX, nodata, ndirs, dist, threads, print, print)