import os, sys

sys.path.append(r'D:\GIS\TOOLBOXES\pysda')

import sdapoly, sdaprop

import geopandas as gpd

# this aoi is too big, more than 10,000 records/ character limit for vertices
# myaoi = r"D:\geodata\project_data\gsp-gsocseq\AOIs\AOI_CONUS_bufbox.shp"

# os.chdir(r'D:\geodata\project_data\gsp-gsocseq\AOIs')

# local copy of eastern pacific islands polygons (with mukey!)
df = gpd.read_file(r'F:\ZBOOK\GIS\NRCS\R\EP_polys.shp')

# getthe column name for clay
sdaprop.numprops()

# send the request
c1 = sdaprop.getprop(df=df, column='MUKEY', prop = 'claytotal_r', method = 'wtd_avg', top = 0, bottom=30)

# should have reprojected to 4326 first
ndf = gpd.read_file(r'F:\ZBOOK\GIS\NRCS\R\EP_polys_4326.shp')

# make sure the cases ,atch
c1['MUKEY'] = c1['mukey']

# join the files
jn = ndf.merge(c1, on='MUKEY')

# write to new file for rasterizing
jn.to_file(r'C:\temp\gsoc_clay.shp', driver = "ESRI Shapefile")
