library(raster)
library(sf)

setwd("D:\\geodata")

epgsoc2 <- raster("D:\\geodata\\project_data2\\gsp-gsocseq\\CONUS_GSOCmap1.5.0.tif")
prvigsoc2 <- raster("D:\\geodata\\project_data2\\gsp-gsocseqp\\CONUS_GSOCmap1.5.0.tif")
wpgsoc2 <- raster("D:\\geodata\\project_data2\\gsp-gsocseqwp\\CONUS_GSOCmap1.5.0.tif")



gdalUtilities::gdalwarp(
  srcfile = "EP_clay_native.tif",
  dstfile = "EP_clay_native_rs.tif",
  t_srs   = "+proj=longlat +datum=WGS84 +no_defs",
  te      = c(bbox(epgsoc2)),
  tr      = res(epgsoc2),
  r       = "average",
  ot      = "Int16",
  co      = "COMPRESS=DEFLATE",
  overwrite = TRUE
)

gdalUtilities::gdalwarp(
  srcfile = "hias_clay_native.tif",
  dstfile = "hias_clay_native_rs.tif",
  t_srs   = "+proj=longlat +datum=WGS84 +no_defs",
  te      = c(bbox(wpgsoc2)),
  tr      = res(wpgsoc2),
  r       = "average",
  ot      = "Int16",
  co      = "COMPRESS=DEFLATE",
  overwrite = TRUE
)

gdalUtilities::gdalwarp(
  srcfile = "prvi_clay_native.tif",
  dstfile = "prvi_clay_native_rs.tif",
  t_srs   = "+proj=longlat +datum=WGS84 +no_defs ",
  te      = c(bbox(prvigsoc2)),
  tr      = res(prvigsoc2),
  r       = "average",
  ot      = "Int16",
  co      = "COMPRESS=DEFLATE",
  overwrite = TRUE
)



