library(sf)
library(raster)


# Set working directory
setwd("D:/geodata/project_data/gsp-gsocseq")


su_sdf <- readRDS(file = "su_sdf.RDS")
su_sf <- st_as_sf(su_sdf)
write_sf(wu_sf, dsn = "test.gpkg", driver = "GPKG")

soc   <- raster("CONUS_GSOCmap1.5.0.tif")
fr_sf <- read_sf(dsn = "rothC_fr_final.gpkg") 

writeRaster(soc, "CONUS_fr_bau.tif")

test <- gdalUtilities::gdal_rasterize(
  src_datasource = "rothC_fr_final.gpkg",
  a              = "f_t_bau",
  dst_filename   = "CONUS_fr_bau.tif",
  of             = "GTiff",
  te             = bbox(soc),
  tr             = res(soc),
  ot             = -99999
)

bau  <- raster("CONUS_fr_bau.tif")
bau  <- mask(bau, soc)
brks <- c(0, 15, 30, 45, 790)
plot(bau, breaks = brks, lab.breaks = brks, col = viridis::viridis(n = 4))

tm_shape(bau) + 
  tm_raster(
    breaks = brks,
    palette = viridis::viridis_pal()(4)
  ) +
  tm_legend(
    legend.outside = TRUE, 
    legend.outside.position = c("right", "top")
  ) +
  tm_layout(
    main.title = paste("Carbon Sequestration T/ha/yr") #, var) #, 
  )
