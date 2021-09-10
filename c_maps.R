rm(list = ls())

library(sf)
library(raster)
library(dplyr)
library(terra)
library(rgdal)


# Set working directory
# setwd("D:/geodata/project_data/gsp-gsocseq")
 
# setwd("D:/geodata/project_data/gsp-gsocseq")
# aoi <- "EPAC"
# fr_df <- readRDS(file = "ep_rothC_fr_final.rds")

setwd("D:/geodata/project_data/gsp-gsocseqwp")
aoi <- "HIAS"
fr_df <- readRDS(file = "HIAS_rothC_fr_final.rds")

# setwd("D:/geodata/project_data/gsp-gsocseqp")
# aoi <- 'PRVI'
# fr_df <- readRDS(file = "prvi_rothC_fr_final.rds")

su_sdf <- readRDS(file = "su_sf.RDS")

# fr_df <- fr_df %>%
#   mutate(x = unlist(map(fr_df$geometry,1)), y = unlist(map(fr_df$geometry,2)))

# Differences, Rates, and Uncertainties,
gsocseq_maps <- within(fr_df, {
  unc_bau         = (f_t_baumax - f_t_baumin) / (2 * f_t_bau)  * 100
  unc_t0_soc      = (SOC_t0_max - SOC_t0_min) / (2 * SOC_t0_r) * 100
  unc_ssm         = (f_t_medmax - f_t_medmin) / (2 * f_t_med)  * 100
  BAU_Uncertainty = unc_bau
  T0_Uncertainty  = unc_t0_soc
  SSM_Uncertainty = unc_ssm
  T0_             = SOC_t0_r
  finalSOC_BAU_   = f_t_bau
  finalSOC_SSM1_  = f_t_low
  finalSOC_SSM2_  = f_t_med
  finalSOC_SSM3_  = f_t_high
  
  # absolute differences (SSM - SOC 2018)
  AbsDiff_BAU_  = f_t_bau  - SOC_t0_r
  AbsDiff_SSM1_ = f_t_low  - SOC_t0_r
  AbsDiff_SSM2_ = f_t_med  - SOC_t0_r
  AbsDiff_SSM3_ = f_t_high - SOC_t0_r
  
  # absolute rate
  ASR_BAU_  = AbsDiff_BAU_  / 20
  ASR_SSM1_ = AbsDiff_SSM1_ / 20
  ASR_SSM2_ = AbsDiff_SSM2_ / 20
  ASR_SSM3_ = AbsDiff_SSM3_ / 20
  
  # relative differences (SSM - SOC BAU)
  RelDiff_SSM1_ = f_t_low  - f_t_bau
  RelDiff_SSM2_ = f_t_med  - f_t_bau
  RelDiff_SSM3_ = f_t_high - f_t_bau
  
  # relative rate
  RSR_SSM1_     = RelDiff_SSM1_ / 20
  RSR_SSM2_     = RelDiff_SSM2_ / 20
  RSR_SSM3_     = RelDiff_SSM3_ / 20
 
  # Uncertainties for the Absolute difference SSM_ - SOC2018
  ASR_BAU_Uncertainty  = sqrt((unc_bau * f_t_bau)^2  + (unc_t0_soc * SOC_t0_r)^2) / 
    abs(SOC_t0_r + f_t_bau)
  ASR_SSM1_Uncertainty = sqrt((unc_ssm * f_t_low)^2  + (unc_t0_soc * SOC_t0_r)^2) / 
    abs(SOC_t0_r + f_t_low)
  ASR_SSM2_Uncertainty = sqrt((unc_ssm * f_t_med)^2  + (unc_t0_soc * SOC_t0_r)^2) / 
    abs(SOC_t0_r + f_t_med)
  ASR_SSM3_Uncertainty = sqrt((unc_ssm * f_t_high)^2 + (unc_t0_soc * SOC_t0_r)^2) / 
    abs(SOC_t0_r + f_t_high)
  
  
  # Uncertainties for the Relative difference  SSM_ - SOCBAU
  RSR_SSM1_Uncertainty = sqrt((unc_ssm * f_t_low)^2  + (unc_bau * f_t_bau)^2) / 
    abs(f_t_bau + f_t_low)
  RSR_SSM2_Uncertainty = sqrt((unc_ssm * f_t_med)^2  + (unc_bau * f_t_bau)^2) / 
    abs(f_t_bau + f_t_med)
  RSR_SSM3_Uncertainty = sqrt((unc_ssm * f_t_high)^2 + (unc_bau * f_t_bau)^2) / 
    abs(f_t_bau + f_t_high)
})

gsocseq_maps <- cbind(gsocseq_maps, st_coordinates(gsocseq_maps))

vars  <- c("id", "X", "Y", "LU", "CLAY", "SOC")
names(gsocseq_maps) <- gsub("\\.", "_", names(gsocseq_maps))
nm    <- names(gsocseq_maps)
vars2 <- nm[grep("^finalSOC|^T0|_Uncertainty$|^AbsDiff|^RelDiff|^ASR|^RSR", nm)]
gsocseq_maps <- gsocseq_maps[c(vars, vars2)]


gsocseq_maps <- gsocseq_maps[gsocseq_maps$finalSOC_BAU_ < 801, ]
gsocseq_maps <- gsocseq_maps[gsocseq_maps$finalSOC_BAU_ > 0, ]
gsocseq_maps <- gsocseq_maps[!is.na(gsocseq_maps$finalSOC_BAU_), ]
gsocseq_maps <- gsocseq_maps[gsocseq_maps$LU %in% c(2, 3, 5, 12, 13), ]
# gsocseq_maps <- gsocseq_maps[ ,!names(gsocseq_maps) %in% c("X", "Y")]


# Convert to points
fr_sf <- st_as_sf(
  gsocseq_maps,
  coords = c("X", "Y"),
  crs    = 4326
)
write_sf(fr_sf, dsn = "gsocseq_maps.gpkg", sf_column_name = "geometry",  driver = "GPKG", overwrite = TRUE)
# test <- read_sf(dsn = "gsocseq_maps.gpkg",  geometry_column = c("geometry"), driver = "GPKG")
test <- readOGR(dsn = "gsocseq_maps.gpkg")

gsoc <- raster("CONUS_GSOCmap1.5.0.tif")
gsoc[!is.na(gsoc)] <- 1
lu   <- raster("CONUS_glc_shv10_DOM.tif")
lu   <- lu %in% c(2, 3, 5, 12, 13)
gsoc <- gsoc * lu


# rasterize points
lapply(vars2, function(x) {
  
  f <- paste0(aoi, "_GSOCseq_", x, "Map030.tif")
  cat("rasterizing ", f, as.character(Sys.time()), "\n")
  writeRaster(gsoc, f, overwrite = TRUE)
  
  test <- gdalUtilities::gdal_rasterize(
    src_datasource = "gsocseq_maps.gpkg",
    a              = x,
    dst_filename   = f,
    of             = "GTiff",
    te             = bbox(gsoc),
    tr             = res(gsoc),
    co             = c("TILED=YES", "COMPRESS=DEFLATE"),
    a_nodata       = -99999
    )
})



# Create map, mask input SSURGO carbon to identified LU ----
lu   <- raster("CONUS_glc_shv10_DOM.tif")
lu   <- lu %in% c(2, 3, 5, 12, 13)
ssurgo_c <- raster("CONUS_GSOCmap1.5.0.tif")
lu_c <- ssurgo_c * lu
lu_c[lu_c == 0] <- NA
writeRaster(lu_c, "SSURGO_Carbon.tif", overwrite = TRUE)


# gsoc <- raster("CONUS_GSOCmap1.5.0.tif")
# gsoc[!is.na(gsoc)] <- 1
# lu   <- raster("CONUS_glc_shv10_DOM.tif")
# 
# # binary response
# lu   <- lu %in% c(2, 3, 5, 12, 13)
# 
# gsoc <- gsoc * lu


# aoi   <- read_sf(dsn = "AOIs/AOI_CONUS_bufbox.shp")
# gsoc  <- rast("CONUS_GSOCmap1.5.0.tif")
# gsoc2 <- crop(gsoc, aoi)
# gsoc2 <- trim(gsoc2)
#

# rasterize AOI
# gdalUtilities::gdal_rasterize(
#   src_datasource = "AOIs/AOI_CONUS_bufbox.shp",
#   dst_filename   = "AOI_RASTER.tif",
#   of             = "GTiff",
#   te             = bbox(gsoc2),
#   tr             = res(gsoc2),
#   a_nodata       = -99999
# )

# multiply AOI (value 1) with ssurgo carbon and project to 5070
conus_r <- rast("AOI_RASTER.tif")
conus_r[!is.na(conus_r)] <- 1
gsoc3 <- gsoc2 * conus_r
# gsoc3 <- project(gsoc3, "epsg:5070")
# writeRaster(gsoc3, "D:/geodata/soils/GSOCmap1.5.0_aea.tif")
# quantile(values(gsoc3), p = seq(0, 1, 0.2), na.rm = TRUE)
# brks <- c(0, 20, 30, 40, 60, 720)


# bau  <- raster("CONUS_fr_bau.tif")
bau  <- raster("PRVI_GSOCseq_finalSOC_BAU_Map030.tif")
bau  <- mask(bau, gsoc)
brks <- c(0, 20, 30, 40, 60, 720)
plot(bau, breaks = brks, lab.breaks = brks, col = viridis::viridis(n = 4))

library(mapview)
mapview(bau)


tm_shape(gsoc3,
         raster.downsample = FALSE
         ) + 
  tm_raster(
    breaks = brks,
    palette = RColorBrewer::brewer.pal(5, "Greys"),
    
  ) +
  tm_legend(
    legend.outside = TRUE, 
    legend.outside.position = c("right", "top")
  ) +
  tm_layout(
    main.title = paste("Soil Organic Carbon (T/ha)") #, var) #, 
  )

