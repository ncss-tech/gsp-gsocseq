library(sf)
library(raster)


# Set working directory
setwd("D:/geodata/project_data/gsp-gsocseq")

aoi <- "CONUS"

su_sdf <- readRDS(file = "su_sdf.RDS")
fr_df  <- readRDS(file = "rothC_fr_final_v2.rds")


# Differences, Rates, and Uncertainties,
gsocseq_maps <- within(fr_df, {
  unc_bau         = (f_t.baumax - f_t.baumin) / (2 * f_t.bau)  * 100
  unc_t0_soc      = (SOC_t0.max - SOC_t0.min) / (2 * SOC_t0.r) * 100
  unc_ssm         = (f_t.medmax - f_t.medmin) / (2 * f_t.med)  * 100
  BAU_Uncertainty = unc_bau
  T0_Uncertainty  = unc_t0_soc
  SSM_Uncertainty = unc_ssm
  T0_             = SOC_t0.r
  finalSOC_BAU_   = f_t.bau
  finalSOC_SSM1_  = f_t.low
  finalSOC_SSM2_  = f_t.med
  finalSOC_SSM3_  = f_t.high
  # absolute differences (SSM - SOC 2018)
  AbsDiff_BAU_  = f_t.bau  - SOC_t0.r
  AbsDiff_SSM1_ = f_t.low  - SOC_t0.r
  AbsDiff_SSM2_ = f_t.med  - SOC_t0.r
  AbsDiff_SSM3_ = f_t.high - SOC_t0.r
  # absolute rate
  ASR_BAU_  = AbsDiff_BAU_  / 20
  ASR_SSM1_ = AbsDiff_SSM1_ / 20
  ASR_SSM2_ = AbsDiff_SSM2_ / 20
  ASR_SSM3_ = AbsDiff_SSM3_ / 20
  # relative differences (SSM - SOC BAU)
  RelDiff_SSM1_ = f_t.low  - f_t.bau
  RelDiff_SSM2_ = f_t.med  - f_t.bau
  RelDiff_SSM3_ = f_t.high - f_t.bau
  # relative rate
  RSR_SSM1_     = RelDiff_SSM1_ / 20
  RSR_SSM2_     = RelDiff_SSM2_ / 20
  RSR_SSM3_     = RelDiff_SSM3_ / 20
  # Uncertainties for the Absolute difference SSM_ - SOC2018
  ASR_BAU_Uncertainty  = sqrt((unc_bau * f_t.bau)^2  + (unc_t0_soc * SOC_t0.r)^2) / 
    abs(SOC_t0.r + f_t.bau)
  ASR_SSM1_Uncertainty = sqrt((unc_ssm * f_t.low)^2  + (unc_t0_soc * SOC_t0.r)^2) / 
    abs(SOC_t0.r + f_t.low)
  ASR_SSM2_Uncertainty = sqrt((unc_ssm * f_t.med)^2  + (unc_t0_soc * SOC_t0.r)^2) / 
    abs(SOC_t0.r + f_t.med)
  ASR_SSM3_Uncertainty = sqrt((unc_ssm * f_t.high)^2 + (unc_t0_soc * SOC_t0.r)^2) / 
    abs(SOC_t0.r + f_t.high)
  # Uncertainties for the Relative difference  SSM_ - SOCBAU
  RSR_SSM1_Uncertainty = sqrt((unc_ssm * f_t.low)^2  + (unc_bau * f_t.bau)^2) / 
    abs(f_t.bau + f_t.low)
  RSR_SSM2_Uncertainty = sqrt((unc_ssm * f_t.med)^2  + (unc_bau * f_t.bau)^2) / 
    abs(f_t.bau + f_t.med)
  RSR_SSM3_Uncertainty = sqrt((unc_ssm * f_t.high)^2 + (unc_bau * f_t.bau)^2) / 
    abs(f_t.bau + f_t.high)
})

vars  <- c("id", "x", "y", "LU", "CLAY", "SOC")
names(gsocseq_maps) <- gsub("\\.", "_", names(gsocseq_maps))
nm    <- names(gsocseq_maps)
vars2 <- nm[grep("^finalSOC|^T0|_Uncertainty$|^AbsDiff|^RelDiff|^ASR|^RSR", nm)]
gsocseq_maps <- gsocseq_maps[c(vars, vars2)]



# Convert to points
fr_sf <- st_as_sf(
  gsocseq_maps,
  coords = c("x", "y"),
  crs    = 4326
)
write_sf(fr_sf, dsn = "gsocseq_maps.gpkg", driver = "GPKG", overwrite = TRUE) 
test <- read_sf(dsn = "gsocseq_maps.gpkg", driver = "GPKG")


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



# Create map ----
aoi   <- read_sf(dsn = "AOIs/CONUS.shp")
gsoc  <- rast("D:/geodata/soils/GSOCmap1.5.0.tif")
gsoc2 <- crop(gsoc, aoi)
gsoc2 <- trim(gsoc2)

gdalUtilities::gdal_rasterize(
  src_datasource = "AOIs/CONUS.shp",
  dst_filename   = "CONUS.tif",
  of             = "GTiff",
  te             = bbox(gsoc2),
  tr             = res(gsoc2),
  a_nodata       = -99999
)
conus_r <- rast("CONUS.tif")
conus_r[!is.na(conus_r)] <- 1
gsoc3 <- gsoc2 * conus_r
gsoc3 <- project(gsoc3, "epsg:5070")
writeRaster(gsoc3, "D:/geodata/soils/GSOCmap1.5.0_aea.tif")
quantile(values(gsoc3), p = seq(0, 1, 0.2), na.rm = TRUE)
brks <- c(0, 20, 30, 40, 60, 720)


bau  <- raster("CONUS_fr_bau.tif")
bau  <- mask(bau, soc)
plot(bau, breaks = brks, lab.breaks = brks, col = viridis::viridis(n = 4))


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

