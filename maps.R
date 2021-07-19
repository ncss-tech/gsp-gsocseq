library(sf)
library(raster)


# Set working directory
setwd("D:/geodata/project_data/gsp-gsocseq")

aoi <- "CONUS"

su_sdf <- readRDS(file = "su_sdf.RDS")
fr_df  <- readRDS(file = "rothC_fr_final_v2.rds")


# Differences, Rates, and Uncertainties,
gsocseq_maps <- within(fr_df, {
  unc_bau     = (f_t.baumax - f_t.baumin) / (2 * f_t.bau)  * 100
  unc_wu_soc  = (SOC_t0.max - SOC_t0.min) / (2 * SOC_t0.r) * 100
  unc_ssm     = (f_t.medmax - f_t.medmin) / (2 * f_t.med)  * 100
  # absolute differences (SSM - SOC 2018)
  bau_soc_dif  = f_t.bau  - SOC_t0.r
  low_soc_dif  = f_t.low  - SOC_t0.r
  med_soc_dif  = f_t.med  - SOC_t0.r
  high_soc_dif = f_t.high - SOC_t0.r
  # absolute rate
  bau_soc_rate  = bau_soc_dif / 20
  low_soc_rate  = low_soc_dif / 20
  med_soc_rate  = med_soc_dif / 20
  high_soc_rate = high_soc_dif / 20
  # relative differences (SSM - SOC BAU)
  low_bau_dif  = f_t.low  - f_t.bau
  med_bau_dif  = f_t.med  - f_t.bau
  high_bau_dif = f_t.high - f_t.bau
  # relative rate
  low_bau_rate  = low_bau_dif / 20
  med_bau_rate  = med_bau_dif / 20
  high_bau_rate = high_bau_dif / 20
  # Uncertainties for the Absolute difference SSM_ - SOC2018
  unc_abs_bau = sqrt((unc_bau  * f_t.bau)^2  + (unc_wu_soc * SOC_t0.r)^2) / 
    abs(SOC_t0.r + f_t.bau)
  unc_abs_low  = sqrt((unc_ssm * f_t.low)^2  + (unc_wu_soc * SOC_t0.r)^2) / 
    abs(SOC_t0.r + f_t.low)
  unc_abs_med  = sqrt((unc_ssm * f_t.med)^2  + (unc_wu_soc * SOC_t0.r)^2) / 
    abs(SOC_t0.r + f_t.med)
  unc_abs_high = sqrt((unc_ssm * f_t.high)^2 + (unc_wu_soc * SOC_t0.r)^2) / 
    abs(SOC_t0.r + f_t.high)
  # Uncertainties for the Relative difference  SSM_ - SOCBAU
  unc_rel_low  = sqrt((unc_ssm * f_t.low)^2  + (unc_bau * f_t.bau)^2) / 
    abs(f_t.bau + f_t.low)
  unc_rel_mod  = sqrt((unc_ssm * f_t.med)^2  + (unc_bau * f_t.bau)^2) / 
    abs(f_t.bau + f_t.med)
  unc_rel_high = sqrt((unc_ssm * f_t.high)^2 + (unc_bau * f_t.bau)^2) / 
    abs(f_t.bau + f_t.high)
})

vars  <- c("id", "x", "y", "LU", "CLAY", "SOC")
names(gsocseq_maps) <- gsub("\\.", "_", names(gsocseq_maps))
nm    <- names(gsocseq_maps)
vars2 <- nm[grep("^f_t|^unc_|_dif$|_rate$", nm)]
gsocseq_maps <- gsocseq_maps[c(vars, vars2)]



# Convert to points
fr_sf <- st_as_sf(
  gsocseq_maps,
  coords = c("x", "y"),
  crs    = 4326
)
write_sf(fr_sf, dsn = "gsocseq_maps.gpkg", driver = "GPKG", overwrite = TRUE) 


gsoc <- raster("CONUS_GSOCmap1.5.0.tif")
gsoc[!is.na(gsoc)] <- 1
lu   <- raster("CONUS_glc_shv10_DOM.tif")
lu   <- lu %in% c(2, 3, 5, 12, 13)
gsoc <- gsoc * lu


# rasterize points
lapply(vars2, function(x) {
  cat("rasterizing ", x, as.character(Sys.time()), "\n")
  writeRaster(gsoc, paste0(aoi, "_", x, ".tif"), overwrite = TRUE)
  test <- gdalUtilities::gdal_rasterize(
    src_datasource = "gsocseq_maps.gpkg",
    a              = x,
    dst_filename   = paste0(aoi, "_", x, ".tif"),
    of             = "GTiff",
    te             = bbox(gsoc),
    tr             = res(gsoc),
    a_nodata       = -99999
    )
})



# Create map ----
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

