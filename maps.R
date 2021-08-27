
# Setup ----
library(sf)
library(terra)

aoi <- "CONUS"
setwd(paste0("D:/geodata/project_data/gsp-gsocseq/", aoi))


fr_df  <- readRDS(file = "rothC_fr_final.rds")
# fr_df  <- cbind(aoi = aoi, fr_df)


# Differences, Rates, and Uncertainties,
gsocseq_maps <- within(fr_df, {
  unc_bau         = (f_t.baumax - f_t.baumin) / (2 * f_t.bau)  * 100
  # unc_t0_soc      = (SOC_t0.max - SOC_t0.min) / (2 * SOC_t0.r) * 100
  unc_t0_soc      = (SOC_t0.max - SOC_t0.min) / (2 * SOC) * 100
  unc_ssm         = (f_t.medmax - f_t.medmin) / (2 * f_t.med)  * 100
  BAU_Uncertainty = unc_bau
  T0_Uncertainty  = unc_t0_soc
  SSM_Uncertainty = unc_ssm
  T0_             = SOC
  # T0_             = SOC_t0.r
  finalSOC_BAU_   = f_t.bau
  finalSOC_SSM1_  = f_t.low
  finalSOC_SSM2_  = f_t.med
  finalSOC_SSM3_  = f_t.high
  # absolute differences (SSM - SOC 2018)
  # AbsDiff_BAU_  = f_t.bau  - SOC_t0.r
  # AbsDiff_SSM1_ = f_t.low  - SOC_t0.r
  # AbsDiff_SSM2_ = f_t.med  - SOC_t0.r
  # AbsDiff_SSM3_ = f_t.high - SOC_t0.r
  AbsDiff_BAU_  = f_t.bau  - SOC
  AbsDiff_SSM1_ = f_t.low  - SOC
  AbsDiff_SSM2_ = f_t.med  - SOC
  AbsDiff_SSM3_ = f_t.high - SOC
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
  # ASR_BAU_Uncertainty  = sqrt((unc_bau * f_t.bau)^2  + (unc_t0_soc * SOC_t0.r)^2) /
  #   abs(SOC_t0.r + f_t.bau)
  # ASR_SSM1_Uncertainty = sqrt((unc_ssm * f_t.low)^2  + (unc_t0_soc * SOC_t0.r)^2) /
  #   abs(SOC_t0.r + f_t.low)
  # ASR_SSM2_Uncertainty = sqrt((unc_ssm * f_t.med)^2  + (unc_t0_soc * SOC_t0.r)^2) /
  #   abs(SOC_t0.r + f_t.med)
  # ASR_SSM3_Uncertainty = sqrt((unc_ssm * f_t.high)^2 + (unc_t0_soc * SOC_t0.r)^2) /
  #   abs(SOC_t0.r + f_t.high)
  ASR_BAU_Uncertainty  = sqrt((unc_bau * f_t.bau)^2  + (unc_t0_soc * SOC)^2) /
    abs(SOC_t0.r + f_t.bau)
  ASR_SSM1_Uncertainty = sqrt((unc_ssm * f_t.low)^2  + (unc_t0_soc * SOC)^2) /
    abs(SOC_t0.r + f_t.low)
  ASR_SSM2_Uncertainty = sqrt((unc_ssm * f_t.med)^2  + (unc_t0_soc * SOC)^2) /
    abs(SOC_t0.r + f_t.med)
  ASR_SSM3_Uncertainty = sqrt((unc_ssm * f_t.high)^2 + (unc_t0_soc * SOC)^2) /
    abs(SOC_t0.r + f_t.high)
  # Uncertainties for the Relative difference  SSM_ - SOCBAU
  RSR_SSM1_Uncertainty = sqrt((unc_ssm * f_t.low)^2  + (unc_bau * f_t.bau)^2) /
    abs(f_t.bau + f_t.low)
  RSR_SSM2_Uncertainty = sqrt((unc_ssm * f_t.med)^2  + (unc_bau * f_t.bau)^2) /
    abs(f_t.bau + f_t.med)
  RSR_SSM3_Uncertainty = sqrt((unc_ssm * f_t.high)^2 + (unc_bau * f_t.bau)^2) /
    abs(f_t.bau + f_t.high)
})

names(gsocseq_maps) <- gsub("\\.", "_", names(gsocseq_maps))
nm    <- names(gsocseq_maps)
vars  <- c("aoi", "x", "y", "cell", "LU", "CLAY", "SOC")
vars2 <- nm[grep("^finalSOC|^T0|_Uncertainty$|^AbsDiff|^RelDiff|^ASR|^RSR", nm)]
gsocseq_maps <- gsocseq_maps[c(vars, vars2)]
saveRDS(gsocseq_maps, file = "gsocseq_maps.rds")


# Convert to points
fr_sf <- st_as_sf(
  gsocseq_maps,
  coords = c("x", "y"),
  crs    = 4326
)

ak1_fr_sf <- fr_sf[fr_sf$aoi == "AK1", ]
ak2_fr_sf <- fr_sf[fr_sf$aoi == "AK2", ]

write_sf(ak1_fr_sf, dsn = paste0("AK1_gsocseq_maps.gpkg"), driver = "GPKG", overwrite = TRUE) 
write_sf(ak2_fr_sf, dsn = paste0("AK2_gsocseq_maps.gpkg"), driver = "GPKG", overwrite = TRUE) 

write_sf(fr_sf, dsn = paste0("CONUS_gsocseq_maps.gpkg"), driver = "GPKG", overwrite = TRUE) 

aoi <- "CONUS"
gsoc <- rast(paste0(aoi, "_GSOCmap1.5.0.tif"))
gsoc[!is.na(gsoc)] <- 1
lu   <- rast(paste0(aoi, "_glc_shv10_DOM.tif"))
lu   <- lu %in% c(2, 3, 5, 12, 13)
gsoc <- gsoc * lu


# rasterize points
lapply(vars2, function(x) {
  
  f <- paste0(aoi, "_GSOCseq_", x, "Map030.tif")
  if (file.exists(f)) file.remove(f)
  
  cat("rasterizing ", f, as.character(Sys.time()), "\n")
  writeRaster(gsoc, f, overwrite = TRUE)
  
  gdalUtilities::gdal_rasterize(
    src_datasource = paste0(aoi, "_gsocseq_maps.gpkg"),
    a              = x,
    dst_filename   = f,
    of             = "GTiff",
    te             = bbox(gsoc),
    tr             = res(gsoc),
    co             = c("COMPRESS=DEFLATE"),
    a_nodata       = -999
    )
})



# QA results ----
gsocseq_maps <- readRDS("gsocseq_maps.rds")
nm    <- names(gsocseq_maps)
vars  <- c("aoi", "x", "y", "cell", "LU", "CLAY", "SOC")
vars2 <- nm[grep("^finalSOC|^T0|_Uncertainty$|^AbsDiff|^RelDiff|^ASR|^RSR", nm)]

summary(gsocseq_maps)

test <- data.frame(t(sapply(vars2, function(x) rbind(quantile(gsocseq_maps[, x], na.rm = TRUE)))))
test <- round(test, 2)
names(test) <- paste0("p", seq(0, 1, 0.25))
test <- cbind(var = row.names(test), test)
row.names(test) <- NULL
View(test)
write.csv(test, "gsocseq_summaries.csv", row.names = FALSE)


aoi <- "CONUS"
f   <- paste0(aoi, "_GSOCseq_", vars2, "Map030.tif")
rs  <- rast(f) 
names(rs) <- vars2


vars <- c(pred     = "final|T0_$", 
          # pred_unc = "^..._Uncertainty$|^.._Uncertainty",
          abs_diff = "AbsDiff", 
          rel_diff = "RelDiff",
          asr      = "ASR_SSM._$|ASR_BAU_$",
          rsr      = "RSR_SSM._$"
          # asr_unc  = "ASR_...._Uncertainty",
          # rsr_unc  = "RSR_...._Uncertainty"
          )

lapply(seq_along(vars), function(i) {
  idx <- grepl(vars[i], names(rs))
  brks <- quantile(values(rs[[idx]]), probs = seq(0, 1, 0.1), na.rm = TRUE)
  png(paste0("plots_", aoi, "_", names(vars)[i], ".png"), units = "in", width = 12, height = 6, res = 300)
  plot(rs[[idx]], breaks = brks, col = viridis::viridis(10))
  dev.off()
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

