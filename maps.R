library(sf)
library(raster)


# Set working directory
setwd("D:/geodata/project_data/gsp-gsocseq")

aoi <- "CONUS"

su_sdf <- readRDS(file = "su_sdf.RDS")
su_sf <- st_as_sf(su_sdf)
write_sf(wu_sf, dsn = "test.gpkg", driver = "GPKG")

gsoc   <- raster("CONUS_GSOCmap1.5.0.tif")
fr_sf <- read_sf(dsn = "rothC_fr_final.gpkg") 


vars <- c("bau", "baumin", "baumax", "low", "med", "high", "medmin", "medmax")

lapply(vars[2:8], function(x) {
  cat("rasterizing ", x, as.character(Sys.time()), "\n")
  writeRaster(gsoc, paste0(aoi, "_fr_", x, ".tif"), overwrite = TRUE)
  test <- gdalUtilities::gdal_rasterize(
    src_datasource = "rothC_fr_final_v2.gpkg",
    a              = paste0("f_t_", x),
    dst_filename   = paste0(aoi, "_fr_", x, ".tif"),
    of             = "GTiff",
    te             = bbox(gsoc),
    tr             = res(gsoc),
    ot             = -99999
    )
})


# Calculate differences ----

# Difference BAU 2040 - SOC 2018
bau_soc_dif <- fr_bau - gsoc
writeRaster(
  bau_soc_dif,
  filename = paste0(aoi, "_GSOCseq_AbsDiff_BAU_Map030"),
  format = "GTiff"
)
writeRaster(
  bau_soc_dif / 20,
  filename = paste0(aoi, "_GSOCseq_ASR_BAU_Map030"),
  format = "GTiff"
)


# Difference Low Scenario - SOC 2018
Diff_Lw_SOC_2018 <- Country_Lwr_2040_Map - Country_SOC_2018_Map
writeRaster(
  Diff_Lw_SOC_2018,
  filename = paste0(aoi, "_GSOCseq_AbsDiff_SSM1_Map030"),
  format = "GTiff"
)
writeRaster(
  Diff_Lw_SOC_2018 / 20,
  filename = paste0(aoi, "_GSOCseq_ASR_SSM1_Map030"),
  format = "GTiff"
)


# Difference Med Scenario - SOC 2018
Diff_Md_SOC_2018 <- Country_Med_2040_Map - Country_SOC_2018_Map
writeRaster(
  Diff_Md_SOC_2018,
  filename = paste0(aoi, "_GSOCseq_AbsDiff_SSM2_Map030"),
  format = "GTiff"
)
writeRaster(
  Diff_Md_SOC_2018 / 20,
  filename = paste0(aoi, "_GSOCseq_ASR_SSM2_Map030"),
  format = "GTiff"
)


# Difference High Scenario - SOC 2018
Diff_Hg_SOC_2018 <- Country_Hgh_2040_Map - Country_SOC_2018_Map
writeRaster(
  Diff_Hg_SOC_2018,
  filename = paste0(aoi, "_GSOCseq_AbsDiff_SSM3_Map030"),
  format = "GTiff"
)
writeRaster(
  Diff_Hg_SOC_2018 / 20,
  filename = paste0(aoi, "_GSOCseq_ASR_SSM3_Map030"),
  format = "GTiff"
)



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

