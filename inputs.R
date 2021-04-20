library(raster)


# TerraClimate ----

tmp8100 <- stack("CONUS_AverageTemperature_1981-2001.tif")
tmp0121 <- stack("CONUS_AverageTemperature_2001-2021.tif")
ppt8100 <- stack("CONUS_Precipitation_1981-2001.tif")
ppt0121 <- stack("CONUS_Precipitation_2001-2021.tif")
pet8100 <- stack("CONUS_PET_1981-2001.tif")
pet0121 <- stack("CONUS_PET_2001-2021.tif")


mos <- formatC(1:12, width = 2, flag = "0")
mos8100 <- lapply(mos, function(x) paste0("X", 1981:2000, x))
mos0121 <- lapply(mos, function(x) paste0("X", 2001:2021, x))

avgRS <- function(rs, var, mos) { 
  rs_l <- lapply(mos, function(x) {
    vars <- paste0(x, "_", var)
    idx  <- which(names(rs) %in% vars)
    cat("averaging", paste(names(rs[[idx]]), collapse = ", "), "\n", "\n")
    mean(rs[[idx]], na.rm = TRUE)
    })
  final_rs <- stack(rs_l)
  names(final_rs) <- month.abb
  return(final_rs)
}

tmp8100_avg <- avgRS(tmp8100, "tmmx", mos)
tmp8100_avg <- tmp8100_avg * 0.1
tmp0121_avg <- avgRS(tmp0121, var = "tmmx", mos0121)
tmp0121_avg <- tmp0121_avg * 0.1
writeRaster(tmp8100_avg, filename = "CONUS_Temp_Stack_81-00_TC.tif")
writeRaster(tmp0121_avg, filename = "CONUS_Temp_Stack_01-21_TC.tif")

ppt8100_avg <- avgRS(ppt8100, "pr", mos8100)
ppt0121_avg <- avgRS(ppt0121, "pr", mos0121)
writeRaster(ppt8100_avg, filename = "CONUS_Prec_Stack_81-00_TC.tif")
writeRaster(ppt0121_avg, filename = "CONUS_Prec_Stack_01-21_TC.tif")

pet8100_avg <- avgRS(pet8100, "pet", mos8100)
pet0121_avg <- avgRS(pet0121, "pet", mos0121)
writeRaster(pet8100_avg, filename = "CONUS_PET_Stack_81-00_TC.tif")
writeRaster(pet0121_avg, filename = "CONUS_PET_Stack_01-21_TC.tif")


yrs8100 <- lapply(1981:2000, function(x) paste0("X", x))
yrs8100 <- lapply(yrs8100, function(x) paste0(x, formatC(1:12, width = 2, flag = "0")))
yrs0121 <- lapply(2001:2021, function(x) paste0("X", x))
yrs0121 <- lapply(yrs0121, function(x) paste0(x, formatC(1:12, width = 2, flag = "0")))


avgRSyr <- function(rs, var, yrs) { 
  rs_l <- lapply(yrs, function(x) {
    vars <- paste0(x, "_", var)
    idx  <- which(names(rs) %in% vars)
    cat("averaging", paste(names(rs[[idx]]), collapse = ", "), "\n", "\n")
    mean(rs[[idx]], na.rm = TRUE)
  })
  final_rs <- stack(rs_l)
  # names(final_rs) <- paste0("X", substr(x[[1]][1], 2, 5))
  return(final_rs)
}

tmp8100_avg_yr <- avgRSyr(tmp8100, "tmmx", yrs8100)
tmp8100_avg_yr <- tmp8100_avg_yr * 0.1
tmp0121_avg_yr <- avgRSyr(tmp0121, "tmmx", yrs0121)
tmp0121_avg_yr <- tmp0121_avg_yr * 0.1


# landcover ----
lc <- raster("D:/geodata/land_use_land_cover/GlcShare_v10_Dominant/glc_shv10_DOM.Tif")
projection(lc) <- "+init=epsg:4326"

lc2 <- crop(lc, ssurgo_st, file = "glc_shv10_DOM_CONUS.tif", progress = "text")



# GSOC ----
gsoc <- raster("D:/geodata/soils/GSOCmap1.5.0.tif")
gsoc2 <- crop(gsoc, ssurgo_st)
gsoc2 <- resample(gsoc2, ssurgo_st, method = "ngb", filename =  "GSOCmap1.5.0_CONUS.tif", progress = "text")



cms_soc <- raster("D:/geodata/soils/CMS_SOC_Mexico_CONUS_1737/data/SOC_prediction_1991_2010.tif")
cms_soc_1km <- aggregate(cms_soc, fact = 4, fun = mean, progress = "text")
cms_soc_1km <- readAll(cms_soc_1km)
cms_soc_1km <- projectRaster(from = cms_soc_1km, to = ssurgo_st, filename = "cms_soc_1km.tif", method = "bilinear", progress = "text", overwrite = TRUE)


