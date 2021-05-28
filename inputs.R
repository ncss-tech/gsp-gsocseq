library(raster)
library(sf)
library(rmapshaper)

setwd("D:/geodata/project_data/gsp-gsocseq")


# AOI ----
# CONUS
conus_sa <- read_sf(dsn = "D:/geodata/soils/gSSURGO_CONUS_FY20_july.gdb", layer = "SAPOLYGON")
conus_sa <- conus_sa %>%
  ms_simplify() %>%
  mutate(asym = "CONUS") %>%
  ms_dissolve(field = "asym")
write_sf(conus_sa, dsn = "CONUS.shp")

aoi_conus_bufbox <- conus_sa %>%
  mutate(asym = "CONUS") %>%
  ms_dissolve(field = "asym") %>%
  st_buffer(dist = 1000 * 100) %>%
  st_transform(crs = 4326) %>%
  st_bbox() %>%
  st_as_sfc()
mapview::mapview(aoi_conus_bufbox)
write_sf(aoi_conus_bufbox, dsn = "AOI_CONUS_bufbox.shp", delete_dsn = TRUE)


# OCONUS
asym <- c("AS", "AK", "FM", "GU", "HI", "MH", "MP", "PR", "PW", "VI")
leg <- lapply(asym, function(x) get_legend_from_SDA(WHERE = paste0("areasymbol LIKE '", x, "%'")))
leg <- do.call("rbind", leg)
leg$asym <- substr(leg$areasymbol, 1, 2)

sapol <- {
  split(leg, leg$asym) ->.;
  lapply(., function(x) {
    cat("fetching ", x$areasymbol, "\n")
    temp = fetchSDA_spatial(x$areasymbol, geom.src = "sapolygon")
    temp = sf::st_as_sf(temp)
  }) ->.;
}
# saveRDS(sapol, file = "sapol_oconus.rds")
sapol <- readRDS("sapol_oconus.rds")

sapol_bndy <- lapply(sapol, function(x) {
  x$asym = substr(x$areasymbol, 1, 2)
  # if (asym == "AK") {
  #   split(x, x$areasymbol) ->.;
  #   lapply(., function(x) {
  #     test = as.data.frame(st_coordinates(x))
  #     x$test = any(test$X > 0)
  #     return(x)
  #     }) ->.;
  #   do.call("rbind", .) ->.;
  # }
  x$asym = ifelse(x$areasymbol %in% paste0("AK", 783:788), "AK2", x$asym)
  temp = ms_simplify(x) %>%
    ms_dissolve(field = "asym") %>%
    ms_explode()
  return(temp)
})

sapol_bufbox <- {
  do.call("rbind", sapol_bndy) ->.;
  split(., .$asym) ->.;
  lapply(., function(x) {
    temp <- x %>%
      st_transform(crs = 6933) %>%
      st_buffer(dist = 1000 * if(all(x$asym == "AK2")) 10 else 10) %>%
      st_transform(crs = 4326) %>%
      st_bbox() %>%
      st_as_sfc() %>%
      st_as_sf()
    temp$asym <- unique(x$asym)
    return(temp)
  }) -> .;
}
test <- sapol_bufbox$PR
mapview::mapview(test)

lapply(sapol_bufbox, function(x) write_sf(x, dsn = paste0("AOI_OCONUS_", x$asym, "_bufbox.shp")))
lapply(sapol_bndy, function(x)   write_sf(x, dsn = paste0("OCONUS_", x$asym, "_boundry.shp")))



# GSOC ----
aoi   <- read_sf(dsn = "AOIs/AOI_CONUS_bufbox.shp")
gsoc  <- raster("D:/geodata/soils/GSOCmap1.5.0.tif")
gsoc2 <- crop(gsoc, aoi)
# writeRaster(gsoc2, filename =  "CONUS_GSOCmap1.5.0.tif", progress = "text", overwrite = TRUE)
gsoc2 <- raster("CONUS_GSOCmap1.5.0.tif")

cms_soc <- raster("D:/geodata/soils/CMS_SOC_Mexico_CONUS_1737/data/SOC_prediction_1991_2010.tif")
cms_soc_1km <- aggregate(cms_soc, fact = 4, fun = mean, progress = "text")
cms_soc_1km <- readAll(cms_soc_1km)
cms_soc_1km <- projectRaster(from = cms_soc_1km, to = ssurgo_st, filename = "cms_soc_1km.tif", method = "bilinear", progress = "text", overwrite = TRUE)



# clay ----
clay <- readAll(raster("gnatsgo_fy20_1km_clay_wt.tif"))
clay_4326 <- projectRaster(
  clay, 
  gsoc2, 
  filename = "CONUS_gnatsgo_fy20_1km_clay_wt.tif", 
  method = "bilinear", 
  progress = "text", 
  overwrite = TRUE, 
  datatype = "INT4S"
)




# TerraClimate ----

# load data
tmp8100 <- stack("CONUS_AverageTemperature_1981-2001.tif")
tmp0121 <- stack("CONUS_AverageTemperature_2001-2021.tif")
ppt8100 <- stack("CONUS_Precipitation_1981-2001.tif")
ppt0121 <- stack("CONUS_Precipitation_2001-2021.tif")
pet8100 <- stack("CONUS_PET_1981-2001.tif")
pet0121 <- stack("CONUS_PET_2001-2021.tif")


# create list of raster names by month
mos <- formatC(1:12, width = 2, flag = "0")
mos8100 <- lapply(mos, function(x) paste0("X", 1981:2000, x))
mos0119 <- lapply(mos, function(x) paste0("X", 2001:2019, x))


# function to average months
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


# average months
tmp8100_avg <- avgRS(tmp8100, "tmmx", mos8100)
tmp8100_avg <- tmp8100_avg * 0.1
tmp0121_avg <- avgRS(tmp0121, var = "tmmx", mos0121)
tmp0121_avg <- tmp0121_avg * 0.1
ppt8100_avg <- avgRS(ppt8100, "pr", mos8100)
ppt0121_avg <- avgRS(ppt0121, "pr", mos0119)
pet8100_avg <- avgRS(pet8100, "pet", mos8100)
pet0121_avg <- avgRS(pet0121, "pet", mos0119)


# resample to GSOC extent
tmp8100_avg <- resample(readAll(tmp8100_avg), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
tmp0121_avg <- resample(readAll(tmp0121_avg), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
ppt8100_avg <- resample(readAll(ppt8100_avg), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
ppt0121_avg <- resample(readAll(ppt0121_avg), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
pet8100_avg <- resample(readAll(pet8100_avg), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
pet0121_avg <- resample(readAll(pet0121_avg), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")


# cache files
# writeRaster(tmp8100_avg, filename = "CONUS_Temp_Stack_81-00_TC.tif", overwrite = TRUE, datatype = "INT2S", options = c("COMPRESS=DEFLATE"))
# writeRaster(tmp0121_avg, filename = "CONUS_Temp_Stack_01-19_TC.tif", overwrite = TRUE, datatype = "INT2S", options = c("COMPRESS=DEFLATE"))
# writeRaster(ppt8100_avg, filename = "CONUS_Prec_Stack_81-00_TC.tif", overwrite = TRUE, datatype = "INT2S", options = c("COMPRESS=DEFLATE"))
# writeRaster(ppt0121_avg, filename = "CONUS_Prec_Stack_01-19_TC.tif", overwrite = TRUE, datatype = "INT2S", options = c("COMPRESS=DEFLATE"))
# writeRaster(pet8100_avg, filename = "CONUS_PET_Stack_81-00_TC.tif", overwrite = TRUE, options = c("COMPRESS=DEFLATE"))
# writeRaster(pet0121_avg, filename = "CONUS_PET_Stack_01-19_TC.tif", overwrite = TRUE, options = c("COMPRESS=DEFLATE"))


# Resample original data to GSOC extent ----
# tmp8100
gdalUtilities::gdalwarp(
  srcfile = "CONUS_AverageTemperature_1981-2001.tif",
  dstfile = "CONUS_Temp_Stack_228_81_00_TC.tif",
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  ot      = "Int16",
  co      = "COMPRESS=DEFLATE"
)
# tmp0119
gdalUtilities::gdalwarp(
  srcfile = "CONUS_AverageTemperature_2001-2021.tif",
  dstfile = "CONUS_Temp_Stack_228_01_19_TC.tif",
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  ot      = "Int16",
  co      = "COMPRESS=DEFLATE"
)
# ppt8100
gdalUtilities::gdalwarp(
  srcfile = "CONUS_Precipitation_1981-2001.tif",
  dstfile = "CONUS_Prec_Stack_228_81_00_TC.tif",
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  ot      = "Int16",
  co      = "COMPRESS=DEFLATE"
)
# ppt0119
gdalUtilities::gdalwarp(
  srcfile = "CONUS_Precipitation_2001-2021.tif",
  dstfile = "CONUS_Prec_Stack_228_01_19_TC.tif",
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  ot      = "Int16",
  co      = "COMPRESS=DEFLATE"
)
# pet8100
gdalUtilities::gdalwarp(
  srcfile = "CONUS_PET_1981-2001.tif",
  dstfile = "CONUS_PET_Stack_228_81-00_TC.tif",
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  ot      = "Int16",
  co      = "COMPRESS=DEFLATE"
)
# pet0119
gdalUtilities::gdalwarp(
  srcfile = "CONUS_PET_2001-2021.tif",
  dstfile = "CONUS_PET_Stack_228_01-19_TC.tif",
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  ot      = "Int16",
  co      = "COMPRESS=DEFLATE"
)


# load cached files
tmp8100_avg <- stack("CONUS_Temp_Stack_81-00_TC.tif")
tmp0121_avg <- stack("CONUS_Temp_Stack_01-19_TC.tif")
ppt8100_avg <- stack("CONUS_Prec_Stack_81-00_TC.tif")
ppt0121_avg <- stack("CONUS_Prec_Stack_01-19_TC.tif")
pet8100_avg <- stack("CONUS_PET_Stack_81-00_TC.tif")
pet0121_avg <- stack("CONUS_PET_Stack_01-19_TC.tif")


# create list of raster names by year
yrs8100 <- lapply(1981:2000, function(x) paste0("X", x))
yrs8100 <- lapply(yrs8100, function(x) paste0(x, formatC(1:12, width = 2, flag = "0")))
yrs0121 <- lapply(2001:2019, function(x) paste0("X", x))
yrs0121 <- lapply(yrs0121, function(x) paste0(x, formatC(1:12, width = 2, flag = "0")))


# function to average years
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


# average by year
tmp8100_avg_yr <- avgRSyr(tmp8100, "tmmx", yrs8100)
tmp8100_avg_yr <- tmp8100_avg_yr * 0.1
tmp0121_avg_yr <- avgRSyr(tmp0121, "tmmx", yrs0121)
tmp0121_avg_yr <- tmp0121_avg_yr * 0.1
ppt8100_avg_yr <- avgRSyr(ppt8100, "pr", yrs8100)
ppt0121_avg_yr <- avgRSyr(ppt0121, "pr", yrs0121)
pet8100_avg_yr <- avgRSyr(pet8100, "pet", yrs8100)
pet0121_avg_yr <- avgRSyr(pet0121, "pet", yrs0121)


# resample to GSOC extent
tmp8100_avg_yr <- resample(readAll(tmp8100_avg_yr), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
tmp0121_avg_yr <- resample(readAll(tmp0121_avg_yr), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")

ppt8100_avg_yr <- resample(readAll(ppt8100_avg_yr), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
ppt0121_avg_yr <- resample(readAll(ppt0121_avg_yr), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")

pet8100_avg_yr <- resample(readAll(pet8100_avg_yr), gsoc2, method = "bilinear", progress = "text")
pet0121_avg_yr <- resample(readAll(pet0121_avg_yr), gsoc2, method = "bilinear", progress = "text")


# cache averages by year
# writeRaster(tmp8100_avg_yr, filename = "CONUS_TEMP_Stack_81-00_TC_yr.tif", overwrite = TRUE)
# writeRaster(tmp0121_avg_yr, filename = "CONUS_TEMP_Stack_01-19_TC_yr.tif", overwrite = TRUE)
# writeRaster(ppt8100_avg_yr, filename = "CONUS_Prec_Stack_81-00_TC_yr.tif", overwrite = TRUE)
# writeRaster(ppt0121_avg_yr, filename = "CONUS_Prec_Stack_01-19_TC_yr.tif", overwrite = TRUE)
# writeRaster(pet8100_avg_yr, filename = "CONUS_PET_Stack_81-00_TC_yr.tif", overwrite = TRUE)
# writeRaster(pet0121_avg_yr, filename = "CONUS_PET_Stack_01-19_TC_yr2.tif", overwrite = TRUE)


# load cached averages by year
tmp8100_avg_yr <- stack("CONUS_TEMP_Stack_81-00_TC_yr.tif")
tmp0121_avg_yr <- stack("CONUS_TEMP_Stack_01-19_TC_yr.tif")
ppt8100_avg_yr <- stack("CONUS_Prec_Stack_81-00_TC_yr.tif")
ppt0121_avg_yr <- stack("CONUS_Prec_Stack_01-19_TC_yr.tif")
pet8100_avg_yr <- stack("CONUS_PET_Stack_81-00_TC_yr.tif")
pet0121_avg_yr <- stack("CONUS_PET_Stack_01-19_TC_yr.tif")


# Calculate eq 1 from MIAMI MODEL (g DM/m2/day)
npp8100_ppt <- 3000 * (1 - exp(-0.000664 * ppt8100_avg_yr))

# Calculate eq 2 from MIAMI MODEL (g DM/m2/day)
npp8100_tmp <- 3000 / (1 + exp(1.315 - 0.119 * tmp8100_avg_yr))

npp8100 <- lapply(1:20, function(i) {
  cat("calculating min from ", i, "\n")
  npp <- min(npp8100_tmp[[i]], npp8100_ppt[[i]], na.rm = TRUE)
})
npp8100 <- stack(npp8100)
npp8100_TnCHaYr <- npp8100 * 1/100 * 0.5
writeRaster(npp8100_TnCHaYr, filename = "CONUS_NPP_MIAMI_tnC_Ha_Year_STACK_81-00.tif", overwrite = TRUE, options = c("COMPRESS=DEFLATE"))

npp8100_TnCHaYr_avg <- mean(npp8100_TnCHaYr) 
writeRaster(npp8100_TnCHaYr_avg, filename = "CONUS_NPP_MIAMI_MEAN_81-00_AOI.tif", overwrite = TRUE, options = c("COMPRESS=DEFLATE"))


# Uncertainties ----

library(terra)

tmp8100 <- rast("CONUS_AverageTemperature_1981-2001.tif")
ppt8100 <- rast("CONUS_Precipitation_1981-2001.tif")


# Min
tmp8100_min <- stack(tmp8100 * 1.02)
ppt8100_min <- stack(ppt8100 * 0.95)

tmp8100_min_avg <- avgRSyr(tmp8100_min, "tmmx", yrs8100) * 0.1
ppt8100_min_avg <- avgRSyr(ppt8100_min, "pr", yrs8100)

npp8100_min_avg_tmp <- 3000 * (1 - exp(-0.000664 * ppt8100_min_avg))
npp8100_min_avg_ppt <- 3000 / (1 + exp(1.315 - 0.119 * tmp8100_min_avg))

npp8100_min <- lapply(1:20, function(i) {
  cat("calculating min from ", i, "\n")
  npp <- min(npp8100_min_avg_tmp[[i]], npp8100_min_avg_tmp[[i]], na.rm = TRUE)
})
npp8100_min <- stack(npp8100_min)
npp8100_min_TnCHaYr <- stack(npp8100_min) * 1/100 * 0.5

npp8100_min_TnCHaYr <- resample(npp8100_min_TnCHaYr, gsoc2, method = "bilinear", progress = "text")

writeRaster(npp8100_min_TnCHaYr, 
            filename = "CONUS_NPP_MIAMI_tnC_Ha_Year_STACK_81-00_MIN.tif", 
            overwrite = TRUE, 
            options = c("COMPRESS=DEFLATE")
)

npp8100_min_TnCHaYr_avg <- mean(npp8100_min_TnCHaYr) 
writeRaster(npp8100_min_TnCHaYr_avg, 
            filename = "CONUS_NPP_MIAMI_MEAN_81-00_AOI_MIN.tif", 
            overwrite = TRUE, 
            options = c("COMPRESS=DEFLATE")
)


# Max
tmp8100_max <- stack(tmp8100 * 0.98)
ppt8100_max <- stack(ppt8100 * 1.05)

tmp8100_max_avg <- avgRSyr(tmp8100_max, "tmmx", yrs8100) * 0.1
ppt8100_max_avg <- avgRSyr(ppt8100_max, "pr", yrs8100)

npp8100_max_avg_tmp <- 3000 * (1 - exp(-0.000664 * ppt8100_max_avg))
npp8100_max_avg_ppt <- 3000 / (1 + exp(1.315 - 0.119 * tmp8100_max_avg))

npp8100_max <- lapply(1:20, function(i) {
  cat("calculating min from ", i, "\n")
  npp <- min(npp8100_max_avg_tmp[[i]], npp8100_max_avg_tmp[[i]], na.rm = TRUE)
})
npp8100_max <- stack(npp8100_max)
npp8100_max_TnCHaYr <- stack(npp8100_max) * 1/100 * 0.5

npp8100_max_TnCHaYr <- resample(readAll(npp8100_max_TnCHaYr), gsoc2, method = "bilinear", progress = "text")

writeRaster(npp8100_max_TnCHaYr, 
            filename = "CONUS_NPP_MIAMI_tnC_Ha_Year_STACK_81-00_MAX.tif", 
            overwrite = TRUE, 
            options = c("COMPRESS=DEFLATE")
)

npp8100_max_TnCHaYr_avg <- mean(npp8100_max_TnCHaYr) 
writeRaster(npp8100_max_TnCHaYr_avg, 
            filename = "CONUS_NPP_MIAMI_MEAN_81-00_AOI_MAX.tif", 
            overwrite = TRUE, 
            options = c("COMPRESS=DEFLATE")
)



# Vegetation Cover ----
vars <- paste0("NDVI_2000-2020_prop_gt_06_CR_MES_", formatC(1:12, width = 2, flag = "0"), "_conus.tif")
veg <- stack(vars)
veg <- (veg * -0.4) + 1
veg2 <- resample(readAll(veg), gsoc2, method = "bilinear", progress = "text")
# writeRaster(veg2, filename = "CONUS_Cov_stack_AOI.tif", options = c("COMPRESS=DEFLATE"), progress = "text", overwrite = TRUE)



# Landcover ----
# http://www.fao.org/geonetwork/srv/en/main.home?uuid=ba4526fd-cdbf-4028-a1bd-5a559c4bff38
lc <- raster("D:/geodata/land_use_land_cover/GlcShare_v10_Dominant/glc_shv10_DOM.Tif")
projection(lc) <- "+init=epsg:4326"

lc2 <- crop(lc, gsoc2, progress = "text")
lc2 <- resample(readAll(lc2), gsoc2, method = "ngb", progress = "text", datatype = "INT1S")
# writeRaster(lc2, file = "CONUS_glc_shv10_DOM.tif", overwrite = TRUE)
LU_AOI <- readAll(raster("CONUS_glc_shv10_DOM.tif"))



# DR layer ----
dr <- (LU_AOI == 2 | LU_AOI == 12| LU_AOI == 13) * 1.44 + (LU_AOI == 4) *  0.25 + (LU_AOI == 3 | LU_AOI == 5 | LU_AOI == 6 | LU_AOI == 8) * 0.67
# writeRaster(dr, file = "CONUS_glc_shv10_DOM_DR.tif", overwrite = TRUE)



# Stack layers ----
# Spin up layers
vars <- c(SOC  = "CONUS_GSOCmap1.5.0.tif",
          CLAY = "CONUS_gnatsgo_fy20_1km_clay_wt.tif",
          TEMP = "CONUS_Temp_Stack_81-00_TC.tif",
          PREC = "CONUS_Prec_Stack_81-00_TC.tif",
          PET  = "CONUS_PET_Stack_81-00_TC.tif",
          LU   = "CONUS_glc_shv10_DOM.tif",
          DR   = "CONUS_glc_shv10_DOM_DR.tif",
          COV  = "CONUS_Cov_stack_AOI.tif"
)
su_rs <- readAll(stack(vars))
# writeRaster(su_rs, filename = "Stack_Set_SPIN_UP_AOI.tif", format = "GTiff", progress = "text", overwrite = TRUE)

su_rs <- readAll(su_rs)
su_pts <- rasterToPoints(su_rs, spatial = TRUE, progress = "text")
su_pts2 <- subset(su_pts, as.character(LU) %in% c(2, 3, 5, 12))
saveRDS(su_pts, file = "su_sdf.RDS")


# Warm up layers
n_wu <- 18

vars <- c(SOC  = "CONUS_GSOCmap1.5.0.tif",
          CLAY = "CONUS_gnatsgo_fy20_1km_clay_wt.tif",
          TEMP = "CONUS_Temp_Stack_228_01_19_TC.tif",
          PREC = "CONUS_Prec_Stack_228_01_19_TC.tif",
          PET  = "CONUS_PET_Stack_228_01-19_TC.tif",
          LU   = "CONUS_glc_shv10_DOM.tif",
          DR   = "CONUS_glc_shv10_DOM_DR.tif",
          COV  = "CONUS_Cov_stack_AOI.tif",
          NPP  = "CONUS_NPP_MIAMI_MEAN_81-00_AOI.tif",
          NPP_MIN = "CONUS_NPP_MIAMI_MEAN_81-00_AOI_MIN.tif",
          NPP_MAX = "CONUS_NPP_MIAMI_MEAN_81-00_AOI_MAX.tif"
)
# LU <- stack(replicate(n_wu, raster("CONUS_glc_shv10_DOM.tif")))
wu_rs <- rast(vars)

tiles <- st_as_sf(st_make_grid(aoi, n = c(6, 6)))
tiles$tile <- 1:nrow(tiles)


# iterate over tiles and extract points
test <- lapply(1:34, function(x) {
  
  cat("cropping ", x, as.character(Sys.time()), "\n")
  
  # crop raster stack
  tiles_sub <- vect(tiles[tiles$tile == x, ])
  # tiles_sub <- as(tiles[tiles$tile == x, ], "Spatial")
  rs_sub    <- crop(wu_rs, tiles_sub)
  
  # mask non-ag landuses
  LU <- c(2, 3, 5, 12)
  rs_sub$CONUS_glc_shv10_DOM[! rs_sub$CONUS_glc_shv10_DOM %in% LU] <- NA
  
  # extract points
  lu_pts <- as.points(rs_sub$CONUS_glc_shv10_DOM)
  wu_pts <- extract(rs_sub, lu_pts, xy = TRUE)
  
  # save points
  saveRDS(wu_pts, file = paste0("wu_pts_tile_", x, ".rds"))
  
  return(dim(wu_pts))
})


f <- paste0("wu_pts_tile_", 1:34, ".rds")
wu_pts <- lapply(f, function(x){
  temp <- readRDS(file = x)
  temp <- temp[complete.cases(temp), ]
  # nm     <- names(temp)
  # idx    <- which(!grepl("_tmmx$|_pr$|_pet$", nm))
  # id_tmp <- which(grepl("_tmmx$", nm))[1]
  # id_ppt <- which(grepl("_pr$", nm))[1]
  # id_pet <- which(grepl("_pet$", nm))[1]
  # 
  # temp <- temp[c(
  #   idx, 
  #   id_tmp:(id_tmp + 11),
  #   id_ppt:(id_ppt + 11),
  #   id_pet:(id_pet + 11)
  # )]
  # dim(temp)
  })
wu_pts <- do.call("rbind", wu_pts)
# saveRDS(wu_pts, file = "wu_pts.rds")
# file.remove(f)


# Forward Stack
vars <- c(SOC  = "CONUS_GSOCmap1.5.0.tif",
          CLAY = "CONUS_gnatsgo_fy20_1km_clay_wt.tif",
          TEMP = "CONUS_Temp_Stack_01-19_TC_yr.tif",
          PREC = "CONUS_Prec_Stack_228_01_19_TC.tif",
          PET  = "CONUS_PET_Stack_228_01-19_TC.tif",
          DR   = "CONUS_glc_shv10_DOM_DR.tif",
          LU   = "CONUS_glc_shv10_DOM.tif",
          COV  = "CONUS_Cov_stack_AOI.tif"
)
fs_rs <- stack(vars)
# writeRaster(fs_rs, filename = "Stack_Set_FOWARD.tif", progress = "text", overwrite = TRUE)




# Define Target Points ----
pts <- raster("CONUS_glc_shv10_DOM.tif")
pts <- rasterToPoints(pts, spatial = TRUE)

pts$CONUS_glc_shv10_DOM <- as.character(pts$CONUS_glc_shv10_DOM)
pts <- subset(pts, CONUS_glc_shv10_DOM %in% c(2, 3, 5, 12))
table(pts$CONUS_glc_shv10_DOM)

pts <- as(pts, 'SpatialPoints')
write_sf(st_as_sfc(pts), dsn = "target_points.shp")

