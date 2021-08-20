library(sf)
library(raster)
library(terra)
library(rmapshaper)

setwd("D:/geodata/project_data/gsp-gsocseq")

source('C:/Users/stephen.roecker/OneDrive - USDA/projects/gsp-gsocseq/code/functions.R')


# AOI ----
## CONUS ----
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


## OCONUS ----
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

# # manually construct AK bufbox
# sapol_bufbox$AK1 <- sapol_bufbox$AK
# st_bbox(sapol_bufbox$AK1)
# st_bbox(sapol_bufbox$AK2)
# sapol_bufbox$AK <- st_as_sf(
#   st_as_sfc(
#     st_bbox(c(
#       xmin = -187.6,
#       xmax = -129.87002,
#       ymin = 51.07513,
#       ymax = 71.68593
#     ))),
#   crs = st_crs(4326)  
#   )
# sapol_bufbox$AK$asym = "AK"

test <- sapol_bufbox$AK2
mapview::mapview(test)

lapply(sapol_bufbox, function(x) write_sf(x, dsn = paste0("AOIs/AOI_OCONUS_", x$asym, "_bufbox.shp")))
lapply(sapol_bndy, function(x)   write_sf(x, dsn = paste0("AOIs/OCONUS_", x$asym, "_boundry.shp")))



# GSOC ----
setwd("D:/geodata/project_data/gsp-gsocseq")
gsoc  <- rast("D:/geodata/soils/GSOCmap1.5.0.tif")

## CONUS ----
aoi   <- read_sf(dsn = "AOIs/AOI_CONUS_bufbox.shp")
gsoc2 <- crop(gsoc, aoi)
# writeRaster(gsoc2, filename =  "CONUS/CONUS_GSOCmap1.5.0.tif", progress = "text", overwrite = TRUE)
gsoc2 <- raster("CONUS/CONUS_GSOCmap1.5.0.tif")


## AK ----
ex <- ext(c(-187.693995463933, -129.8666668672, 51.0744147105421, 71.6827479614154))
aoi_ak1 <- read_sf(dsn = "AOIs/AOI_OCONUS_AK_bufbox.shp")
aoi_ak2 <- read_sf(dsn = "AOIs/AOI_OCONUS_AK2_bufbox.shp")

gsoc_ak1 <- crop(gsoc, vect(aoi_ak1), filename = "AK/AK1_GSOCmap1.5.0.tif")
gsoc_ak2 <- crop(gsoc, vect(aoi_ak2), filename = "AK/AK2_GSOCmap1.5.0.tif")

# # construct extent manually
# ak_ext <- ext(c(-187.693995463933, -129.8666668672, 51.0744147105421, 71.6827479614154))
# ak_r   <- rast(ak_ext, res = res(gsoc), crs = crs(gsoc))
# 
# gsoc_ak1 <- resample(gsoc_ak1, ak_r)
# gsoc_ak2 <- resample(gsoc_ak2, ak_r)
# 
# gsoc_ak <- merge(gsoc_ak1, gsoc_ak2, filename =  "AK/AK_GSOCmap1.5.0.tif", overwrite = TRUE)
# gsoc2 <- rast("AK_GSOCmap1.5.0.tif")


## CMS ----
cms_soc <- raster("D:/geodata/soils/CMS_SOC_Mexico_CONUS_1737/data/SOC_prediction_1991_2010.tif")
cms_soc_1km <- aggregate(cms_soc, fact = 4, fun = mean, progress = "text")
cms_soc_1km <- readAll(cms_soc_1km)
cms_soc_1km <- projectRaster(from = cms_soc_1km, to = ssurgo_st, filename = "cms_soc_1km.tif", method = "bilinear", progress = "text", overwrite = TRUE)



# clay ----
## CONUS ----
gdalUtilities::gdalwarp(
  srcfile = "gnatsgo_fy20_1km_clay_wt.tif",
  dstfile = "CONUS_gnatsgo_fy20_1km_clay_wt.tif",
  t_srs   = crs(gsoc2),
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear"
)
# clay <- raster("gnatsgo_fy20_1km_clay_wt.tif")
# clay_4326 <- projectRaster(
#   clay, 
#   gsoc2, 
#   filename = "CONUS_gnatsgo_fy20_1km_clay_wt.tif", 
#   method = "bilinear", 
#   progress = "text", 
#   overwrite = TRUE, 
#   datatype = "INT4S"
# )


## AK ----
aoi <- "AK1"
aoi <- "AK2"

gsoc_ak <- rast(paste0(aoi, "_GSOCmap1.5.0.tif"))
clay_ak <- rast("gnatsgo_ak_fy20_1km_clay_wt.tif")
clay_ak <- resample(
  clay_ak, 
  gsoc_ak, 
  filename = paste0(aoi, "_gnatsgo_fy20_1km_clay_wt.tif"),
  method   = "bilinear",
  overwrite = TRUE, 
  datatype = "INT4S"
)



# TerraClimate ----
setwd("D:/geodata/project_data/gsp-gsocseq/CONUS")

aoi <- "CONUS"

## Resample original data to GSOC extent ----
gsoc2 <- rast(paste0(aoi, "_GSOCmap1.5.0.tif"))


# tmp8100
gdalUtilities::gdalwarp(
  srcfile = paste0(aoi, "_AverageTemperature_1981-2001.tif"),
  dstfile = paste0(aoi, "_Temp_Stack_228_81_00_TC.tif"),
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  srcnodata = -3.4e+38,
  dstnodata = -99999
)
# tmp0120
gdalUtilities::gdalwarp(
  srcfile = paste0(aoi, "_AverageTemperature_2001-2021.tif"),
  dstfile = paste0(aoi, "_Temp_Stack_228_01_19_TC.tif"),
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  srcnodata = -3.4e+38,
  dstnodata = -99999
)
# ppt8100
gdalUtilities::gdalwarp(
  srcfile = paste0(aoi, "_Precipitation_1981-2001.tif"),
  dstfile = paste0(aoi, "_Prec_Stack_228_81_00_TC.tif"),
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  srcnodata = -3.4e+38,
  dstnodata = -99999
)
# ppt0120
gdalUtilities::gdalwarp(
  srcfile = paste0(aoi, "_Precipitation_2001-2021.tif"),
  dstfile = paste0(aoi, "_Prec_Stack_228_01_19_TC.tif"),
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  srcnodata = -3.4e+38,
  dstnodata = -99999
)
# pet8100
gdalUtilities::gdalwarp(
  srcfile = paste0(aoi, "_PET_1981-2001.tif"),
  dstfile = paste0(aoi, "_PET_Stack_228_81-00_TC.tif"),
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  srcnodata = -3.4e+38,
  dstnodata = -99999
)
# pet0120
gdalUtilities::gdalwarp(
  srcfile = paste0(aoi, "_PET_2001-2021.tif"),
  dstfile = paste0(aoi, "_PET_Stack_228_01-19_TC.tif"),
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear",
  srcnodata = -3.4e+38,
  dstnodata = -99999
)


# ## merge rasters ----
# aoi  <- "AK"
# 
# aois <- c("AK1", "AK2")
# vars <- list(
#   tmp8100 = "Temp_Stack_228_81_00_TC",
#   tmp0121 = "Temp_Stack_228_01_20_TC",
#   ppt8100 = "Prec_Stack_228_81_00_TC",
#   ppt0121 = "Prec_Stack_228_01_20_TC",
#   pet8100 = "PET_Stack_228_81-00_TC",
#   pet0121 = "PET_Stack_228_01-20_TC"
# )
# 
# fn   <- lapply(vars, function(x) paste0(paste(aois, x, sep = "_"), ".tif"))
# 
# lapply(fn, function(x) {
#   cat("merging", x, as.character(Sys.time()), "\n")
#   
#   # r1 <- readAll(raster(x[1]))
#   # r2 <- readAll(raster(x[2]))
#   r1 <- rast(x[1])
#   r2 <- rast(x[2])
#   
#   test2 = merge(
#     r1, r2,
#     filename = gsub("AK1", "AK", x[1]),
#     overwrite = TRUE
#   )
#   #  gdalUtilities::gdalwarp(
#   #   srcfile = x,
#   #   dstfile = gsub("AK1", "AK", x[1])
#   # )
# })



## load data ----
aoi <- "CONUS"
tmp8100 <- rast(paste0(aoi, "_Temp_Stack_228_81_00_TC.tif"))
tmp0119 <- rast(paste0(aoi, "_Temp_Stack_228_01_19_TC.tif"))
ppt8100 <- rast(paste0(aoi, "_Prec_Stack_228_81_00_TC.tif"))
ppt0119 <- rast(paste0(aoi, "_Prec_Stack_228_01_19_TC.tif"))
pet8100 <- rast(paste0(aoi, "_PET_Stack_228_81-00_TC.tif"))
pet0119 <- rast(paste0(aoi, "_PET_Stack_228_01-19_TC.tif"))



# create list of raster names by month
mos <- formatC(1:12, width = 2, flag = "0")
names(mos) <- month.abb
mos8100 <- lapply(mos, function(x) paste0(1981:2000, x))
mos0119 <- lapply(mos, function(x) paste0(2001:2019, x))



## average months ----
tmp8100_avg <- avgRS(tmp8100, "tmmx", mos8100) * 0.1
tmp0119_avg <- avgRS(tmp0119, "tmmx", mos0119) * 0.1

ppt8100_avg <- avgRS(ppt8100, "pr", mos8100)
ppt0119_avg <- avgRS(ppt0119, "pr", mos0119)

pet8100_avg <- avgRS(pet8100, "pet", mos8100) * 0.1
pet0119_avg <- avgRS(pet0119, "pet", mos0119) * 0.1



# cache files
writeRaster(
  tmp8100_avg, datatype = "INT2S",  gdal = c("COMPRESS=DEFLATE"),
  filename = paste0(aoi, "_Temp_Stack_81-00_TC.tif"), overwrite = TRUE
  )
writeRaster(
  tmp0119_avg, datatype = "INT2S", gdal = c("COMPRESS=DEFLATE"),
  filename = paste0(aoi, "_Temp_Stack_01-19_TC.tif"), overwrite = TRUE
  )
writeRaster(
  ppt8100_avg, datatype = "INT2S", gdal = c("COMPRESS=DEFLATE"), 
  filename = paste0(aoi, "_Prec_Stack_81-00_TC.tif"), overwrite = TRUE
  )
writeRaster(
  ppt0119_avg, datatype = "INT2S", gdal = c("COMPRESS=DEFLATE"), 
  filename = paste0(aoi, "_Prec_Stack_01-19_TC.tif"), overwrite = TRUE
  )
writeRaster(
  pet8100_avg, gdal = c("COMPRESS=DEFLATE"), 
  filename = paste0(aoi, "_PET_Stack_81-00_TC.tif"), overwrite = TRUE
  )
writeRaster(
  pet0119_avg, gdal = c("COMPRESS=DEFLATE"),
  filename = paste0(aoi, "_PET_Stack_01-19_TC.tif"), overwrite = TRUE
  )


# load cached files
tmp8100_avg <- rast(paste0(aoi, "_Temp_Stack_81-00_TC.tif"))
tmp0121_avg <- rast(paste0(aoi, "_Temp_Stack_01-19_TC.tif"))
ppt8100_avg <- rast(paste0(aoi, "_Prec_Stack_81-00_TC.tif"))
ppt0121_avg <- rast(paste0(aoi, "_Prec_Stack_01-19_TC.tif"))
pet8100_avg <- rast(paste0(aoi, "_PET_Stack_81-00_TC.tif"))
pet0120_avg <- rast(paste0(aoi, "_PET_Stack_01-19_TC.tif"))



## average by year ----
# create list of raster names by year
yrs8100 <- lapply(1981:2000, function(x) paste0(x))
yrs8100 <- lapply(yrs8100, function(x) paste0(x, formatC(1:12, width = 2, flag = "0")))
yrs0119 <- lapply(2001:2019, function(x) paste0(x))
yrs0119 <- lapply(yrs0119, function(x) paste0(x, formatC(1:12, width = 2, flag = "0")))


tmp8100_avg_yr <- avgRS(tmp8100, "tmmx", yrs8100) * 0.1
tmp0119_avg_yr <- avgRS(tmp0119, "tmmx", yrs0119) * 0.1
ppt8100_avg_yr <- avgRS(ppt8100, "pr",   yrs8100)
ppt0119_avg_yr <- avgRS(ppt0119, "pr",   yrs0119)
pet8100_avg_yr <- avgRS(pet8100, "pet",  yrs8100) * 0.1
pet0119_avg_yr <- avgRS(pet0119, "pet",  yrs0119) * 0.1


# # resample to GSOC extent
# tmp8100_avg_yr <- resample(readAll(tmp8100_avg_yr), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
# tmp0121_avg_yr <- resample(readAll(tmp0121_avg_yr), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
# 
# ppt8100_avg_yr <- resample(readAll(ppt8100_avg_yr), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
# ppt0121_avg_yr <- resample(readAll(ppt0121_avg_yr), gsoc2, method = "bilinear", progress = "text", datatype = "INT2S")
# 
# pet8100_avg_yr <- resample(readAll(pet8100_avg_yr), gsoc2, method = "bilinear", progress = "text")
# pet0120_avg_yr <- resample(readAll(pet0120_avg_yr), gsoc2, method = "bilinear", progress = "text")


# cache averages by year
writeRaster(
  tmp8100_avg_yr, datatype = "INT2S",  gdal = c("COMPRESS=DEFLATE"),
  filename = paste0(aoi, "_Temp_Stack_81-00_TC_yr.tif"), overwrite = TRUE
)
writeRaster(
  tmp0119_avg_yr, datatype = "INT2S", gdal = c("COMPRESS=DEFLATE"),
  filename = paste0(aoi, "_Temp_Stack_01-19_TC_yr.tif"), overwrite = TRUE
)
writeRaster(
  ppt8100_avg_yr, datatype = "INT2S", gdal = c("COMPRESS=DEFLATE"), 
  filename = paste0(aoi, "_Prec_Stack_81-00_TC_yr.tif"), overwrite = TRUE
)
writeRaster(
  ppt0119_avg_yr, datatype = "INT2S", gdal = c("COMPRESS=DEFLATE"), 
  filename = paste0(aoi, "_Prec_Stack_01-19_TC_yr.tif"), overwrite = TRUE
)
writeRaster(
  pet8100_avg_yr, gdal = c("COMPRESS=DEFLATE"), 
  filename = paste0(aoi, "_PET_Stack_81-00_TC_yr.tif"), overwrite = TRUE
)
writeRaster(
  pet0119_avg_yr, gdal = c("COMPRESS=DEFLATE"),
  filename = paste0(aoi, "_PET_Stack_01-19_TC_yr.tif"), overwrite = TRUE
)



# load cached averages by year
tmp8100_avg_yr <- rast(paste0(aoi, "_Temp_Stack_81-00_TC_yr.tif"))
tmp0121_avg_yr <- rast(paste0(aoi, "_Temp_Stack_01-19_TC_yr.tif"))
ppt8100_avg_yr <- rast(paste0(aoi, "_Prec_Stack_81-00_TC_yr.tif"))
ppt0121_avg_yr <- rast(paste0(aoi, "_Prec_Stack_01-19_TC_yr.tif"))
pet8100_avg_yr <- rast(paste0(aoi, "_PET_Stack_81-00_TC_yr.tif"))
pet0120_avg_yr <- rast(paste0(aoi, "_PET_Stack_01-19_TC_yr.tif"))


## calculate NPP ----
# Calculate eq 1 from MIAMI MODEL (g DM/m2/day)
npp8100_ppt <- 3000 * (1 - exp(-0.000664 * ppt8100_avg_yr))

# Calculate eq 2 from MIAMI MODEL (g DM/m2/day)
npp8100_tmp <- 3000 / (1 + exp(1.315 - 0.119 * tmp8100_avg_yr))

npp8100 <- lapply(1:20, function(i) {
  cat("calculating min from ", i, "\n")
  npp <- min(npp8100_tmp[[i]], npp8100_ppt[[i]], na.rm = TRUE)
})
npp8100 <- rast(npp8100)
npp8100_TnCHaYr <- npp8100 * 1/100 * 0.5
writeRaster(npp8100_TnCHaYr, filename = paste0(aoi, "_NPP_MIAMI_tnC_Ha_Year_STACK_81-00.tif"), overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))

npp8100_TnCHaYr_avg <- mean(npp8100_TnCHaYr) 
names(npp8100_TnCHaYr_avg) <- "NPP_MIAMI_MEAN_81-00_AOI"
writeRaster(npp8100_TnCHaYr_avg, filename = paste0(aoi, "_NPP_MIAMI_MEAN_81-00_AOI.tif"), overwrite = TRUE, gdal = c("COMPRESS=DEFLATE"))



# Uncertainties ----

## Min ----
# tmp8100 <- rast(paste0(aoi, "_Temp_Stack_228_81_00_TC.tif"))
# ppt8100 <- rast(paste0(aoi, "_Prec_Stack_228_81_00_TC.tif"))
# 
# tmp8100_min <- tmp8100 * 1.02
# ppt8100_min <- ppt8100 * 0.95
# 
# tmp8100_min_avg <- avgRS(tmp8100_min, "tmmx", yrs8100) * 0.1
# ppt8100_min_avg <- avgRS(ppt8100_min, "pr", yrs8100)

tmp8100_avg_yr <- rast(paste0(aoi, "_Temp_Stack_81-00_TC_yr.tif"))
ppt8100_avg_yr <- rast(paste0(aoi, "_Prec_Stack_81-00_TC_yr.tif"))

tmp8100_min_avg <- tmp8100_avg_yr * 1.02
ppt8100_min_avg <- ppt8100_avg_yr * 0.95

npp8100_min_avg_ppt <- 3000 * (1 - exp(-0.000664     * ppt8100_min_avg))
npp8100_min_avg_tmp <- 3000 / (1 + exp(1.315 - 0.119 * tmp8100_min_avg))

npp8100_min <- lapply(1:20, function(i) {
  cat("calculating min from ", i, "\n")
  npp <- min(npp8100_min_avg_ppt[[i]], npp8100_min_avg_tmp[[i]], na.rm = TRUE)
})
npp8100_min <- rast(npp8100_min)
npp8100_min_TnCHaYr <- npp8100_min * 1/100 * 0.5

writeRaster(npp8100_min_TnCHaYr, 
            filename = paste0(aoi, "_NPP_MIAMI_tnC_Ha_Year_STACK_81-00_MIN.tif"), 
            overwrite = TRUE, 
            gdal = c("COMPRESS=DEFLATE")
)

npp8100_min_TnCHaYr_avg <- mean(npp8100_min_TnCHaYr) 
names(npp8100_min_TnCHaYr_avg) <- "NPP_MIAMI_MEAN_81-00_AOI_MIN"
writeRaster(npp8100_min_TnCHaYr_avg, 
            filename = paste0(aoi, "_NPP_MIAMI_MEAN_81-00_AOI_MIN.tif"), 
            overwrite = TRUE, 
            gdal = c("COMPRESS=DEFLATE")
)


## Max ----
# tmp0119 <- rast(paste0(aoi, "_Temp_Stack_228_01_19_TC.tif"))
# ppt0119 <- rast(paste0(aoi, "_Prec_Stack_228_01_19_TC.tif"))
# 
# tmp8100_max <- tmp8100 * 0.98
# ppt8100_max <- ppt8100 * 1.05
# 
# tmp8100_max_avg <- avgRS(tmp8100_max, "tmmx", yrs8100) * 0.1
# ppt8100_max_avg <- avgRS(ppt8100_max, "pr",   yrs8100)

tmp8100_avg_yr <- rast(paste0(aoi, "_Temp_Stack_81-00_TC_yr.tif"))
ppt8100_avg_yr <- rast(paste0(aoi, "_Prec_Stack_81-00_TC_yr.tif"))

tmp8100_max_avg <- tmp8100_avg_yr * 0.98
ppt8100_max_avg <- ppt8100_avg_yr * 1.05

npp8100_max_avg_ppt <- 3000 * (1 - exp(-0.000664     * ppt8100_max_avg))
npp8100_max_avg_tmp <- 3000 / (1 + exp(1.315 - 0.119 * tmp8100_max_avg))

npp8100_max <- lapply(1:20, function(i) {
  cat("calculating min from ", i, "\n")
  npp <- min(npp8100_max_avg_ppt[[i]], npp8100_max_avg_tmp[[i]], na.rm = TRUE)
})
npp8100_max         <- rast(npp8100_max)
npp8100_max_TnCHaYr <- npp8100_max * 1/100 * 0.5

# npp8100_max_TnCHaYr <- resample(readAll(npp8100_max_TnCHaYr), gsoc2, method = "bilinear", progress = "text")

writeRaster(npp8100_max_TnCHaYr, 
            filename = paste0(aoi, "_NPP_MIAMI_tnC_Ha_Year_STACK_81-00_MAX.tif"), 
            overwrite = TRUE, 
            gdal = c("COMPRESS=DEFLATE")
)

npp8100_max_TnCHaYr_avg <- mean(npp8100_max_TnCHaYr) 
names(npp8100_max_TnCHaYr_avg) <- "NPP_MIAMI_MEAN_81-00_AOI_MAX"
writeRaster(npp8100_max_TnCHaYr_avg, 
            filename = paste0(aoi, "_NPP_MIAMI_MEAN_81-00_AOI_MAX.tif"), 
            overwrite = TRUE, 
            gdal = c("COMPRESS=DEFLATE")
)



# Vegetation Cover ----
aoi <- "CONUS"
gsoc2 <- rast(paste0(aoi, "_GSOCmap1.5.0.tif"))

# vars <- paste0(aoi, "_NDVI_2000-2020_prop_gt_06_CR_MES_", formatC(1:12, width = 2, flag = "0"), ".tif")
vars <- paste0("NDVI_2000-2020_prop_gt_06_CR_MES_", formatC(1:12, width = 2, flag = "0"), "_conus.tif")
veg <- rast(vars)
veg[is.na(veg)] <- 0
veg <- (veg * -0.4) + 1
names(veg) <- paste0("COV_", formatC(1:12, width = 2, flag = "0"))

writeRaster(veg, 
            filename = paste0(aoi, "_Cov_stack_AOI.tif"), 
            overwrite = TRUE, 
            gdal = c("COMPRESS=DEFLATE")
)

fn <- paste0(aoi, "_Cov_stack_AOI2.tif")
if (file.exists(fn)) file.remove(fn)
gdalUtilities::gdalwarp(
  srcfile = paste0(aoi, "_Cov_stack_AOI.tif"),
  dstfile = paste0(aoi, "_Cov_stack_AOI2.tif"),
  te      = c(bbox(gsoc2)),
  tr      = res(gsoc2),
  r       = "bilinear"
)

veg <- rast(paste0(aoi, "_Cov_stack_AOI2.tif"))
veg2 <- mask(veg, gsoc2)
writeRaster(veg2, filename = paste0(aoi, "_Cov_stack_AOI.tif"), gdal = c("COMPRESS=DEFLATE"), overwrite = TRUE)



# Landcover ----
# http://www.fao.org/geonetwork/srv/en/main.home?uuid=ba4526fd-cdbf-4028-a1bd-5a559c4bff38
aoi <- "CONUS"
gsoc2 <- rast(paste0(aoi, "_GSOCmap1.5.0.tif"))

lc <- rast("D:/geodata/land_use_land_cover/GlcShare_v10_Dominant/glc_shv10_DOM.Tif")
crs(lc) <- "+init=epsg:4326"

lc2 <- crop(lc, gsoc2) #)
lc2 <- resample(lc2, gsoc2, method = "ngb", datatype = "INT1S")
writeRaster(lc2, filename = paste0(aoi, "_glc_shv10_DOM.tif"), overwrite = TRUE)
LU_AOI <- rast(paste0(aoi, "_glc_shv10_DOM.tif"))



# DR layer ----
dr <- (LU_AOI == 2 | LU_AOI == 12| LU_AOI == 13) * 1.44 + (LU_AOI == 4) *  0.25 + (LU_AOI == 3 | LU_AOI == 5 | LU_AOI == 6 | LU_AOI == 8) * 0.67
writeRaster(dr, file = paste0(aoi, "_glc_shv10_DOM_DR.tif"), overwrite = TRUE)



# Stack layers ----

aoi <- "CONUS"
setwd("D:/geodata/project_data/gsp-gsocseq/CONUS")

## Spin up layers ----
vars <- c(SOC  = paste0(aoi, "_GSOCmap1.5.0.tif"),
          CLAY = paste0(aoi, "_gnatsgo_fy20_1km_clay_wt.tif"),
          TEMP = paste0(aoi, "_Temp_Stack_81-00_TC.tif"),
          PREC = paste0(aoi, "_Prec_Stack_81-00_TC.tif"),
          PET  = paste0(aoi, "_PET_Stack_81-00_TC.tif"),
          LU   = paste0(aoi, "_glc_shv10_DOM.tif"),
          DR   = paste0(aoi, "_glc_shv10_DOM_DR.tif"),
          COV  = paste0(aoi, "_Cov_stack_AOI.tif")
)
su_rs <- lapply(names(vars), function(x) {
  r <- rast(vars[x])
  n <- nlyr(r)
  names(r) <- if (n > 1) {
    paste0(x, "_", formatC(1:nlyr(r), width = 2, flag = "0"))
  } else x
  return(r)
  })
su_rs <- rast(su_rs)
# writeRaster(su_rs, filename = "Stack_Set_SPIN_UP_AOI.tif", format = "GTiff", progress = "text", overwrite = TRUE)

su_df <- as.data.frame(su_rs, xy = TRUE, cells = TRUE, na.rm = FALSE)
su_df <- su_df[
  as.character(su_df$LU) %in% c(2, 3, 5, 12, 13) 
  & complete.cases(su_df), 
]
su_sf <- st_as_sf(
  su_df,
  coords = c("x", "y"),
  crs    = 4326
)
saveRDS(su_sf, file = paste0(aoi, "_su_sf.RDS"))



## Warm up layers ----
su_sf <- readRDS(file = paste0(aoi, "_su_sf.RDS"))

n_wu <- 18

vars <- c(SOC  = paste0(aoi, "_GSOCmap1.5.0.tif"),
          CLAY = paste0(aoi, "_gnatsgo_fy20_1km_clay_wt.tif"),
          TEMP = paste0(aoi, "_Temp_Stack_228_01_19_TC.tif"),
          PREC = paste0(aoi, "_Prec_Stack_228_01_19_TC.tif"),
          PET  = paste0(aoi, "_PET_Stack_228_01-19_TC.tif"),
          LU   = paste0(aoi, "_glc_shv10_DOM.tif"),
          DR   = paste0(aoi, "_glc_shv10_DOM_DR.tif"),
          COV  = paste0(aoi, "_Cov_stack_AOI.tif"),
          NPP  = paste0(aoi, "_NPP_MIAMI_MEAN_81-00_AOI.tif"),
          NPP_MIN = paste0(aoi, "_NPP_MIAMI_MEAN_81-00_AOI_MIN.tif"),
          NPP_MAX = paste0(aoi, "_NPP_MIAMI_MEAN_81-00_AOI_MAX.tif")
)
# LU <- stack(replicate(n_wu, raster("CONUS_glc_shv10_DOM.tif")))
wu_rs  <- rast(vars)
# wu_rs  <- stack(vars)


### partition points (if too large) ----
su_sf$idx <- as.integer(
  cut(su_sf$id, 
      breaks = quantile(su_sf$id, p = seq(0, 1, 0.1)), 
      include.lowest = TRUE)
)

test <- lapply(1:10, function(x) {
  cat("extracting part", x, as.character(Sys.time()), "\n")
  temp  <- vect(su_pts[which(su_pts$idx == x), ])
  # temp  <- su_pts[1:100, ]
  wu_ex <- extract(wu_rs, temp, xy = TRUE)
  # wu_ex <- as.data.frame(extract(wu_rs, temp, sp = TRUE))
  wu_ex <- cbind(idx = x, wu_ex)
  saveRDS(wu_ex, file = paste0("wu_pts_sub_", x, "_v2.rds"))
})

f_p <- paste0("wu_pts_sub_", 1:10, "_v2.rds")

wu_pts_p1 <- lapply(f_p[1], function(x){
  temp <- readRDS(file = x)
})
wu_pts_p1 <- data.table::rbindlist(wu_pts_p1)
data.table::fwrite(wu_pts_p1, file = "wu_pts_p1.csv")

wu_pts_p2 <- lapply(f_p[6:10], function(x) {
  temp <- readRDS(file = x)
})
wu_pts_p2 <- data.table::rbindlist(wu_pts_p2)
data.table::fwrite(wu_pts_p2, file = "wu_pts_p2.csv")

wu_pts_p1 <- data.table::fread(file = "wu_pts_p1.csv")
wu_pts_p2 <- data.table::fread(file = "wu_pts_p2.csv")
wu_pts <- rbind(wu_pts_p1, wu_pts_p2)
data.table::fwrite(wu_pts, file = "wu_pts_v2.csv")


### don't partition points (if # of points is small) ----
su_v  <- vect(su_sf)
wu_ex <- extract(wu_rs, su_v, xy = TRUE, cells = TRUE)
saveRDS(wu_ex, file = paste0(aoi, "_wu_df.rds"))



# tiles <- st_as_sf(st_make_grid(aoi, n = c(6, 6)))
# tiles$tile <- 1:nrow(tiles)
# 
# # iterate over tiles and extract points
# test <- lapply(1:34, function(x) {
#   
#   cat("cropping ", x, as.character(Sys.time()), "\n")
#   
#   # crop raster stack
#   tiles_sub <- vect(tiles[tiles$tile == x, ])
#   # tiles_sub <- as(tiles[tiles$tile == x, ], "Spatial")
#   rs_sub    <- crop(wu_rs, tiles_sub)
#   
#   # mask non-ag landuses
#   LU <- c(2, 3, 5, 12)
#   rs_sub$CONUS_glc_shv10_DOM[! rs_sub$CONUS_glc_shv10_DOM %in% LU] <- NA
#   
#   # extract points
#   lu_pts <- as.points(rs_sub$CONUS_glc_shv10_DOM)
#   wu_pts <- extract(rs_sub, lu_pts, xy = TRUE)
#   
#   # save points
#   saveRDS(wu_pts, file = paste0("wu_pts_tile_", x, "_v2.rds"))
#   
#   return(dim(wu_pts))
# })
# 
# 
# f_p1 <- paste0("wu_pts_tile_", 1:21, "_v2.rds")
# f_p2 <- paste0("wu_pts_tile_", 22:34, "_v2.rds")
# 
# wu_pts_p1 <- lapply(f_p1, function(x){
#   temp <- readRDS(file = x)
#   temp <- temp[complete.cases(temp), ]
#   })
# # saveRDS(do.call("rbind", wu_pts_p1), file = "wu_pts_v2_part1.rds")
# wu_pts_p2 <- lapply(f_p2, function(x){
#   temp <- readRDS(file = x)
#   temp <- temp[complete.cases(temp), ]
# })
# # saveRDS(do.call("rbind", wu_pts_p2), file = "wu_pts_v2_part2.rds")
# # file.remove(f)


## Forward Run ----
aoi <- "AK2"

su_sf <- readRDS(file = paste0(aoi, "_su_sf.RDS"))

vars <- c(SOC  = paste0(aoi, "_GSOCmap1.5.0.tif"),
          CLAY = paste0(aoi, "_gnatsgo_fy20_1km_clay_wt.tif"),
          TEMP = paste0(aoi, "_Temp_Stack_01-19_TC.tif"),
          PREC = paste0(aoi, "_Prec_Stack_01-19_TC.tif"),
          PET  = paste0(aoi, "_PET_Stack_01-19_TC.tif"),
          LU   = paste0(aoi, "_glc_shv10_DOM.tif"),
          DR   = paste0(aoi, "_glc_shv10_DOM_DR.tif"),
          COV  = paste0(aoi, "_Cov_stack_AOI.tif")
)
fr_rs <- lapply(names(vars), function(x) {
  r <- rast(vars[x])
  n <- nlyr(r)
  names(r) <- if (n > 1) {
    paste0(x, "_", formatC(1:nlyr(r), width = 2, flag = "0"))
  } else x
  return(r)
})
fr_rs <- rast(fr_rs)

su_v  <- vect(su_sf)
fr_df <- extract(fr_rs, su_v, xy = TRUE, cells = TRUE)
saveRDS(fr_df, file = paste0(aoi, "_fr_df.RDS"))



## AK ----
su_ak1 <- readRDS("AK1_su_sf.RDS")
su_ak2 <- readRDS("AK2_su_sf.RDS")
wu_ak1 <- readRDS("AK1_wu_df.rds")
wu_ak2 <- readRDS("AK2_wu_df.rds")
fr_ak1 <- readRDS("AK1_fr_df.rds")
fr_ak2 <- readRDS("AK2_fr_df.rds")

su_ak <- rbind(
    cbind(aoi = "AK1", su_ak1),
    cbind(aoi = "AK2", su_ak2)
  )
su_ak <- cbind(sf::st_coordinates(su_ak), su_ak)
saveRDS(su_ak, file = "su_sf.RDS")

wu_ak <- rbind(
    cbind(aoi = "AK1", wu_ak1),
    cbind(aoi = "AK2", wu_ak2)
)
saveRDS(wu_ak, file = "wu_df.RDS")
data.table::fwrite(wu_ak, file = "wu_df.csv")

fr_ak <- rbind(
  cbind(aoi = "AK1", fr_ak1),
  cbind(aoi = "AK2", fr_ak2)
)
saveRDS(fr_ak, file = "fr_df.RDS")


# inspect outputs
write_sf(test2, dsn = "test.gpkg", driver = "GPKG", overwrite = TRUE)
gdalUtilities::gdal_rasterize(
  src_datasource = "test.gpkg",
  a              = "CLAY",
  dst_filename   = "test_clay.tif",
  of             = "GTiff",
  te             = bbox(gsoc),
  tr             = res(gsoc),
  co             = c("COMPRESS=DEFLATE"),
  a_srs          = "EPSG:4326",
  a_nodata       = -999
)



# Define Target Points ----
pts <- raster("CONUS_glc_shv10_DOM.tif")
pts <- rasterToPoints(pts, spatial = TRUE)

pts$CONUS_glc_shv10_DOM <- as.character(pts$CONUS_glc_shv10_DOM)
pts <- subset(pts, CONUS_glc_shv10_DOM %in% c(2, 3, 5, 12))
table(pts$CONUS_glc_shv10_DOM)

pts <- as(pts, 'SpatialPoints')
write_sf(st_as_sfc(pts), dsn = "target_points.shp")

