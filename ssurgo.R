library(aqp)
library(soilDB)
library(dplyr)
library(sf)
library(rmapshaper)

setwd("D:/geodata/project_data/gsp-gsocseq")


# fetch gNATSGO ----

# f_us <- fetchGDB(dsn ="D:/geodata/soils/gNATSGO_CONUS_FY20.gdb", WHERE = "areasymbol LIKE '%'")
mu_us <- get_mapunit_from_GDB(dsn = "D:/geodata/soils/gNATSGO_CONUS_FY20.gdb", stats = TRUE)
mu_us$source <- "gnatsgo"
# save(f_us, file = "D:/geodata/project_data/gsp-bs/data/gnatsgo.RData")
load(file = "D:/geodata/project_data/gsp-bs/data/gnatsgo.RData")


f_ak  <- fetchGDB(dsn ="D:/geodata/soils/gNATSGO_AK_july2020.gdb", WHERE = "areasymbol LIKE '%'")
mu_ak <- get_mapunit_from_GDB(dsn = "D:/geodata/soils/gNATSGO_AK_july2020.gdb", stats = TRUE)
mu_ak$source <- "AK"
# save(f_ak, mu_ak, file = "D:/geodata/project_data/gsp-bs/data/gnatsgo_ak.RData")
load(file = "D:/geodata/project_data/gsp-bs/data/gnatsgo_ak.RData")



# fetchSSURGO ----
mu_us$idx <- rep(1:100, length.out = nrow(mu_us))
f_ssurgo <- {
  split(mu_us, mu_us$idx) ->.;
  lapply(., function(x) {
    cat("getting ", unique(x$idx), "\n")
    temp <- fetchSDA(WHERE = paste("areasymbol != 'US' AND ", "mukey IN ('", paste0(x$mukey, collapse = "', '"), "')"), childs = FALSE, duplicates = TRUE)
  }) ->.;
}
# save(f_ssurgo, file = "f_ssurgo_list.RData")
load(file = "f_ssurgo_list.RData")

f_s <- do.call("rbind", lapply(f_ssurgo, site))
f_h <- do.call("rbind", lapply(f_ssurgo, horizons))

f_s <- f_s[!duplicated(paste(f_s$nationalmusym, f_s$cokey)), ]
f_h <- f_h[!duplicated(paste(f_h$cokey, f_h$chkey)), ]

load(file = "D:/geodata/project_data/gsp-bs/data/gnatsgo.RData")


# fetch STATSGO2 ----
# some horizon data seems to be missing from gNATSGO
mu_statsgo <- get_mapunit_from_SDA(WHERE = "areasymbol = 'US'")
mu_statsgo$source <- "statsgo"
mu_statsgo$idx <- rep(1:4, length.out = nrow(mu_statsgo))
f_statsgo <- {
  split(mu_statsgo, mu_statsgo$idx) ->.;
  lapply(., function(x) {
    cat("getting ", unique(x$idx), "\n")
    temp <- fetchSDA(WHERE = paste("areasymbol = 'US' AND ", "mukey IN ('", paste0(x$mukey, collapse = "', '"), "')"), childs = FALSE, duplicates = TRUE)
  }) ->.;
}
f_statsgo <- aqp::combine(f_statsgo)
# saveRDS(f_statsgo, file = "statsgo_sda.rds")
f_statsgo <- readRDS(file = "statsgo_sda.rds")



# compute weighted averages ----

f_us_seg      <- segment(f_us, intervals = c(0, 30))
f_statsgo_seg <- segment(f_statsgo, intervals = c(0, 30))
f_ak_seg      <- segment(f_ak, intervals = c(0, 30))

h_us      <- horizons(f_us_seg)
h_statsgo <- horizons(f_statsgo_seg)
h_ak      <- horizons(f_ak_seg)



# combine ssurgo and statsgo
idx <- names(mu_us)[names(mu_us) %in% names(mu_statsgo)]
mu <- rbind(mu_us[idx], mu_statsgo[idx])

idx <- siteNames(f_us)[siteNames(f_us) %in% siteNames(f_statsgo)]
s <- rbind(cbind(source = "gnatsgo", site(f_us)[idx]), 
           cbind(source = "statsgo", site(f_statsgo)[idx])
)

idx <- sapply(h_us, is.factor)
h_us[idx] <- lapply(h_us[idx], as.character)
idx <- sapply(h_statsgo, is.factor)
h_statsgo[idx] <- lapply(h_statsgo[idx], as.character)

nm <- names(h_us)[names(h_us) %in% names(h_statsgo)]

h <- rbind(
  cbind(source = "gnatsgo", h_us[nm]),
  cbind(source = "statsgo", h_statsgo[nm])
)


mu_agg <- h %>%
  mutate(soc   = om_r / 1.724,
         thk   = hzdepb_r - hzdept_r,
         stock = thk * (soc * dbthirdbar_r / 100) * ((100 - fragvol_r) / 100)
  ) %>%
  group_by(source, cokey) %>%
  summarize(soc_wt  = weighted.mean(soc,         w = thk, na.rm = TRUE),
            clay_wt = weighted.mean(claytotal_r, w = thk, na.rm = TRUE) 
  ) %>%
  right_join(s,  by = c("cokey", "source")) %>%
  right_join(mu, by = c("mukey", "source")) %>%
  group_by(source, mukey, musym, muname) %>%
  summarize(
    soc_wt  = round(weighted.mean(soc_wt,  w = comppct_r, na.rm = TRUE), 2),
    clay_wt = round(weighted.mean(clay_wt, w = comppct_r, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  as.data.frame()
# write.csv(mu_agg, file = "gnatsgo_gsoc.csv", row.names = FALSE)
mu_agg <- read.csv(file = "gnatsgo_gsoc.csv", stringsAsFactors = FALSE)


s_ak <- cbind(source = "AK", site(f_ak))
mu_ak_agg <- h_ak %>%
  mutate(source = "AK", 
         soc   = om_r / 1.724,
         thk   = hzdepb_r - hzdept_r
  ) %>%
  group_by(source, cokey) %>%
  summarize(clay_wt = weighted.mean(claytotal_r, w = thk, na.rm = TRUE)) %>%
  right_join(s_ak,  by = c("cokey", "source")) %>%
  right_join(mu_ak, by = c("mukey", "source")) %>%
  group_by(source, mukey, musym, muname) %>%
  summarize(
    clay_wt = round(weighted.mean(clay_wt, w = comppct_r, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  as.data.frame()
# write.csv(mu_ak_agg, file = "gnatsgo_ak_gsoc.csv", row.names = FALSE)
mu_ak_agg <- read.csv(file = "gnatsgo_ak_gsoc.csv", stringsAsFactors = FALSE)



# rasterize ----

## STATSGO2 ----

# create temp
r  <- raster("D:/geodata/soils/gnatsgo_fy20_30m_compress.tif")
writeRaster(r, filename = "D:/geodata/soils/gstatsgo2_fy20_30m.tif", datatype = "INT4U", progress = "text", overwrite = TRUE)


# load polygons and tidy
gsm <- read_sf(dsn = "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us.shp")
gsm$mukey2 <- as.integer(gsm$MUKEY)
gsm <- st_transform(gsm, crs = proj4string(r))
write_sf(gsm, dsn = "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us2.shp", delete_dsn = TRUE)


# rasterize
gdalUtilities::gdal_rasterize(
  src_datasource = "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us2.shp",
  a              = "mukey2",
  dst_filename   = "D:/geodata/soils/gstatsgo2_fy20_30m.tif",
  of             = "GTiff",
  te             = bbox(r),
  tr             = res(r),
  co             = c("COMPRESS=DEFLATE"),
  a_nodata       = -99999
)



# deratify ----

setwd("D:/geodata/project_data/gsp-gsocseq")

library(raster)
library(terra)

mu_agg <- read.csv(file = "gnatsgo_gsoc.csv", stringsAsFactors = FALSE)
mu_agg$mukey2 <- as.character(mu_agg$mukey)


## gNATSGO ----

# upscale
gnatsgo   <- "D:/geodata/soils/gnatsgo_fy20_30m.tif"
gnatsgo90 <- gsub("30m", "90m", gnatsgo)
# if (file.exists(gnatsgo90)) file.remove(gnatsgo90)
r <- raster(gnatsgo)
gdalUtilities::gdalwarp(
  srcfile   = gnatsgo,
  dstfile   = gnatsgo90,
  tr        = c(90, 90),
  ot        = "Int32",
  r         = "mode",
  overwrite = TRUE
)


# attach attributes
r <- raster("D:/geodata/soils/gnatsgo_fy20_90m.tif")
# r2 <- r[1:100, 1:100, drop = FALSE]
mu_agg2 <- subset(mu_agg, source == "gnatsgo")[2:6]
names(mu_agg2)[1] <- "ID"
levels(r) <- mu_agg2
r <- readAll(r)


# deratify
vars <- c("soc_wt", "clay_wt")
lapply(vars, function(x) {
  cat(x, as.character(Sys.time()), "\n")
  # beginCluster(type = "SOCK")
  deratify(r, att = x, 
           filename = paste0("D:/geodata/project_data/gsp-gsocseq/ssurgo_fy20_90m_", x, ".tif"),
           options = c("COMPRESS=DEFLATE"), 
           overwrite = TRUE, 
           progress = "text" 
  )
  # endCluster()
})


# upscale again
gdalUtilities::gdalwarp(
  srcfile = "D:/geodata/project_data/gsp-gsocseq/gnatsgo_ak_fy20_90m_clay_wt.tif", 
  dstfile = "D:/geodata/project_data/gsp-gsocseq/gnatsgo_ak_fy20_1km_clay_wt.tif",
  tr = c(1000, 1000),
  ot = "Int32",
  r = "average",
  srcnodata = -2147483648,
  dstnodata = -9999,
  overwrite = TRUE
)



## AK ----

# upscale
ak   <- "D:/geodata/soils/gnatsgo_ak_july2020_30m.tif"
ak90 <- gsub("30m", "90m", gnatsgo)
# if (file.exists(ak90)) file.remove(akgo90)
r  <- raster(ak)
gdalUtilities::gdalwarp(
  srcfile   = ak, 
  dstfile   = ak90,
  tr        = c(90, 90),
  ot        = "Int32",
  r         = "mode",
  overwrite = TRUE
)

# attach attributes
r2 <- raster("D:/geodata/soils/gnatsgo_ak_july2020_90m.tif")
mu_agg2 <- subset(mu_ak_agg, source == "AK")[2:5]
names(mu_agg2)[1] <- "ID"
levels(r2) <- mu_agg2
r2 <- readAll(r2)


# deratify
vars <- c("clay_wt")
lapply(vars, function(x) {
  cat(x, as.character(Sys.time()), "\n")
  # beginCluster(type = "SOCK")
  deratify(r2, att = x, 
           filename = paste0("D:/geodata/project_data/gsp-gsocseq/gnatsgo_ak_fy20_90m_", x, ".tif"),
           options = c("COMPRESS=DEFLATE"), 
           overwrite = TRUE, 
           progress = "text" 
  )
  # endCluster()
})


# upscale again
gdalUtilities::gdalwarp(
  srcfile = "D:/geodata/project_data/gsp-gsocseq/gnatsgo_ak_fy20_90m_clay_wt.tif", 
  dstfile = "D:/geodata/project_data/gsp-gsocseq/gnatsgo_ak_fy20_1km_clay_wt.tif",
  tr = c(1000, 1000),
  ot = "Int32",
  r = "average",
  srcnodata = -2147483648,
  dstnodata = -9999,
  overwrite = TRUE
)



# combine SSURGO & STATSGO2

lf <- list.files(getwd(), pattern = ".tif$")

ssurgo_st  <- readAll(raster("D:/geodata/project_data/gsp-gsocseq/ssurgo_fy20_1km_clay_wt.tif"))
statsgo_st <- readAll(raster("D:/geodata/project_data/gsp-gsocseq/gstatsgo2_fy20_1km_clay_wt.tif"))

gnatsgo_clay <- merge(
  ssurgo_st, 
  statsgo_st,
  filename = "gnatsgo_fy20_1km_clay_wt.tif",
  overwrite = TRUE
)


gnatsgo_soc <- merge(
  ssurgo_st$gnatsgo_fy20_1km_soc_wt, 
  statsgo_st$gstatsgo2_fy20_30m_soc_wt,
  filename = "gnatsgo_fy20_1km_soc_wt.tif"
)



