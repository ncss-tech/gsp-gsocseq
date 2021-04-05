library(aqp)
library(soilDB)
library(dplyr)

setwd("D:/geodata/project_data/gsp-gsocseq")


# fetch SSURGO ----

f_us <- fetchGDB(dsn ="D:/geodata/soils/gNATSGO_CONUS_FY20.gdb", WHERE = "areasymbol LIKE '%'")
mu_us <- get_mapunit_from_GDB(dsn = "D:/geodata/soils/gNATSGO_CONUS_FY20.gdb", stats = TRUE)
# save(f_us, mu_us, file = "D:/geodata/project_data/gsp-bs/data/gnatsgo.RData")
load(file = "D:/geodata/project_data/gsp-bs/data/gnatsgo.RData")


# fetch STATSGO2, horizon data seems to be missing from gNATSGO
mu_statsgo <- get_mapunit_from_SDA(WHERE = "areasymbol = 'US'")
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

h_us      <- horizons(f_us_seg)
h_statsgo <- horizons(f_statsgo_seg)

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
  mutate(soc = om_r / 1.724,
         thk = hzdepb_r - hzdept_r
         ) %>%
  group_by(source, cokey) %>%
  summarize(soc_wt  = weighted.mean(soc,         w = thk, na.rm = TRUE),
            clay_wt = weighted.mean(claytotal_r, w = thk, na.rm = TRUE) 
            ) %>%
  right_join(site(f_us), by = "cokey") %>%
  right_join(mu_us,      by = "mukey") %>%
  group_by(source, mukey, musym, muname) %>%
  summarize(
    soc_wt  = round(weighted.mean(soc_wt,  w = comppct_r, na.rm = TRUE), 2),
    clay_wt = round(weighted.mean(clay_wt, w = comppct_r, na.rm = TRUE))
    ) %>%
  ungroup() %>%
  as.data.frame()
mu_agg <- mu_agg[c(2:6, 1)]
# write.csv(mu_agg, file = "gnatsgo_gsoc.csv", row.names = FALSE)
mu_agg <- read.csv(file = "gnatsgo_gsoc.csv", stringsAsFactors = FALSE)



# rasterize ----

library(raster)

r <- raster("D:/geodata/soils/gnatsgo_fy20_1km.tif")
# r2 <- r[1:100, 1:100, drop = FALSE]
# r2 <- ratify(r2)
# rat <- levels(r2)[[1]]
names(mu_agg)[1] <- "ID"
levels(r) <- subset(mu_agg, source == "gnatsgo")

vars <- c("soc_wt", "clay_wt")
lapply(vars, function(x) {
  cat(x, as.character(Sys.time()), "\n")
  # beginCluster(type = "SOCK")
  deratify(r, att = x, 
           filename = paste0("D:/geodata/project_data/gsp-gsocseq/gnatsgo_fy20_1km_", x, ".tif"),
           options = c("COMPRESS=DEFLATE"), 
           overwrite = TRUE, 
           progress = "text" 
  )
  # endCluster()
})


gnatsgo <- "D:/geodata/soils/gnatsgo_fy20_30m.tif"
r <- raster(gnatsgo)

gdalUtils::gdalwarp(
  srcfile = gnatsgo,
  dstfile = gsub("30m.tif", "120m.tif", gnatsgo),
  t_srs = proj4string(r), 
  te = bbox(r),
  r = "mode",
  tr = c(120, 120),
  ot = "Int32",
  verbose   = TRUE,
  overwrite = TRUE
)


# rasterize

library(sf)
library(raster)
library(gdalUtils)

r  <- raster("D:/geodata/soils/gnatsgo_fy20_30m_compress.tif")
writeRaster(r, filename = "D:/geodata/soils/gstatsgo2_fy20_30m.tif", datatype = "INT4U", progress = "text", overwrite = TRUE)

gsm <- read_sf(dsn = "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us.shp")
gsm$mukey2 <- as.integer(gsm$MUKEY)
gsm <- st_transform(gsm, crs = proj4string(r))
write_sf(gsm, dsn = "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us2.shp", delete_dsn = TRUE)

system('"C:/OSGeo4W64/bin/gdal_rasterize.exe" -a "mukey2" -l "gsmsoilmu_a_us2" -te -2356155 270015 2263815 3172635 -tr 30 30 -ot Int32 "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us2.shp" "D:/geodata/soils/gstatsgo2_fy20_30m.tif"')

system('"C:/OSGeo4W64/bin/gdalwarp.exe" -overwrite  -te -2356155 270015 2263815 3172635 -tr 1000 1000 -t_srs "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs" -ot "Int32" -r "near" -of "GTiff" "D:/geodata/soils/gstatsgo2_fy20_30m.tif" "D:/geodata/soils/gstatsgo2_fy20_1km.tif"')

