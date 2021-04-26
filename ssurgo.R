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

h_us      <- horizons(f_us_seg)
h_statsgo <- horizons(f_statsgo_seg)




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



# rasterize SSURGO ----

library(raster)
# library(terra)

r <- raster("D:/geodata/soils/gnatsgo_fy20_90m.tif")
# r2 <- r[1:100, 1:100, drop = FALSE]
mu_agg2 <- subset(mu_agg, source == "gnatsgo")[2:6]
names(mu_agg2)[1] <- "ID"
levels(r) <- mu_agg2
r <- readAll(r)

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


system('"C:/OSGeo4W64/bin/gdalwarp.exe" -overwrite  -te -2356155 270045 2263815 3172635 -tr 1000 1000 -t_srs "EPSG:5070" -ot "Int32" -r "average" -of "GTiff" "D:/geodata/project_data/gsp-gsocseq/ssurgo_fy20_90m_clay_wt.tif" "D:/geodata/project_data/gsp-gsocseq/ssurgo_fy20_1km_clay_wt.tif"')


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


# rasterize STATSGO2

library(sf)
library(raster)
library(gdalUtils)

r  <- raster("D:/geodata/soils/gnatsgo_fy20_30m_compress.tif")
writeRaster(r, filename = "D:/geodata/soils/gstatsgo2_fy20_30m.tif", datatype = "INT4U", progress = "text", overwrite = TRUE)

mu_agg$mukey2 <- as.character(mu_agg$mukey)

gsm <- read_sf(dsn = "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us.shp")
gsm$mukey2 <- as.integer(gsm$MUKEY)
gsm <- left_join(gsm, mu_agg[mu_agg$source == "statsgo", ], by = c("MUKEY" = "mukey2"))
gsm <- st_transform(gsm, crs = proj4string(r))
# write_sf(gsm, dsn = "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us2.shp", delete_dsn = TRUE)
gsm <- read_sf(dsn = "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us2.shp")



# mukey
vars <- c("mukey", "soc_wt", "clay_wt")[3]
lapply(vars, function(x) {
  cat("rasterizing ", x, "\n")
  
  prec <- ifelse(x == "soc_wt", "Float32", "Int32")
  
  system(paste0('"C:/OSGeo4W64/bin/gdal_rasterize.exe" -a "', x, '" -l "gsmsoilmu_a_us2" -te -2356155 269635 2263845 3172635 -tr 1000 1000 -ot ', prec, ' -a_nodata -9999 "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us2.shp" "D:/geodata/project_data/gsp-gsocseq/gstatsgo2_fy20_1km_"', x, '".tif"'))
})

system('"C:/OSGeo4W64/bin/gdal_rasterize.exe" -a "mukey2" -l "gsmsoilmu_a_us2" -te -2356155 270015 2263815 3172635 -tr 30 30 -ot Int32 "D:/geodata/soils/wss_gsmsoil_US_[2016-10-13]/spatial/gsmsoilmu_a_us2.shp" "D:/geodata/soils/gstatsgo2_fy20_30m.tif"')



system('"C:/OSGeo4W64/bin/gdalwarp.exe" -overwrite  -te -2356155 270015 2263815 3172635 -tr 1000 1000 -t_srs "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs" -ot "Int32" -r "near" -of "GTiff" "D:/geodata/soils/gstatsgo2_fy20_30m.tif" "D:/geodata/soils/gstatsgo2_fy20_1km.tif"')



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



