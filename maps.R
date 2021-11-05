
# Setup ----
library(sf)
library(terra)

aoi <- "CONUS"
setwd(paste0("D:/geodata/project_data/gsp-gsocseq/", aoi))


fr_df  <- readRDS(file = "_gsocseq_maps_epac.rds")
fr_df  <- readRDS(file = "rothC_fr_final_analytical.rds")
# fr_df  <- cbind(aoi = aoi, fr_df)


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
  # # Uncertainties for the Relative difference  SSM_ - SOCBAU
  RSR_SSM1_Uncertainty = sqrt((unc_ssm * f_t.low)^2  + (unc_bau * f_t.bau)^2) /
    abs(f_t.bau + f_t.low)
  RSR_SSM2_Uncertainty = sqrt((unc_ssm * f_t.med)^2  + (unc_bau * f_t.bau)^2) /
    abs(f_t.bau + f_t.med)
  RSR_SSM3_Uncertainty = sqrt((unc_ssm * f_t.high)^2 + (unc_bau * f_t.bau)^2) /
    abs(f_t.bau + f_t.high)
})

vars <- c("finalSOC_BAU_", "T0_", "f_t.low", "f_t.med", "f_t.high", "f_t.medmin", "f_t.medmax")
vars2 <- grepl("^RSR_|Uncertainty", names(gsocseq_maps))
idx <- (
  rowSums(gsocseq_maps[vars]  > 800) + 
    rowSums(gsocseq_maps[vars]  < 0)   +
    rowSums(gsocseq_maps[vars2] < 0)
) > 0
sum(idx)
# idx2 <- apply(gsocseq_maps[vars], 1, function(x) any(x > 800) | any(x < 0))
gsocseq_maps <- gsocseq_maps[!idx, ]

names(gsocseq_maps) <- gsub("\\.", "_", names(gsocseq_maps))
nm    <- names(gsocseq_maps)
vars  <- c("aoi", "x", "y", "cell", "LU", "CLAY", "SOC", "CinputFORWARD_r")
vars2 <- nm[grep("^finalSOC|^T0|_Uncertainty$|^AbsDiff|^RelDiff|^ASR|^RSR", nm)]
gsocseq_maps <- gsocseq_maps[c(vars, vars2)]
saveRDS(gsocseq_maps, file = "gsocseq_maps.rds")


# convert to points ----
fr_sf <- st_as_sf(
  gsocseq_maps,
  coords = c("x", "y"),
  crs    = 4326
)

ak1_fr_sf <- fr_sf[fr_sf$aoi == "AK1", ]
ak2_fr_sf <- fr_sf[fr_sf$aoi == "AK2", ]

write_sf(ak1_fr_sf, dsn = paste0("AK1_gsocseq_maps.gpkg"), driver = "GPKG", overwrite = TRUE) 
write_sf(ak2_fr_sf, dsn = paste0("AK2_gsocseq_maps.gpkg"), driver = "GPKG", overwrite = TRUE) 

write_sf(fr_sf, dsn = paste0("CONUS_gsocseq_maps_as_final.gpkg"), driver = "GPKG", overwrite = TRUE) 



# rasterize points ----
aoi <- "AK1"
gsoc <- rast(paste0(aoi, "_GSOCmap1.5.0.tif"))
gsoc[!is.na(gsoc)] <- 1
lu   <- rast(paste0(aoi, "_glc_shv10_DOM.tif"))
lu   <- lu %in% c(2, 3, 5, 12, 13)
gsoc <- gsoc * lu


lapply(vars2, function(x) {
  
  f  <- paste0(aoi, "_GSOCseq_", x, "Map030.tif")
  
  cat("rasterizing ", f, as.character(Sys.time()), "\n")
  
  if (file.exists(f)) file.remove(f)
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


## tabulate ----
summary(gsocseq_maps)

gm <- gsocseq_maps

idx <- grepl("^ASR_SSM._$|^ASR_BAU_$", names(gm))
asr_avg <- aggregate(x = gm[idx], by = list(LU = gm$LU), mean)
asr_avg2 <- aggregate(x = gm[idx], by = list(LU = rep(0, length(gm$LU))), mean)
asr_avg  <- rbind(asr_avg, asr_avg2)

idx <- grepl("^RSR_SSM._$|^RSR_BAU_$", names(gm))
rsr_avg <- aggregate(x = gm[idx], by = list(LU = gm$LU), mean)
rsr_avg2 <- aggregate(x = gm[idx], by = list(LU = rep(0, length(gm$LU))), mean)
rsr_avg  <- rbind(rsr_avg, rsr_avg2)

avg <- cbind(asr_avg[c(1, 5:2)], rsr_avg[4:2])
names(avg) <- gsub("_$", "", names(avg))
signif(avg, 2)


test <- data.frame(t(sapply(vars2, function(x) rbind(quantile(gsocseq_maps[, x], na.rm = TRUE)))))
test <- round(test, 2)
names(test) <- paste0("p", seq(0, 1, 0.25))
test <- cbind(var = row.names(test), test)
row.names(test) <- NULL
View(test)
write.csv(test, "gsocseq_summaries_final.csv", row.names = FALSE)



### LRR ----
mlra_sf <- read_sf("D:/geodata/soils/MLRA_52.shp")
mlra_sf$MLRA_ID <- as.character(mlra_sf$MLRA_ID)

gdalUtilities::gdal_rasterize(
  src_datasource = "D:/geodata/soils/MLRA_52.shp",
  dst_filename   = "D:/geodata/soils/MLRA_52.tif",
  a              = "MLRA_ID",
  of             = "GTiff",
  te             = ext(gsoc2)[c(1, 3, 2, 4)],
  tr             = res(gsoc2),
  a_nodata       = -99999
  )
mlra_r  <- rast("D:/geodata/soils/MLRA_52.tif")
mlra_df <- terra::extract(mlra_r, gsocseq_maps[, c("x", "y")])

gsocseq_maps <- cbind(gsocseq_maps, mlra_df)
gsocseq_maps$MLRA_52 <- as.character(gsocseq_maps$MLRA_52)
vars <- c("MLRA_ID", "LRRSYM")
gsocseq_maps <- merge(gsocseq_maps, mlra_sf[vars], by.x = "MLRA_52", by.y = "MLRA_ID", all.x = TRUE)



## plot ----
aoi <- "AK1"
f   <- paste0(aoi, "_GSOCseq_", vars2, "Map030.tif")
rs  <- rast(f) 
names(rs) <- vars2


vars <- c(pred     = "final|T0_$", 
          pred_unc = "^..._Uncertainty$|^.._Uncertainty",
          abs_diff = "AbsDiff", 
          rel_diff = "RelDiff",
          asr      = "ASR_SSM._$|ASR_BAU_$",
          rsr      = "RSR_SSM._$",
          asr_unc  = "ASR_...._Uncertainty",
          rsr_unc  = "RSR_...._Uncertainty"
)

aoi <- 'AK1'
lapply(seq_along(vars), function(i) {
  idx <- grepl(vars[i], names(rs))
  brks <- quantile(values(rs[[idx]]), probs = seq(0, 1, 0.1), na.rm = TRUE)
  png(paste0("plots_", aoi, "_", names(vars)[i], ".png"), units = "in", width = 12, height = 6, res = 300)
  plot(rs[[idx]], breaks = brks, col = viridis::viridis(10))
  dev.off()
})



# OCONUS ----

setwd("D:/geodata/project_data/gsp-gsocseq/OCONUS")

lf <- list.files()
lf <- lf[grepl(".tif$", lf)]

vars <- substr(lf, 1, 7)
lf <- by(lf,  vars, rast)


## filter extreme values ----
lf2 <- lapply(unique(vars), function(x) {
  # subset list
  gu  <- lf[[x]]
  idx <- grep("Uncertainty", names(gu))
  
  # tally extreme values
  idx2 <- sum(gu[[idx]] < 0) | sum(gu[[idx]] > 200)
  t1 <- sum(values(idx2) > 0, na.rm = TRUE)
  t2 <- sum(values(idx2) < 1, na.rm = TRUE)
  print(c(t1, t2))
  
  # # filter extreme values
  # idx2 <- ifel(idx2 < 1, 1, NA)
  # gu2 <- gu * idx2
  # 
  # lapply(1:nlyr(gu2), function(i) {
  #   writeRaster(gu2[[i]], filename = paste0("./v2/", names(gu2)[i], ".tif"), overwrite = TRUE, NAflag = -9999)
  # })
})



## QA statistics ----
setwd("D:/geodata/project_data/gsp-gsocseq/OCONUS/v2")

lf <- list.files()
lf <- lf[grepl(".tif$", lf)]
aois <- substr(lf, 1, 7)
lf <- by(lf,  aois, rast)

test <- lapply(lf, function(x) {
  temp <- apply(values(x), 2, function(x2) quantile(x2, na.rm = TRUE))
  round(t(temp), 2)
})

test <- as.data.frame(do.call("rbind", test))
var <- row.names(test)
test <- cbind(aoi = substr(var, 1, 7), var = var, test)
row.names(test) <- NULL
View(test)
write.csv(test, "gsocseq_summaries_final.csv", row.names = FALSE)



## plot ----
vars <- c(pred     = "final|T0_$", 
          pred_unc = "BAU_Uncertainty|T0_Uncertainty|SSM_Uncertainty",
          abs_diff = "AbsDiff", 
          rel_diff = "RelDiff",
          asr      = "ASR_SSM._Map030|ASR_BAU_Map030",
          rsr      = "RSR_SSM._Map030",
          asr_unc  = "ASR_...._Uncertainty",
          rsr_unc  = "RSR_...._Uncertainty"
)

lapply(unique(aois), function(x) {
  
  rs <- lf[[x]]
  
  lapply(seq_along(vars), function(i) {
    idx <- grep(vars[i], names(rs))
    brks <- quantile(values(rs[[idx]]), probs = seq(0, 1, 0.1), na.rm = TRUE)
    png(paste0("plots_", x, "_", names(vars)[i], ".png"), units = "in", width = 12, height = 6, res = 300)
    plot(rs[[idx]], breaks = brks, col = viridis::viridis(10))
    dev.off()
  })
})



# Create maps ----

library(tmap)


setwd("D:/geodata/project_data/gsp-gsocseq/CONUS")


# read state boundaries
aoi <- read_sf(dsn = "../AOIs/CONUS.shp") %>%
  rmapshaper::ms_simplify(drop_null_geometries = FALSE) %>%
  st_make_valid()
# mutate(CONUS = TRUE)
# rmapshaper::ms_dissolve(field = "CONUS")
plot(aoi)


# construct dataframe of filenames, variables, senarios, and units
lf <- {
  list.files() ->.;
  .[grepl("GSOCseq", .) & grepl(".tif$", .)]
}

lapply(lf, function(x) {
  gdalUtilities::gdalwarp(srcfile = x, dstfile = gsub(".tif", "_5070.tif", x), t_srs = "EPSG:5070")
})

lf <- {
  list.files() ->.;
  .[grepl("GSOCseq", .) & grepl("5070.tif$", .)]
}

lf_s <- lapply(lf, function(x) {
  strsplit(x, "_|Map030_5070.tif")[[1]][5:6]
})
lf_s <- do.call("rbind", lf_s)

df <- data.frame(
  lf, 
  vars = lf_s[, 1], 
  cond = lf_s[, 2], 
  unc  = ifelse(grepl("Uncertainty", lf), "Uncertainty", "")
)
df <- within(df, {
  cond = ifelse(cond == "Uncertainty", vars, cond)
  vars = ifelse(vars == cond, "Uncertainty", vars)
  comb = paste(vars, unc)
  comb = gsub(" $", "", comb)
  # comb = ifelse(comb %in% c("T0", "BAU", "SSM"), "cond", comb)
  unit = factor(
    comb, 
    levels = c("T0", "finalSOC", "AbsDiff", "RelDiff", "ASR", "RSR", "ASR Uncertainty", "RSR Uncertainty", "Uncertainty Uncertainty"),
    labels = c("t/ha", "t/ha", "t/ha", "t/ha", "t/ha/yr", "t/ha/yr", "%", "%", "%")
  )
  label = factor(
    comb, 
    levels = c("T0", "finalSOC", "AbsDiff", "RelDiff", "ASR", "RSR", "ASR Uncertainty", "RSR Uncertainty", "Uncertainty Uncertainty"),
    labels = c("Intial Soil Organic Carbon", "Soil Organic Carbon", "Absolute Difference", "Relative Difference", "Absolute Sequestration Rate", "Relative Sequestration Rate", "Absolute Sequestration Rate Uncertainty", "Relative Sequestration Rate Uncertainty", "Uncertainty")
  )
  brks = NA
  brks = ifelse(unit == "t/ha/yr", data.frame(x = c(0, 0.05, 0.1, 0.15, 0.2)), brks)
  brks = ifelse(unit == "t/ha", data.frame(x = c(0, 15, 30, 50, 800)), brks)
  brks = ifelse(unit == "t/ha" & comb %in% c("AbsDiff", "RelDiff"), data.frame(x = c(0, 0.2, 0.5, 1, 170)), brks)
  brks = ifelse(unit == "%", data.frame(x = c(0, 25, 40, 60, 100)), brks)
  
})


# load rasters into a list
l_rs <- {
  split(df, df$comb) ->.;
  lapply(., function(x) {
    rs <- rast(x$lf)
    names(rs) <- x$cond
    return(rs)
  }) ->.;
}


# plot
# idx <- c(13, 27, 29)
# by(df[idx, ], df$comb[idx], function(x) {
by(df, df$comb, function(x) {
  
  cat(paste("plotting", x$comb[1]), "\n")
      
  tm_p <- tm_shape(l_rs[[x$comb[1]]], raster.downsample = FALSE) + 
    tm_raster(
      title = x$unit[1], 
      palette = "viridis", 
      # style = "quantile", n = 4
      breaks = x$brks[[1]]
    ) + 
    tm_shape(aoi) + tm_borders(col = "black", lwd = 0.2) +
    tm_grid(projection = 4326) +
    tm_scale_bar(position = c("left", "bottom")) +
    tm_layout(
      main.title = x$label[1],
      panel.labels = x$cond, 
      legend.format = list(digits = ifelse(x$unit[1] == "%", 0, 2))
    )
  
  png(paste0("plots_", x$comb[1], ".png"), units = "in", width = 12, height = 6, res = 300)
  print(tm_p)
  dev.off()
})



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

