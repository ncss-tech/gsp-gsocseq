
# SPATIAL SOIL R  for VECTORS

###### SPIN UP ################

# MSc Ing Agr Luciano E Di Paolo
# Dr Ing Agr Guillermo E Peralta
###################################
# SOilR from Sierra, C.A., M. Mueller, S.E. Trumbore (2012).
#Models of soil organic matter decomposition: the SoilR package, version 1.0 Geoscientific Model Development, 5(4),
#1045--1060. URL http://www.geosci-model-dev.net/5/1045/2012/gmd-5-1045-2012.html.
#####################################

# Setup ----
library(SoilR)
library(raster)
library(rgdal)
library(soilassessment)


source('C:/Users/stephen.roecker/OneDrive - USDA/projects/gsp-gsocseq/code/functions.R')


# set to aoi
setwd("D:/geodata/project_data/gsp-gsocseq/AK")


# Load inputs ----

# Stack_Set_1 is a stack that contains the spatial variables 
su_df <- readRDS(file = "su_df.RDS")


# Extract the layers from the Vector
SOC_r    <- su_df$SOC
SOC_min  <- su_df$SOC * 0.8
SOC_max  <- su_df$SOC * 1.2

pClay_r   <- su_df$CLAY
pClay_min <- su_df$CLAY * 0.9
pClay_max <- su_df$CLAY * 1.1

DR   <- su_df$DR
LU   <- su_df$LU
TEMP <- su_df[grepl("^TEMP_", names(su_df))]
PREC <- su_df[grepl("^PREC_", names(su_df))]
PET  <- su_df[grepl("^PET_",  names(su_df))]
COV  <- su_df[grepl("^COV_",  names(su_df))]


# Paddy fields coefficent  ----
# fPR = 0.4 if the target point is class = 13 , else fPR=1
# From Shirato and Yukozawa 2004
fPR <- (LU == 13) * 0.4 + (LU != 13) * 1


# Temperature effects per month ----
fT_r   <- as.data.frame(lapply(TEMP,        function(x) {
  temp <- fT.RothC(x)
  temp <- ifelse(is.na(temp), 0, temp)
}))
fT_min <- as.data.frame(lapply(TEMP * 1.02, function(x) {
  temp <- fT.RothC(x)
  temp <- ifelse(is.na(temp), 0, temp)
}))
fT_max <- as.data.frame(lapply(TEMP * 0.98, function(x) {
  temp <- fT.RothC(x)
  temp <- ifelse(is.na(temp), 0, temp)
}))

    
# Moisture effects per month ----
fW <- function(pClay, PREC, PET, COV, s_thk = 30, pE = 1) {
  
  M     <- PREC - PET * pE
  bare  <- as.data.frame(lapply(COV, function(x) x > 0.8))
  
  B <- as.data.frame(lapply(bare, function(x) ifelse(x == FALSE, 1, 1.8)))
  Max.TSMD <- -(20 + 1.3 * pClay - 0.01 * (pClay ^ 2)) * (s_thk / 23)
  Max.TSMD <- as.data.frame(lapply(1 / B, function(x) x * Max.TSMD))
  
  Acc.TSMD <- M
  Acc.TSMD[, 1] <- ifelse(M[, 1] > 0, 0, M[, 1])
  
  for (i in 2:ncol(M)) {
    Acc.TSMD_i <- Acc.TSMD[, i - 1] + M[, i]
    Acc.TSMD[, i] <- ifelse(
      Acc.TSMD_i < 0,
      Acc.TSMD_i,
      0
      )
    Acc.TSMD[, i] <- ifelse(
      Acc.TSMD[, i] <= Max.TSMD[, i],
      Max.TSMD[, i],
      Acc.TSMD[, i]
      )
  }
   
  fW <- Acc.TSMD
    
  for (i in 1:ncol(fW)) {
    fW[, i] <- ifelse(
      Acc.TSMD[, i] > 0.444 * Max.TSMD[, i], 
      1, 
      (0.2 + 0.8 * ((Max.TSMD[, i] - Acc.TSMD[, i]) / (Max.TSMD[, i] - 0.444 * Max.TSMD[, i])))
      )
    fW[, i] <- raster::clamp(fW[, i], lower = 0.2)
    }
    
  return(fW)
  }

fW_r   <- fW(pClay_r      , PREC,        PET, COV, s_thk = 30, pE = 1)
fW_min <- fW(pClay_r * 0.9, PREC * 0.95, PET, COV, s_thk = 30, pE = 1)
fW_max <- fW(pClay_r * 1.1, PREC * 1.05, PET, COV, s_thk = 30, pE = 1)



# Vegetation Cover effects ----
fC <- COV
    

# Set the factors frame for Model calculations ----
xi_r   <- fT_r
xi_min <- fT_min
xi_max <- fT_max

for (i in 1:ncol(fT_r)) {
  xi_r[, i] <- fT_r[, i]       * fW_r[, i]   * fC[, i] * fPR
}
for (i in 1:ncol(fT_min)) {
  xi_min[, i] <- fT_min[, i] * fW_min[, i] * fC[, i] * fPR
}
for (i in 1:ncol(fT_max)) {
  xi_max[, i] <- fT_max[, i] * fW_max[, i] * fC[, i] * fPR
}

saveRDS(xi_r,   "su_effcts_r.rds")
saveRDS(xi_min, "su_effcts_min.rds")
saveRDS(xi_max, "su_effcts_max.rds")



# # Roth C outputs
# ro <- list(data.frame(id = NULL, C1 = NULL, C2 = NULL, C3 = NULL, C4 = NULL, C5 = NULL))[rep(1, nrow(ri))]
# 
# # RUN THE MODEL from soilassessment
# # Roth C soilassesment
# 
# iter <- as.integer(seq(200, nrow(ri), 200))
# for (i in 1:nrow(ri)) {
#    if (any(i %in% iter)) cat("iteration", i, as.character(Sys.time()), "\n")
#   temp <- carbonTurnover(
#     tt   = years,
#     C0   = c(DPMptf, RPMptf, BIOptf, HUMptf, ri[i, 3]),
#     In   = Cinputs,
#     Dr   = ri[i, 4],
#     clay = ri[i, 5],
#     effcts = data.frame(years, rep(unlist(ri[i, 5:16]), length.out = length(years))),
#     solver = "euler"
#     )
#  ro[[i]] <- c(ri[1], unlist(temp[6000, 2:6]))
# }


# RUN THE MODEL from soilassessment ----
# Roth C soilassesment in parallel
su_df <- readRDS(file = "su_df.RDS")

years <- seq(1 / 12, 500, by = 1 / 12)

Cinputs <- 1
DR      <- su_df$DR

pClay_r   <- su_df$CLAY
pClay_min <- su_df$CLAY * 0.9
pClay_max <- su_df$CLAY * 1.1

SOC_r    <- su_df$SOC
SOC_min  <- su_df$SOC * 0.8
SOC_max  <- su_df$SOC * 1.2

FallIOM_r   <- 0.049 * SOC_r^(1.139)
FallIOM_min <- 0.049 * SOC_min^(1.139)
FallIOM_max <- 0.049 * SOC_max^(1.139)

fractI <- cbind(
  DPM = DR / (DR + 1), 
  RPM = 1 - DR / (DR + 1), 
  HUM = 0
  )

xi_r   <- readRDS("su_effcts_r.rds")
xi_min <- readRDS("su_effcts_min.rds")
xi_max <- readRDS("su_effcts_max.rds")

xi_r2   <- rowMeans(xi_r)
xi_min2 <- rowMeans(xi_min)
xi_max2 <- rowMeans(xi_max)


library(parallel)

clus <- makeCluster(15)

# C input equilibrium. (Ceq) ----
# clusterExport(clus, list("DR", "pClay_r", "ri_r", "FallIOM_r", "years", "carbonTurnover")) # , "rothC"))
# 
# Sys.time()
# rothC_r <- parLapply(clus, 1:nrow(su_df), function(i) {
#   
#   temp <- carbonTurnover(
#     tt   = years,
#     C0   = c(0, 0, 0, 0, FallIOM_r[i]),
#     In   = 1,
#     Dr   = DR[i],
#     clay = pClay_r[i],
#     effcts = data.frame(years, rep(unlist(ri_r[i, 1:12]), length.out = length(years))),
#     solver = "euler"
#   )
#   temp <- tail(temp, 1)
#   
#   return(temp)
# })
# Sys.time()
# # saveRDS(rothC_r, file = "rothC_r_v3.rds")
# stopCluster(clus)

clus <- makeCluster(15)
clusterExport(clus, list("SOC_r", "pClay_r", "xi_r2", "fractI", "fget_equilibrium_fractions.RothC_input", "fIOM.Falloon.RothC"))

Sys.time()
rothC_r <- parLapply(clus, 1:nrow(su_df), function(i) {
  
  temp <- fget_equilibrium_fractions.RothC_input(
    xi     = xi_r2[i], 
    C.tot  = SOC_r[i], 
    clay   = pClay_r[i], 
    fractI = fractI[i, ]
  )
  
  return(temp)
})
Sys.time()
# saveRDS(rothC_r, file = "rothC_r_v3_analytical.rds")
stopCluster(clus)



# # UNCERTAINTIES C input equilibrium (MINIMUM)
# # C input equilibrium. (Ceq) ----
# clusterExport(clus, list("DR", "pClay_min", "ri_min", "FallIOM_min", "years", "carbonTurnover")) # , "rothC"))
# 
# Sys.time()
# rothC_min <- parLapply(clus, 1:nrow(su_df), function(i) {
#   
#   temp <- carbonTurnover(
#     tt   = years,
#     C0   = c(0, 0, 0, 0, FallIOM_min[i]),
#     In   = 1,
#     Dr   = DR[i],
#     clay = pClay_min[i],
#     effcts = data.frame(years, rep(unlist(ri_min[i, 1:12]), length.out = length(years))),
#     solver = "euler"
#   )
#   temp <- tail(temp, 1)
#   
#   return(temp)
# })
# Sys.time()
# # saveRDS(rothC_min, file = "rothC_min_v2.rds")
# stopCluster(clus)


clus <- makeCluster(15)
clusterExport(clus, list("SOC_min", "pClay_min", "xi_min2", "fractI", "fget_equilibrium_fractions.RothC_input", "fIOM.Falloon.RothC"))

Sys.time()
rothC_min <- parLapply(clus, 1:nrow(su_df), function(i) {
  
  temp <- fget_equilibrium_fractions.RothC_input(
    xi     = xi_min2[i], 
    C.tot  = SOC_min[i], 
    clay   = pClay_min[i], 
    fractI = fractI[i, ]
  )
  
  return(temp)
})
Sys.time()
# saveRDS(rothC_min, file = "rothC_min_v3_analytical.rds")
stopCluster(clus)


# UNCERTAINTIES C input equilibrium (MAXIMUM)
# clusterExport(clus, list("su_df", "DR", "pClay_max", "ri_max", "FallIOM_max", "years", "carbonTurnover"))
# 
# Sys.time()
# rothC_max <- parLapply(clus, 1:nrow(su_df), function(i) {
#   
#   temp <- carbonTurnover(
#     tt   = years,
#     C0   = c(0, 0, 0, 0, FallIOM_max[i]),
#     In   = 1,
#     Dr   = DR[i],
#     clay = pClay_max[i],
#     effcts = data.frame(years, rep(unlist(ri_max[i, 1:12]), length.out = length(years))),
#     solver = "euler"
#   )
#   temp <- tail(temp, 1)
#   
#   return(temp)
# })
# Sys.time()
# # saveRDS(rothC_max, file = "rothC_max_v2.rds")
# stopCluster(clus)


clus <- makeCluster(15)
clusterExport(clus, list("SOC_max", "pClay_max", "xi_max2", "fractI", "fget_equilibrium_fractions.RothC_input", "fIOM.Falloon.RothC"))

Sys.time()
rothC_max <- parLapply(clus, 1:nrow(su_df), function(i) {
  
  temp <- fget_equilibrium_fractions.RothC_input(
    xi     = xi_max2[i], 
    C.tot  = SOC_max[i], 
    clay   = pClay_max[i], 
    fractI = fractI[i, ]
  )
  
  return(temp)
})
Sys.time()
# saveRDS(rothC_max, file = "rothC_max_v3_analytical.rds")
stopCluster(clus)



# Outputs ----
rothC_r <- as.data.frame(
  cbind(su_df[c("aoi", "cell", "x", "y", "LU", "DR")], SOC = SOC_r, FallIOM = FallIOM_r, pClay = pClay_r,
        source = "r",
  do.call(
    "rbind", 
    readRDS("rothC_r_v3_analytical.rds")
    )))

rothC_min <- as.data.frame(
  cbind(su_df[c("aoi", "cell", "x", "y", "LU", "DR")], SOC = SOC_min, FallIOM = FallIOM_min, pClay = pClay_min, 
        source = "min",
  do.call(
    "rbind", 
    readRDS("rothC_min_v3_analytical.rds")
  )))

rothC_max <- as.data.frame(
  cbind(su_df[c("aoi", "cell", "x", "y", "LU", "DR")], SOC = SOC_max, FallIOM = FallIOM_max,  pClay = pClay_max, 
        source = "max",
  do.call(
    "rbind", 
    readRDS("rothC_max_v3_analytical.rds")
  )))

rothC_df <- rbind(rothC_r, rothC_min, rothC_max)

rothC_df$fract.sum <- rowSums(rothC_df[grepl("^fract", names(rothC_df))])

# rothC_df <- within(rothC_df, {
#   fb_t = C1 + C2 + C3 + C4 + C5
#   m    = (fb_t - FallIOM) / 1
#   Ceq  = (SOC  - FallIOM) / m
# })
# 
# 
# # crops, tree crops, and rice
# idx <- ifelse(rothC_r$LU %in% c(2, 12, 13), TRUE, FALSE)
# rothC_crops <- within(rothC_df[idx, ], {
#   RPM_p = ((0.184  * SOC + 0.1555) * (pClay + 1.275) ^ (-0.1158)) * 0.9902 +
#     0.4788
#   BIO_p = ((0.014  * SOC + 0.0075) * (pClay + 8.8473) ^ (0.0567)) * 1.09038 +
#     0.04055
#   HUM_p = ((0.7148 * SOC + 0.5069) * (pClay + 0.3421) ^ (0.0184)) * 0.9878 -
#     0.3818
#   DPM_p = SOC - FallIOM - RPM_p - HUM_p - BIO_p
#   feq_t = RPM_p + BIO_p + HUM_p + DPM_p + FallIOM
# })
# 
# 
# # trees
# idx <- ifelse(rothC_r$LU == 4, TRUE, FALSE)
# rothC_trees <- within(rothC_df[idx, ], {
#   RPM_p = ((0.184  * SOC + 0.1555) * (pClay + 1.275) ^ (-0.1158)) * 1.7631 + 0.4043
#   BIO_p = ((0.014  * SOC + 0.0075) * (pClay + 8.8473) ^ (0.0567)) * 0.9757 + 0.0209
#   HUM_p = ((0.7148 * SOC + 0.5069) * (pClay + 0.3421) ^ (0.0184)) * 0.8712 - 0.2904
#   DPM_p = SOC - FallIOM - RPM_p - HUM_p - BIO_p
#   feq_t = RPM_p + BIO_p + HUM_p + DPM_p + FallIOM
# })
# 
# 
# # grass and shrubs
# idx <- ifelse(rothC_df$LU %in% c(3, 5, 6, 8), TRUE, FALSE)
# rothC_grasses <- within(rothC_df[idx, ], {
#   RPM_p = ((0.184  * SOC + 0.1555) * (pClay + 1.275) ^ (-0.1158)) * 1.3837 + 0.4692
#   BIO_p = ((0.014  * SOC + 0.0075) * (pClay + 8.8473) ^ (0.0567)) * 1.03401 + 0.02531
#   HUM_p = ((0.7148 * SOC + 0.5069) * (pClay + 0.3421) ^ (0.0184)) * 0.9316 - 0.5243
#   DPM_p = SOC - FallIOM - RPM_p - HUM_p - BIO_p
#   feq_t = RPM_p + BIO_p + HUM_p + DPM_p + FallIOM
# })
# 
# 
# # combine
# rothC_df <- rbind(rothC_crops, rothC_trees, rothC_grasses)
# rothC_df[c("time", "m")] <- NULL


# reshape
# nm   <- names(rothC_df)
# vars <- c("id", "x", "y", "source", "LU")
# vars2 <- nm[! nm %in% vars]
# 
# system.time(rothC_dfw <- reshape(rothC_df, direction = "wide",
#                      idvar = c("id"),
#                      timevar = "source", 
#                      v.names = vars2
#                      ))

library(data.table)

rothC_dfw2 <- dcast(
  as.data.table(rothC_df), 
  cell + x + y + aoi + LU + DR ~ source, 
  value.var = c("SOC", "FallIOM", "pClay", "fract.dpm", "fract.rpm", "fract.bio", "fract.hum", "fract.iom", "Cin", "fract.sum"),
  sep = "."
  )
 
# # convert to sf
# rothC_sf <- st_as_sf(
#   rothC_dfw,
#   coords = c("x", "y"),
#   crs = 4326
# )
rc_df <- as.data.frame(rothC_dfw2)
saveRDS(rc_df, file = "su_results_v3_analytical.rds")


# inspect outputs
rc_sf <- st_as_sf(
  test,
  coords = c("X", "Y"),
  crs    = 4326
)
rc_sf <- rc_sf[rc_sf$aoi == "AK1", ]
names(rc_sf) <- gsub("\\.", "_", names(rc_sf))
write_sf(rc_sf, dsn = "test.gpkg", driver = "GPKG", overwrite = TRUE)

aoi <- "AK1"
gsoc <- rast(paste0(aoi, "_GSOCmap1.5.0.tif"))
gsoc[!is.na(gsoc)] <- 1
lu   <- rast(paste0(aoi, "_glc_shv10_DOM.tif"))
lu   <- lu %in% c(2, 3, 5, 12, 13)
gsoc <- gsoc * lu
writeRaster(gsoc, "test_socfracsum.tif", overwrite = TRUE)

gdalUtilities::gdal_rasterize(
  src_datasource = "test.gpkg",
  a              = "SOC_t0_r",
  dst_filename   = "test_socfracsum.tif",
  of             = "GTiff",
  te             = bbox(gsoc),
  tr             = res(gsoc),
  co             = c("COMPRESS=DEFLATE"),
  a_srs          = "EPSG:4326",
  a_nodata       = -999
)
