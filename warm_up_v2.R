library(SoilR)
library(sf)
library(raster)
library(rgdal)
library(soilassessment)


# Set working directory
setwd("D:/geodata/project_data/gsp-gsocseq")


# Load warm up data
wu_df <- readRDS(file = "wu_pts.rds")
wu_df <- wu_df[complete.cases(wu_df), ]
wu_sf <- st_as_sf(
  wu_df,
  coords = c("x", "y"),
  crs    = 4326
  ) %>%
  st_transform(crs = 5070)
wu_df <- cbind(round(st_coordinates(wu_sf)), st_drop_geometry(wu_sf))
wu_df$xy <- paste0(wu_df$X, wu_df$Y, sep = "_")

# Load spin up results
su_sf <- readRDS(file = "conus_su_results.rds")
su_sf <- st_transform(su_sf, crs = 5070)
su_df <- cbind(round(st_coordinates(su_sf)), st_drop_geometry(su_sf))
su_df$xy <- paste0(su_df$X, su_df$Y, sep = "_")


# combine
wu_df <- merge(wu_df, su_df, by = "xy", all.x = TRUE)



# Extract the layers from the Vector
SOC_r    <- wu_df$CONUS_GSOCmap1.5.0
SOC_min  <- wu_df$SOC * 0.8
SOC_max  <- wu_df$SOC * 1.2

pClay_r   <- wu_df$CONUS_gnatsgo_fy20_1km_clay_wt
pClay_min <- pClay_r * 0.9
pClay_max <- pClay_r * 1.1

DR   <- wu_df$CONUS_glc_shv10_DOM_DR
LU   <- wu_df$CONUS_glc_shv10_DOM
TEMP <- wu_df[grepl("_tmmx$", names(wu_df))]
PREC <- wu_df[grepl("_pr$", names(wu_df))]
PET  <- wu_df[grepl("_pet$",  names(wu_df))]
COV  <- wu_df[grepl("^CON_",  names(wu_df))]
NPP  <- wu_df$`CONUS_NPP_MIAMI_MEAN_81-00_AOI`
NPP_min <- wu_df$`CONUS_NPP_MIAMI_MEAN_81-00_AOI_MIN`
NPP_max <- wu_df$`CONUS_NPP_MIAMI_MEAN_81-00_AOI_MAX`


#Apply NPP coeficientes
NPP <- (
   LU == 2 | LU == 12 | LU == 13) * NPP * 0.53 + 
  (LU == 4) * NPP * 0.88 + 
  (LU == 3 | LU == 5 | LU == 6 | LU == 8) * NPP * 0.72
NPP_min <- (
   LU == 2 | LU == 12 | LU == 13) * NPP_min * 0.53 + 
  (LU == 4) * NPP_min * 0.88 + 
  (LU == 3 | LU == 5 | LU == 6 | LU == 8) * NPP_min * 0.72
NPP_max <- (
  LU == 2 | LU == 12 | LU == 13) * NPP_max * 0.53 + 
  (LU == 4) * NPP_max * 0.88 + 
  (LU == 3 | LU == 5 | LU == 6 | LU == 8) * NPP_max * 0.72


# Define Years for Cinputs calculations ----
year <- seq(1 / 12, 1, by = 1 / 12)
nSim <- ncol(TEMP) / 12

# Paddy fields coefficent  ----
# fPR = 0.4 if the target point is class = 13 , else fPR=1
# From Shirato and Yukozawa 2004
fPR <- (LU == 13) * 0.4 + (LU != 13) * 1


# Temperature effects per month ----
fT <- as.data.frame(lapply(TEMP, function(x) fT.RothC(x)))


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

fW_r   <- fW(pClay_r,   PREC, PET, COV, s_thk = 30, pE = 1)
fW_min <- fW(pClay_min, PREC, PET, COV, s_thk = 30, pE = 1)
fW_max <- fW(pClay_max, PREC, PET, COV, s_thk = 30, pE = 1)


# Vegetation Cover effects ----
fC <- COV


# Set the factors frame for Model calculations ----
xi_r   <- fT
xi_min <- fT
xi_max <- fT

for (i in 1:ncol(fT)) {
  xi_r[, i] <- fT[, i] * fW_r[, i] * fC[, i] * fPR
}
for (i in 1:ncol(fT)) {
  xi_min[, i] <- fT[, i] * fW_min[, i] * fC[, i] * fPR
}
for (i in 1:ncol(fT)) {
  xi_max[, i] <- fT[, i] * fW_max[, i] * fC[, i] * fPR
}



