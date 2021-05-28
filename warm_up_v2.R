library(SoilR)
library(sf)
library(raster)
library(rgdal)
library(soilassessment)


# Vectorized NPPmodel
source("C:/workspace2/github/ncss-tech/gsp-gsocseq/NPPmodel.R")


# Set working directory
setwd("D:/geodata/project_data/gsp-gsocseq")


# Load warm up data
wu_df <- readRDS(file = "wu_pts.rds")
# wu_df <- wu_df[complete.cases(wu_df), ]
wu_sf <- st_as_sf(
  wu_df,
  coords = c("x", "y"),
  crs    = 4326
  ) %>%
  st_transform(crs = 5070)
wu_df <- cbind(round(st_coordinates(wu_sf)), st_drop_geometry(wu_sf))
wu_df$xy <- paste0(wu_df$X, wu_df$Y, sep = "_")

idx <- 1:nrow(wu_df)
idx <- as.integer(cut(idx, breaks = quantile(idx, probs = seq(0, 1, 0.2)), include.lowest = TRUE))
wu_df <- wu_df[idx == 2, ]


# Load spin up results
su_sf <- readRDS(file = "conus_su_results.rds")
su_sf <- st_transform(su_sf, crs = 5070)
su_df <- cbind(round(st_coordinates(su_sf)), st_drop_geometry(su_sf))
su_df$xy <- paste0(su_df$X, su_df$Y, sep = "_")


# combine
wu_df <- merge(wu_df, su_df, by = "xy", all.x = TRUE)
rm(wu_sf, su_sf, su_df)


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
NPP_rv  <- wu_df$`CONUS_NPP_MIAMI_MEAN_81-00_AOI`
NPP_min <- wu_df$`CONUS_NPP_MIAMI_MEAN_81-00_AOI_MIN`
NPP_max <- wu_df$`CONUS_NPP_MIAMI_MEAN_81-00_AOI_MAX`


# Apply NPP coeficientes
NPP_rv <- (
   LU == 2 | LU == 12 | LU == 13)          * NPP_rv    * 0.53 + 
  (LU == 4) * NPP_rv * 0.88 + 
  (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP_rv    * 0.72
NPP_min <- (
   LU == 2 | LU == 12 | LU == 13)          * NPP_min * 0.53 + 
  (LU == 4) * NPP_min * 0.88 + 
  (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP_min * 0.72
NPP_max <- (
  LU == 2  | LU == 12 | LU == 13)          * NPP_max * 0.53 + 
  (LU == 4) * NPP_max * 0.88 + 
  (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP_max * 0.72


# Cinputs
yrs0119  <- lapply(2001:2019, function(x) paste0(x, formatC(1:12, width = 2, flag = "0")))

TEMP_avg <- sapply(yrs0119, function(x) {
  x2   <- paste0(x, collapse = "|")
  idx  <- which(grepl(x2, names(TEMP)))
  rowMeans(TEMP[idx])
})
# TEMP_avg     <- rowMeans(TEMP)
TEMP_avg_min <- TEMP_avg * 1.02
TEMP_avg_max <- TEMP_avg * 0.98

PREC_sum <- sapply(yrs0119, function(x) {
  x2   <- paste0(x, collapse = "|")
  idx  <- which(grepl(x2, names(PREC)))
  rowSums(PREC[idx])
})
# PREC_sum     <- rowSums(PREC)
PREC_sum_min <- PREC_sum * 0.95
PREC_sum_max <- PREC_sum * 1.05

NPP_M     <- sapply(1:ncol(PREC_sum), function(i) {
  NPPmodel(PREC_sum[, i]    , TEMP_avg[, i]    , "miami") * (1/100) * 0.5
})
NPP_M_min <- sapply(1:ncol(PREC_sum), function(i) {
  NPPmodel(PREC_sum_min[, i], TEMP_avg_min[, i], "miami") * (1/100) * 0.5
})
NPP_M_max <- sapply(1:ncol(PREC_sum), function(i) {
  NPPmodel(PREC_sum_max[, i], TEMP_avg_max[, i], "miami") * (1/100) * 0.5
})

npp <- rbind(
  cbind(source = "rv",  npp = as.data.frame(NPP_M)),
  cbind(source = "min", npp = as.data.frame(NPP_M_min)),
  cbind(source = "max", npp = as.data.frame(NPP_M_max))
)
names(npp)[-1] <- paste0("X", 2001:2019) 
lu <- c(LU, LU, LU)

npp_coef <- function(LU, NPP) {
  (LU == 2 | LU == 12 | LU == 13) * NPP * 0.53 +
  (LU == 4) * NPP * 0.88 + 
  (LU ==3  | LU == 5  | LU == 6 | LU == 8) * NPP * 0.72
}

npp[2:20] <- lapply(2:20, function(i) {
  npp_coef(lu, npp[, i])
})

NPP_M     <- as.matrix(subset(npp, source == "rv",  select = - source))
NPP_M_min <- as.matrix(subset(npp, source == "min", select = - source))
NPP_M_max <- as.matrix(subset(npp, source == "max", select = - source))

C_rv  <- cbind(replicate(19, wu_df$Ceq.r))
C_min <- cbind(replicate(19, wu_df$Ceq.min))
C_max <- cbind(replicate(19, wu_df$Ceq.max))

C_rv  <- C_rv[, 1]  / NPP_M[, 1]     * NPP_M[, 1]
C_min <- C_min[, 1] / NPP_M_min[, 1] * NPP_M_min[, 1]
C_max <- C_max[, 1] / NPP_M_max[, 1] * NPP_M_max[, 1]

idx <- 2:19
C_rv[, idx]  <- sapply(idx, function(i) C_rv[,  i - 1] / NPP_M[, (i - 1)]     * NPP_M[, i])
C_min[, idx] <- sapply(idx, function(i) C_min[, i - 1] / NPP_M_min[, i - 1] * NPP_M_min[, i])
C_max[, idx] <- sapply(idx, function(i) C_max[, i - 1] / NPP_M_max[, i - 1] * NPP_M_max[, i])



# Define Years for Cinputs calculations ----
year <- seq(1 / 12, 1, by = 1 / 12)
nSim <- nlayers(stack("CONUS_Temp_Stack_228_01_19_TC.tif")) / 12


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



