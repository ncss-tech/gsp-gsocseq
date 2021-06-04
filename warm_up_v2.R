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
# wu_df <- readRDS(file = "wu_pts.rds")
# # wu_df <- wu_df[complete.cases(wu_df), ]
# wu_sf <- st_as_sf(
#   wu_df,
#   coords = c("x", "y"),
#   crs    = 4326
#   ) %>%
#   st_transform(crs = 5070)
# wu_df <- cbind(round(st_coordinates(wu_sf)), st_drop_geometry(wu_sf))
# wu_df$xy <- paste0(wu_df$X, wu_df$Y, sep = "_")
# 
# # idx <- 1:nrow(wu_df)
# # idx <- as.integer(cut(idx, breaks = quantile(idx, probs = seq(0, 1, 0.2)), include.lowest = TRUE))
# # wu_df <- wu_df[idx == 2, ]
# rm(wu_sf)

# Load spin up results
# su_sf <- readRDS(file = "conus_su_results.rds")
# su_vars <- names(su_sf)
# su_sf <- st_transform(su_sf, crs = 5070)
# su_df <- cbind(round(st_coordinates(su_sf)), st_drop_geometry(su_sf))
# su_df$xy <- paste0(su_df$X, su_df$Y, sep = "_")
# saveRDS(su_df, "conus_su_results_v2.rds")
# su_df <- readRDS("conus_su_results_v2.rds")
# rm(su_sf)

# combine
# wu_df <- merge(wu_df, su_df, by = "xy", all.x = TRUE)
# rm(su_df)
# saveRDS(wu_df, "wu_df.rds")
# 
# vars <- c("X.y", "Y.y", "ID.y", "CONUS_gnatsgo_fy20_1km_clay_wt", "CONUS_glc_shv10_DOM_DR", su_vars[c(-1, -48)])
# saveRDS(wu_df[vars], "su_df.rds")

wu_df <- readRDS("wu_df.rds")


# Extract the layers from the Vector
LU   <- wu_df$CONUS_glc_shv10_DOM
TEMP <- wu_df[grepl("_tmmx$", names(wu_df))]
PREC <- wu_df[grepl("_pr$", names(wu_df))]
NPP_rv  <- wu_df$`CONUS_NPP_MIAMI_MEAN_81-00_AOI`
NPP_min <- wu_df$`CONUS_NPP_MIAMI_MEAN_81-00_AOI_MIN`
NPP_max <- wu_df$`CONUS_NPP_MIAMI_MEAN_81-00_AOI_MAX`


# Calculate Cinputs
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


# Apply NPP coeficientes
NPP_rv <- (
  LU == 2 | LU == 12 | LU == 13)           * NPP_rv  * 0.53 + 
  (LU == 4)                                * NPP_rv  * 0.88 + 
  (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP_rv  * 0.72
NPP_min <- (
  LU == 2 | LU == 12 | LU == 13)           * NPP_min * 0.53 + 
  (LU == 4)                                * NPP_min * 0.88 + 
  (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP_min * 0.72
NPP_max <- (
  LU == 2  | LU == 12 | LU == 13)          * NPP_max * 0.53 + 
  (LU == 4)                                * NPP_max * 0.88 + 
  (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP_max * 0.72


NPP_M     <- sapply(1:ncol(PREC_sum), function(i) {
  NPPmodel(PREC_sum[, i]    , TEMP_avg[, i]    , "miami") * (1/100) * 0.5
})
NPP_M_min <- sapply(1:ncol(PREC_sum), function(i) {
  NPPmodel(PREC_sum_min[, i], TEMP_avg_min[, i], "miami") * (1/100) * 0.5
})
NPP_M_max <- sapply(1:ncol(PREC_sum), function(i) {
  NPPmodel(PREC_sum_max[, i], TEMP_avg_max[, i], "miami") * (1/100) * 0.5
})


npp_coef <- function(LU, NPP) {
  (LU == 2 | LU == 12 | LU == 13)          * NPP * 0.53 +
  (LU == 4)                                * NPP * 0.88 + 
  (LU ==3  | LU == 5  | LU == 6 | LU == 8) * NPP * 0.72
}

NPP_M     <- sapply(1:19, function(i) npp_coef(LU, NPP_M[, i]))
NPP_M_min <- sapply(1:19, function(i) npp_coef(LU, NPP_M_min[, i]))
NPP_M_max <- sapply(1:19, function(i) npp_coef(LU, NPP_M_max[, i]))


C_rv  <- cbind(replicate(19, wu_df$Ceq.r))
C_min <- cbind(replicate(19, wu_df$Ceq.min))
C_max <- cbind(replicate(19, wu_df$Ceq.max))


C_rv[, 1]  <- C_rv[, 1]  / NPP_M[, 1]     * NPP_M[, 1]
C_min[, 1] <- C_min[, 1] / NPP_M_min[, 1] * NPP_M_min[, 1]
C_max[, 1] <- C_max[, 1] / NPP_M_max[, 1] * NPP_M_max[, 1]


idx <- 2:19
C_rv[, idx]  <- sapply(idx, function(i) C_rv[,  i - 1] / NPP_M[, (i - 1)]   * NPP_M[, i])
C_min[, idx] <- sapply(idx, function(i) C_min[, i - 1] / NPP_M_min[, i - 1] * NPP_M_min[, i])
C_max[, idx] <- sapply(idx, function(i) C_max[, i - 1] / NPP_M_max[, i - 1] * NPP_M_max[, i])


idx <- 1:19
C_rv[, idx]  <- sapply(idx, function(i) ifelse(is.na(C_rv[,  i]),  0,  C_rv[,  i]))
C_min[, idx] <- sapply(idx, function(i) ifelse(is.na(C_min[,  i]), 0,  C_min[,  i]))
C_max[, idx] <- sapply(idx, function(i) ifelse(is.na(C_max[,  i]), 0,  C_max[,  i]))


saveRDS(C_rv,  "wu_C_rv.rds")
saveRDS(C_min, "wu_C_min.rds")
saveRDS(C_max, "wu_C_max.rds")


rm(wu_df)
rm(PREC_sum, PREC_sum_min, PREC_sum_max)
rm(TEMP_avg, TEMP_avg_min, TEMP_avg_max)
rm(NPP_M, NPP_M_min, NPP_M_max)
rm(C_rv, C_min, C_max)



# Moisture effects per month ----
fW <- function(pClay, PREC, PET, COV, s_thk = 30, pE = 1) {
  
  M     <- PREC - PET * pE
  
  n    <- ncol(M)  / 12
  idx  <- rep(1:12, n)
  bare  <- as.data.frame(lapply(COV[, idx], function(x) x > 0.8))
  
  B <- as.data.frame(lapply(bare, function(x) ifelse(x == FALSE, 1, 1.8)))
  Max.TSMD <- -(20 + 1.3 * pClay - 0.01 * (pClay ^ 2)) * (s_thk / 23)
  Max.TSMD <- as.data.frame(lapply(1 / B, function(x) x * Max.TSMD))
  rm(B, bare)
  
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


wu_df <- readRDS("wu_df.rds")
pClay_r   <- wu_df$CONUS_gnatsgo_fy20_1km_clay_wt
TEMP <- wu_df[grepl("_tmmx$", names(wu_df))]
PREC <- wu_df[grepl("_pr$", names(wu_df))]
PET  <- wu_df[grepl("_pet$",  names(wu_df))]
COV  <- wu_df[grepl("^CON_",  names(wu_df))]
rm(wu_df)


fW_r   <- fW(pClay_r      , PREC,        PET, COV, s_thk = 30, pE = 1)
fW_min <- fW(pClay_r * 0.9, PREC * 0.95, PET, COV, s_thk = 30, pE = 1)
fW_max <- fW(pClay_r * 1.1, PREC * 1.05, PET, COV, s_thk = 30, pE = 1)


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


# Vegetation Cover effects ----
fC <- COV[, rep(1:12, length.out = ncol(PREC))]


# Paddy fields coefficent  ----
# fPR = 0.4 if the target point is class = 13 , else fPR=1
# From Shirato and Yukozawa 2004
fPR <- (LU == 13) * 0.4 + (LU != 13) * 1


# Set the factors frame for Model calculations ----
xi_r   <- fT_r
xi_min <- fT_min
xi_max <- fT_max

for (i in 1:ncol(fT_r)) {
  xi_r[, i]   <- fT_r[, i]  * fW_r[, i]    * fC[, i] * fPR
}
for (i in 1:ncol(fT_min)) {
  xi_min[, i] <- fT_min[, i] * fW_min[, i] * fC[, i] * fPR
}
for (i in 1:ncol(fT_max)) {
  xi_max[, i] <- fT_max[, i] * fW_max[, i] * fC[, i] * fPR
}


xi_r   <- cbind(id = 1:nrow(xi_r),   xi_r)
xi_min <- cbind(id = 1:nrow(xi_min), xi_min)
xi_max <- cbind(id = 1:nrow(xi_max), xi_max)


saveRDS(xi_r,   "wu_effcts_r.rds");   rm(xi_r,   fW_r,   fT_r,   fC, fPR)
saveRDS(xi_min, "wu_effcts_min.rds"); rm(xi_min, fW_min, fT_min, fC, fPR)
saveRDS(xi_max, "wu_effcts_max.rds"); rm(xi_max, fW_max, fT_max, fC, fPR)



# RUN THE MODEL from soilassessment ----
# Roth C soilassesment in parallel
nSim  <- raster::nlayers(raster::stack("CONUS_Temp_Stack_228_01_19_TC.tif")) / 12
years <- seq(1 / 12, 1, by = 1 / 12)

su_df <- readRDS("su_df.rds")

xi_r   <- readRDS("wu_effcts_r.rds")
xi_min <- readRDS("wu_effcts_min.rds")
xi_max <- readRDS("wu_effcts_max.rds")

C_rv  <- readRDS("wu_C_rv.rds")
C_min <- readRDS("wu_C_min.rds")
C_max <- readRDS("wu_C_max.rds")

DR   <- su_df$CONUS_glc_shv10_DOM_DR


library(parallel)

clus <- makeCluster(7)

# C input equilibrium. (Ceq) ----
clusterExport(clus, list("su_df", "C_rv", "DR", "xi_r", "years", "carbonTurnover"))

Sys.time()
rothC_r <- parLapply(clus, 1:nrow(su_df), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(su_df[i, ]$DPM_p.r, su_df[i, ]$RPM_p.r, su_df[i, ]$BIO_p.r, su_df[i, ]$HUM_p.r, su_df[i, ]$FallIOM.r),
    In   = C_rv[i, 1],
    Dr   = DR[i],
    clay = su_df[i, ]$CONUS_gnatsgo_fy20_1km_clay_wt,
    effcts = data.frame(years, rep(unlist(xi_r[i, 2:229]), length.out = length(years))),
    solver = "euler"
  )
  fp <- list(tail(temp, 1)[-1])
  
  for (j in 2:19) {
    temp <- carbonTurnover(
      tt   = years,
      C0   = c(fp[[1]][1], fp[[1]][2], fp[[1]][3], fp[[1]][4], fp[[1]][5]),
      In   = C_rv[i, j],
      Dr   = DR[i],
      clay = su_df[i, ]$CONUS_gnatsgo_fy20_1km_clay_wt,
      effcts = data.frame(years, rep(unlist(xi_r[i, 2:229]), length.out = length(years))),
      solver = "euler"
    )
    fp[[1]] <- tail(temp, 1)[-1]
  }
  
  fp <- unlist(fp)
  
  return(fp)
})
Sys.time()
saveRDS(rothC_r, file = "rothC_r_wu.rds")
stopCluster(clus)

rc_wu <- do.call("rbind", rothC_r)

