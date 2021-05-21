
# SPATIAL SOIL R  for VECTORS

###### SPIN UP ################

# MSc Ing Agr Luciano E Di Paolo
# Dr Ing Agr Guillermo E Peralta
###################################
# SOilR from Sierra, C.A., M. Mueller, S.E. Trumbore (2012).
#Models of soil organic matter decomposition: the SoilR package, version 1.0 Geoscientific Model Development, 5(4),
#1045--1060. URL http://www.geosci-model-dev.net/5/1045/2012/gmd-5-1045-2012.html.
#####################################


library(SoilR)
library(raster)
library(rgdal)
library(soilassessment)


# Set working directory
setwd("D:/geodata/project_data/gsp-gsocseq")


# Stack_Set_1 is a stack that contains the spatial variables 
su_sdf <- readRDS(file = "su_sdf.RDS")
su_df  <- as.data.frame(su_sdf)
su_df  <- su_df[complete.cases(su_df), ]
su_df  <- cbind(ID = 1:nrow(su_df), su_df)


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


# Define Years for Cinputs calculations ----
years <- seq(1 / 12, 500, by = 1 / 12)


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

ri_r   <- cbind(id = 1:nrow(xi_r),   xi_r)
ri_min <- cbind(id = 1:nrow(xi_min), xi_min)
ri_max <- cbind(id = 1:nrow(xi_max), xi_max)

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

library(parallel)

clus <- makeCluster(15)


# C input equilibrium. (Ceq) ----
clusterExport(clus, list("wu_df", "ri_r", "years", "carbonTurnover")) # , "rothC"))

Sys.time()
rothC_r <- parLapply(clus, 1:200, function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(DPMptf, RPMptf, BIOptf, HUMptf, ri_r[i, 3]),
    In   = Cinputs,
    Dr   = ri_r[i, 4],
    clay = ri_r[i, 5],
    effcts = data.frame(years, rep(unlist(ri_r[i, 6:17]), length.out = length(years))),
    solver = "euler"
  )
  temp <- c(ri_r[i, 1], unlist(temp[6000, 2:6]))
  
  return(temp)
})
Sys.time()
# saveRDS(rothC_r, file = "rothC_r.rds")
stopCluster(clus)


# UNCERTAINTIES C input equilibrium (MINIMUM)
clusterExport(clus, list("ri_min", "years", "DPMptf", "RPMptf", "BIOptf", "HUMptf", "Cinputs", "carbonTurnover"))

Sys.time()
rothC_min <- parLapply(clus, 1:nrow(ri_min), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(DPMptf, RPMptf, BIOptf, HUMptf, ri_min[i, 3]),
    In   = Cinputs,
    Dr   = ri_min[i, 4],
    clay = ri_min[i, 5],
    effcts = data.frame(years, rep(unlist(ri_min[i, 6:17]), length.out = length(years))),
    solver = "euler"
  )
  temp <- c(ri_min[i, 1], unlist(temp[6000, 2:6]))
  
  return(temp)
})
Sys.time()
# saveRDS(rothC_min, file = "rothC_min.rds")
stopCluster(clus)


# UNCERTAINTIES C input equilibrium (MAXIMUM)
clusterExport(clus, list("ri_max", "years", "DPMptf", "RPMptf", "BIOptf", "HUMptf", "Cinputs", "carbonTurnover"))

Sys.time()
rothC_max <- parLapply(clus, 1:nrow(ri_max), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(DPMptf, RPMptf, BIOptf, HUMptf, ri_max[i, 3]),
    In   = Cinputs,
    Dr   = ri_max[i, 4],
    clay = ri_max[i, 5],
    effcts = data.frame(years, rep(unlist(ri_max[i, 6:17]), length.out = length(years))),
    solver = "euler"
  )
  temp <- c(ri_max[i, 1], unlist(temp[6000, 2:6]))
  
  return(temp)
})
Sys.time()
# saveRDS(rothC_max, file = "rothC_max.rds")
stopCluster(clus)




rothC_r <- as.data.frame(
  cbind(su_df[c("ID", "x", "y")], SOC = SOC_r, FallIOM = FallIOM_r, pClay = pClay_r, LU,
        source = "r",
  do.call(
    "rbind", 
    readRDS("rothC_r.rds")
    )))

rothC_min <- as.data.frame(
  cbind(su_df[c("ID", "x", "y")], SOC = SOC_min, FallIOM = FallIOM_min, pClay = pClay_min, LU, 
        source = "min",
  do.call(
    "rbind", 
    readRDS("rothC_min.rds")
  )))

rothC_max <- as.data.frame(
  cbind(su_df[c("ID", "x", "y")], SOC = SOC_max, FallIOM = FallIOM_max,  pClay = pClay_max, LU, 
        source = "max",
  do.call(
    "rbind", 
    readRDS("rothC_max.rds")
  )))

rothC_df <- rbind(rothC_r, rothC_min, rothC_max)

rothC_df <- within(rothC_df, {
  fb_t = C1 + C2 + C3 + C4 + C5
  m    = (fb_t - FallIOM) / 1
  Ceq  = (SOC  - FallIOM) / m
})


# crops, tree crops, and rice
idx <- ifelse(rothC_r$LU %in% c(2, 12, 13), TRUE, FALSE)
rothC_crops <- within(rothC_df[idx, ], {
  RPM_p = ((0.184  * SOC + 0.1555) * (pClay + 1.275) ^ (-0.1158)) * 0.9902 +
    0.4788
  BIO_p = ((0.014  * SOC + 0.0075) * (pClay + 8.8473) ^ (0.0567)) * 1.09038 +
    0.04055
  HUM_p = ((0.7148 * SOC + 0.5069) * (pClay + 0.3421) ^ (0.0184)) * 0.9878 -
    0.3818
  DPM_p = SOC - FallIOM - RPM_p - HUM_p - BIO_p
  feq_t = RPM_p + BIO_p + HUM_p + DPM_p + FallIOM
})


# trees
idx <- ifelse(rothC_r$LU == 4, TRUE, FALSE)
rothC_trees <- within(rothC_df[idx, ], {
  RPM_p = ((0.184  * SOC + 0.1555) * (pClay + 1.275) ^ (-0.1158)) * 1.7631 + 0.4043
  BIO_p = ((0.014  * SOC + 0.0075) * (pClay + 8.8473) ^ (0.0567)) * 0.9757 + 0.0209
  HUM_p = ((0.7148 * SOC + 0.5069) * (pClay + 0.3421) ^ (0.0184)) * 0.8712 - 0.2904
  DPM_p = SOC - FallIOM - RPM_p - HUM_p - BIO_p
  feq_t = RPM_p + BIO_p + HUM_p + DPM_p + FallIOM
})


# grass and shrubs
idx <- ifelse(rothC_df$LU %in% c(3, 5, 6, 8), TRUE, FALSE)
rothC_grasses <- within(rothC_df[idx, ], {
  RPM_p = ((0.184  * SOC + 0.1555) * (pClay + 1.275) ^ (-0.1158)) * 1.3837 + 0.4692
  BIO_p = ((0.014  * SOC + 0.0075) * (pClay + 8.8473) ^ (0.0567)) * 1.03401 + 0.02531
  HUM_p = ((0.7148 * SOC + 0.5069) * (pClay + 0.3421) ^ (0.0184)) * 0.9316 - 0.5243
  DPM_p = SOC - FallIOM - RPM_p - HUM_p - BIO_p
  feq_t = RPM_p + BIO_p + HUM_p + DPM_p + FallIOM
})


# combine
rothC_df <- rbind(rothC_crops, rothC_trees, rothC_grasses)
rothC_df[c("V1", "m")] <- NULL


# reshape
nm   <- names(rothC_df)
vars <- c("ID", "x", "y", "source", "LU")
vars2 <- nm[! nm %in% vars]

rothC_dfw <- reshape(rothC_df, direction = "wide",
                     idvar = c("ID"),
                     timevar = "source", 
                     v.names = vars2)
 
# convert to sf
rothC_sf <- st_as_sf(
  rothC_dfw,
  coords = c("x", "y"),
  crs = 4326
)

saveRDS(rothC_sf, file = "conus_su_results.rds")
