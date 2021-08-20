library(SoilR)
library(sf)
library(raster)
library(rgdal)
library(soilassessment)

rm(list = ls())

# Vectorized NPPmodel
source("D:/GIS/TOOLBOXES/gsp-gsocseq/functions.R")


# Set working directory
setwd("D:\\geodata\\fin_project_data4\\fin_project_data5\\gsp-gsocseq_pe")


# Load warm up data
wu_df <- data.table::fread(file = "wu_pts_v2.csv")
# wu_sf <- st_as_sf(
#   wu_df,
#   coords = c("x", "y"),
#   crs    = 4326
#   ) %>%
#   st_transform(crs = 5070)
# wu_df    <- cbind(round(st_coordinates(wu_sf)), st_drop_geometry(wu_sf))
# wu_df    <- wu_df[complete.cases(wu_df), ]
# wu_df$xy <- paste(wu_df$X, wu_df$Y, sep = "_")
# rm(wu_sf)
# 
# idx <- 1:nrow(wu_df)
# idx <- as.integer(cut(idx, breaks = quantile(idx, probs = seq(0, 1, 0.2)), include.lowest = TRUE))
# wu_df <- wu_df[idx == 2, ]


# Load spin up results
su_df <- readRDS(file = "conus_su_results_v3_analytical.rds")
# su_df <- st_transform(su_df, crs = 5070)
# su_df <- cbind(round(st_coordinates(su_df)), st_drop_geometry(su_df))
# su_df$xy <- paste(su_df$X, su_df$Y, sep = "_")
# saveRDS(su_df, "conus_su_results_v3.rds")
# su_df <- readRDS("conus_su_results_v3.rds")
# su_vars <- names(su_df)


# combine
# wu_df <- merge(wu_df, su_df, by = "xy", all.x = TRUE)
# saveRDS(wu_df, "wu_df.rds")
# 
# vars <- c("X.y", "Y.y", "ID.y", "CONUS_gnatsgo_fy20_1km_clay_wt", "CONUS_glc_shv10_DOM_DR", su_vars[c(-1, -2, -3, -50)])
# su_df <- wu_df[vars]
# saveRDS(su_df, "su_df.rds")
# rm(su_df)
# wu_df <- readRDS("wu_df.rds")
wu_df <- cbind(su_df, wu_df)
rm(su_df)


# Extract variables
wu_df <- as.data.frame(wu_df)

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
# NPP_rv <- (
#   LU == 2 | LU == 12 | LU == 13)           * NPP_rv  * 0.53 + 
#   (LU == 4)                                * NPP_rv  * 0.88 + 
#   (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP_rv  * 0.72
# NPP_min <- (
#   LU == 2 | LU == 12 | LU == 13)           * NPP_min * 0.53 + 
#   (LU == 4)                                * NPP_min * 0.88 + 
#   (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP_min * 0.72
# NPP_max <- (
#   LU == 2  | LU == 12 | LU == 13)          * NPP_max * 0.53 + 
#   (LU == 4)                                * NPP_max * 0.88 + 
#   (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP_max * 0.72


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
  (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP * 0.72
}

NPP_M     <- sapply(1:19, function(i) npp_coef(LU, NPP_M[, i]))
NPP_M_min <- sapply(1:19, function(i) npp_coef(LU, NPP_M_min[, i]))
NPP_M_max <- sapply(1:19, function(i) npp_coef(LU, NPP_M_max[, i]))


C_rv  <- cbind(replicate(19, wu_df$Cin.r))
C_min <- cbind(replicate(19, wu_df$Cin.min))
C_max <- cbind(replicate(19, wu_df$Cin.max))

# C_rv  <- apply(C_rv,  2, function(x) ifelse(is.na(x), 0, x))
# C_min <- apply(C_min, 2, function(x) ifelse(is.na(x), 0, x))
# C_max <- apply(C_max, 2, function(x) ifelse(is.na(x), 0, x))

# this line
C_rv[, 1]  <- C_rv[, 1]  / NPP_rv     * NPP_M[, 1]
C_min[, 1] <- C_min[, 1] / NPP_min * NPP_M_min[, 1]
C_max[, 1] <- C_max[, 1] / NPP_max * NPP_M_max[, 1]


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
rm(PREC, PREC_sum, PREC_sum_min, PREC_sum_max)
rm(TEMP, TEMP_avg, TEMP_avg_min, TEMP_avg_max)
rm(NPP_M, NPP_M_min, NPP_M_max, NPP_rv, NPP_min, NPP_max)
rm(C_rv, C_min, C_max)



# Moisture effects per month ----
xi <- lapply(1:10, function(x) {
  
  fn <- paste0("wu_pts_sub_", x, "_v2.rds")
  
  cat("computing ", fn , as.character(Sys.time()), "\n")
  
  wu_df <- readRDS(file = fn)
  pClay_r <- wu_df$CONUS_gnatsgo_fy20_1km_clay_wt
  TEMP <- wu_df[grepl("_tmmx$", names(wu_df))]
  PREC <- wu_df[grepl("_pr$", names(wu_df))]
  PET  <- wu_df[grepl("_pet$",  names(wu_df))]
  COV  <- wu_df[grepl("^CON_",  names(wu_df))]
  LU   <- wu_df$CONUS_glc_shv10_DOM
  
  
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
  
  
  saveRDS(xi_r,   paste0("wu_effcts_r_sub_",   x, ".rds"))
  saveRDS(xi_min, paste0("wu_effcts_min_sub_", x, ".rds"))
  saveRDS(xi_max, paste0("wu_effcts_max_sub_", x, ".rds"))
  
})


xi_r <- do.call(
  "rbind", 
  lapply(1:10, function(x){
    temp <- readRDS(paste0("wu_effcts_r_sub_", x, ".rds"))
  })
)

xi_min <- do.call(
  "rbind", 
  lapply(1:10, function(x){
    temp <- readRDS(paste0("wu_effcts_min_sub_", x, ".rds"))
  })
)

xi_max <- do.call(
  "rbind", 
  lapply(1:10, function(x){
    temp <- readRDS(paste0("wu_effcts_max_sub_", x, ".rds"))
  })
)

data.table::fwrite(xi_r,   "wu_effcts_r.csv")
data.table::fwrite(xi_min, "wu_effcts_min.csv")
data.table::fwrite(xi_max, "wu_effcts_max.csv")


# RUN THE MODEL from soilassessment ----
# Roth C soilassesment in parallel
# load inputs ----
nSim  <- raster::nlayers(raster::stack("CONUS_Temp_Stack_228_01_19_TC.tif")) / 12
years <- seq(1 / 12, 1, by = 1 / 12)

su_df <- as.data.frame(readRDS("conus_su_results_v3_analytical.rds"))

# pClay_r   <- su_df$CONUS_gnatsgo_fy20_1km_clay_wt
# pClay_min <- su_df$CONUS_gnatsgo_fy20_1km_clay_wt * 0.9
# pClay_max <- su_df$CONUS_gnatsgo_fy20_1km_clay_wt * 1.1

pClay_r   <- su_df$pClay.r
pClay_min <- su_df$pClay.min * 0.9
pClay_max <- su_df$pClay.max* 1.1

xi_r   <- as.data.frame(data.table::fread("wu_effcts_r.csv"))
xi_min <- as.data.frame(data.table::fread("wu_effcts_min.csv"))
xi_max <- as.data.frame(data.table::fread("wu_effcts_max.csv"))

C_rv  <- readRDS("wu_C_rv.rds")
C_min <- readRDS("wu_C_min.rds")
C_max <- readRDS("wu_C_max.rds")


library(parallel)

# C input equilibrium. (Ceq) ----
# su_df2 <- su_df[1:2000, ]; C_rv2 <- C_rv[1:2000, ]; xi_r2 <- xi_r[1:2000, ]

clus <- makeCluster(4)
clusterExport(clus, list("su_df", "C_rv", "xi_r", "years", "carbonTurnover", "rothC_wu"))

Sys.time()
rothC_r <- rothC_wu(
  time = years,
  su_df = su_df, pClay_var = "pClay.r", C0_vars = c("fract.dpm.r", "fract.rpm.r", "fract.bio.r", "fract.hum.r", "fract.iom.r"), DR_var = "DR", 
  C_m   = C_rv,
  xi_df = xi_r
)
Sys.time()
saveRDS(rothC_r, file = "rothC_r_wu_v2.rds")
stopCluster(clus)

# ==============================================================================


# rerun on negative values using a different solver
rc_wu <- as.data.frame(do.call("rbind", readRDS(file = "rothC_r_wu_v2.rds")))
idx <- which(apply(rc_wu, 1, function(x) any(x < 0)))

su_df <- as.data.frame(readRDS("conus_su_results_v3_analytical.rds"))
su_df <- su_df[idx, ]

C_rv  <- readRDS("wu_C_rv.rds")[idx, ]
# xi_r  <- readRDS("wu_effcts_r.rds")[idx, ]
x <- read.csv("wu_effcts_r.csv")
xi_r <- x[idx, ]
rm(x)

library(parallel)
clus <- makeCluster(4)
clusterExport(clus, list("idx", "su_df", "C_rv", "xi_r", "years", "RothCModel", "getC", "rothC_wu"))

Sys.time()
rothC_r_nn <- rothC_wu_nn(
  time  = years,
  su_df = su_df, pClay_var = "pClay.r", C0_vars = c("fract.dpm.r", "fract.rpm.r", "fract.bio.r", "fract.hum.r", "fract.iom.r"), DR_var = "DR", 
  C_m   = C_rv,
  xi_df = xi_r
)
Sys.time()
saveRDS(rothC_r_nn, file = "rothC_r_wu_nonneg.rds")
stopCluster(clus)

# ==============================================================================

# Cmin input equilibrium. (Ceq) ----
clus <- makeCluster(4)
clusterExport(clus, list("su_df", "C_min", "xi_min", "years", "carbonTurnover", "rothC_wu"))

su_df <- as.data.frame(readRDS("conus_su_results_v3_analytical.rds"))

Sys.time()
rothC_min <- rothC_wu(
  time = years,
  su_df = su_df, pClay_var = "pClay.min",  C0_vars = c("fract.dpm.min", "fract.rpm.min", "fract.bio.min", "fract.hum.min", "fract.iom.min"), DR_var = "DR", 
  C_m = C_min,
  xi_df = xi_min
)
Sys.time()
saveRDS(rothC_min, file = "rothC_min_wu.rds")
stopCluster(clus)

# ==============================================================================

# rerun on negative values using a different solver
rc_wu_min <- as.data.frame(do.call("rbind", readRDS(file = "rothC_min_wu.rds")))
idx <- which(apply(rc_wu_min, 1, function(x) any(x < 0)))

su_df  <- readRDS("conus_su_results_v3_analytical.rds")[idx, ]
C_min  <- readRDS("wu_C_min.rds")[idx, ]
xi_min <- as.data.frame(data.table::fread("wu_effcts_min.csv"))[idx, ]


library(parallel)
clus <- makeCluster(4)
clusterExport(clus, list("idx", "su_df", "C_min", "xi_min", "years", "RothCModel", "getC", "rothC_wu_nn"))

Sys.time()
rothC_min_nn <- rothC_wu_nn(
  time = years,
  su_df = su_df, pClay_var = "pClay.min", C0_vars = c("fract.dpm.min", "fract.rpm.min", "fract.bio.min", "fract.hum.min", "fract.iom.min"), DR_var = "DR", 
  C_m = C_min,
  xi_df = xi_min
)
Sys.time()
saveRDS(rothC_min_nn, file = "rothC_min_wu_nonneg_v2.rds")
stopCluster(clus)



# ==============================================================================

# Cmax input equilibrium. (Ceq) ----

su_df  <- as.data.frame(readRDS("conus_su_results_v3_analytical.rds"))

clus <- makeCluster(4)
clusterExport(clus, list("su_df", "C_max", "xi_max", "years", "carbonTurnover", "rothC_wu"))

Sys.time()
rothC_max <- rothC_wu(
  time = years,
  su_df = su_df, pClay_var = "pClay.max", C0_vars = c("fract.dpm.max", "fract.rpm.max", "fract.bio.max", "fract.hum.max", "fract.iom.max"), DR_var = "DR", 
  C_m = C_max,
  xi_df = xi_max
)
Sys.time()
saveRDS(rothC_max, file = "rothC_max_wu_v2.rds")
stopCluster(clus)


# ==============================================================================

# rerun on negative values using a different solver
rc_max_wu <- as.data.frame(do.call("rbind", readRDS(file = "rothC_max_wu_v2.rds")))
idx <- which(apply(rc_max_wu, 1, function(x) any(x < 0)))


# reload the following data from above if necessary
su_df  <- as.data.frame(readRDS("conus_su_results_v3_analytical.rds"))
su_df <- su_df[idx, ]

C_max  <- readRDS("wu_C_max.rds")[idx, ]
xi_max <- as.data.frame(data.table::fread("wu_effcts_max.csv"))[idx, ]


library(parallel)
clus <- makeCluster(4)
clusterExport(clus, list("idx", "su_df", "C_max", "xi_max", "years", "RothCModel", "getC", "rothC_wu_nn"))

Sys.time()
rothC_max_nn <- rothC_wu_nn(
  time  = years,
  su_df = su_df, pClay_var = "pClay.max", C0_vars = c("fract.dpm.max", "fract.rpm.max", "fract.bio.max", "fract.hum.max", "fract.iom.max"), DR_var = "DR",
  C_m   = C_max,
  xi_df = xi_max
)
Sys.time()
saveRDS(rothC_max_nn, file = "rothC_max_wu_nonneg.rds")
stopCluster(clus)


# ==============================================================================
# load rothC outputs, C inputs and spin-up results ----
rc_wu_r    <- as.data.frame(do.call("rbind", readRDS(file = "rothC_r_wu_v2.rds")))
rc_wu_r_nn <- as.data.frame(do.call("rbind", readRDS(file = "rothC_r_wu_nonneg.rds")))

rc_wu_min     <- as.data.frame(do.call("rbind", readRDS(file = "rothC_min_wu.rds")))
rc_wu_min_nn  <- as.data.frame(do.call("rbind", readRDS(file = "rothC_min_wu_nonneg_v2.rds")))

rc_wu_max     <- as.data.frame(do.call("rbind", readRDS(file = "rothC_max_wu_v2.rds")))
rc_wu_max_nn  <- as.data.frame(do.call("rbind", readRDS(file = "rothC_max_wu_nonneg.rds")))

C_rv  <- readRDS("wu_C_rv.rds")
C_min <- readRDS("wu_C_min.rds")
C_max <- readRDS("wu_C_max.rds")

# su_df    <- readRDS("su_sdf_v2.rds")
su_df  <- as.data.frame(readRDS("conus_su_results_v3_analytical.rds"))


# replace
idx <- which(apply(rc_wu_r, 1,   function(x) any(x < 0)))
rc_wu_r[idx, -1] <- rc_wu_r_nn

# rc_wu_min <- cbind(19, rc_wu_min); names(rc_wu_min) <- paste0("V", 1:6) # fix
idx <- which(apply(rc_wu_min, 1, function(x) any(x < 0)))
rc_wu_min[idx, -1] <- rc_wu_min_nn

idx <- which(apply(rc_wu_max, 1, function(x) any(x < 0)))
rc_wu_max[idx, -1] <- rc_wu_max_nn


# check
format(summary(rc_wu_r), big.mark = ",", scientific = FALSE)
sum(apply(rc_wu_r, 1, function(x) any(x < 0)))


# combine
# names(su_df)[1:3] <- c("X", "Y", "ID")
vars <- c("x", "y", "id", "SOC.r", "Cin.r")

rc_wu_all <- rbind(
  cbind(source = "r",   id = 1:nrow(rc_wu_r),   Cinput = C_rv[19],  CinputFORWARD = rowMeans(C_rv),  rc_wu_r),
  cbind(source = "min", id = 1:nrow(rc_wu_min), Cinput = C_min[19], CinputFORWARD = rowMeans(C_min), rc_wu_min),
  cbind(source = "max", id = 1:nrow(rc_wu_max), Cinput = C_max[19], CinputFORWARD = rowMeans(C_max), rc_wu_max)
)
rc_wu_all$SOC_t0 <- rowSums(rc_wu_all[5:10])
names(rc_wu_all)[5:10] <- c("time", "DPM_wu", "RPM_wu", "BIO_wu", "HUM_wu", "IOM_wu")

rc_wu_all <- reshape(rc_wu_all, direction = "wide",
                     idvar   = "id",
                     timevar = "source", 
                     v.names = names(rc_wu_all[c(3, 4, 6:11)])
)
rc_wu_all <- cbind(su_df[vars], rc_wu_all)


# save final results
saveRDS(rc_wu_all, file = "rothC_r_wu_final.rds")

