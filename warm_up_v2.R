
# Setup ----
library(SoilR)
library(sf)
library(raster)
library(rgdal)
library(soilassessment)


source('C:/Users/stephen.roecker/OneDrive - USDA/projects/gsp-gsocseq/code/functions.R')


setwd("D:/geodata/project_data/gsp-gsocseq/CONUS")


# Load inputs ----

# warm up data
wu_df <- data.table::fread(file = "wu_pts.csv", data.table = FALSE)
wu_df <- as.data.frame(wu_df)

wu_df <- readRDS("wu_df.rds")
wu_df <- wu_df[order(wu_df$cell), ]


# Load spin up results
su_df <- readRDS("su_results_v3_analytical.rds")

all(su_df$cell == wu_df$cell)
wu_df <- cbind(su_df[-1], wu_df)
# rm(su_df)

idx <- complete.cases(wu_df)
sum(idx)
wu_df <- wu_df[idx, ]
su_df <- su_df[idx, ]

vv <- readRDS("Vector_points.rds")

wu_df <- wu_df[wu_df$cell %in% vv$cell, ]
su_df <- su_df[su_df$cell %in% vv$cell, ]




# Extract variables
LU   <- wu_df$LU
TEMP <- wu_df[grepl("_tmmx$", names(wu_df))]
PREC <- wu_df[grepl("_pr$",   names(wu_df))]
NPP_rv  <- wu_df$`NPP_MIAMI_MEAN_81-00_AOI`
NPP_min <- wu_df$`NPP_MIAMI_MEAN_81-00_AOI_MIN`
NPP_max <- wu_df$`NPP_MIAMI_MEAN_81-00_AOI_MAX`



# Calculate Cinputs
yrs0119  <- lapply(2001:2019, function(x) paste0(x, formatC(1:12, width = 2, flag = "0")))

TEMP_avg <- sapply(yrs0119, function(x) {
  x2   <- paste0(x, collapse = "|")
  idx  <- which(grepl(x2, names(TEMP)))
  rowMeans(TEMP[idx])
})
TEMP_avg_min <- TEMP_avg * 1.02
TEMP_avg_max <- TEMP_avg * 0.98

PREC_sum <- sapply(yrs0119, function(x) {
  x2   <- paste0(x, collapse = "|")
  idx  <- which(grepl(x2, names(PREC)))
  rowSums(PREC[idx])
})
PREC_sum_min <- PREC_sum * 0.95
PREC_sum_max <- PREC_sum * 1.05


# Apply NPP coeficientes
npp_coef <- function(LU, NPP) {
  (LU == 2 | LU == 12 | LU == 13)          * NPP * 0.53 +
  (LU == 4)                                * NPP * 0.88 + 
  (LU == 3 | LU == 5  | LU == 6 | LU == 8) * NPP * 0.72
}

NPP_rv <- npp_coef(LU, NPP_rv)
NPP_min <- npp_coef(LU, NPP_min)
NPP_max <- npp_coef(LU, NPP_max)


NPP_M     <- sapply(1:ncol(PREC_sum), function(i) {
  NPPmodel(PREC_sum[, i]    , TEMP_avg[, i]    , "miami") * (1/100) * 0.5
})
NPP_M_min <- sapply(1:ncol(PREC_sum), function(i) {
  NPPmodel(PREC_sum_min[, i], TEMP_avg_min[, i], "miami") * (1/100) * 0.5
})
NPP_M_max <- sapply(1:ncol(PREC_sum), function(i) {
  NPPmodel(PREC_sum_max[, i], TEMP_avg_max[, i], "miami") * (1/100) * 0.5
})


NPP_M     <- sapply(1:19, function(i) npp_coef(LU, NPP_M[, i]))
NPP_M_min <- sapply(1:19, function(i) npp_coef(LU, NPP_M_min[, i]))
NPP_M_max <- sapply(1:19, function(i) npp_coef(LU, NPP_M_max[, i]))


C_rv  <- cbind(replicate(19, wu_df$Cin.r))
C_min <- cbind(replicate(19, wu_df$Cin.min))
C_max <- cbind(replicate(19, wu_df$Cin.max))

C_rv[, 1]  <- C_rv[, 1]  / NPP_rv  * NPP_M[, 1]
C_min[, 1] <- C_min[, 1] / NPP_min * NPP_M_min[, 1]
C_max[, 1] <- C_max[, 1] / NPP_max * NPP_M_max[, 1]


idx <- 2:19
C_rv[, idx]  <- sapply(idx, function(i) C_rv[,  i - 1] / NPP_M[, (i - 1)]   * NPP_M[, i])
C_min[, idx] <- sapply(idx, function(i) C_min[, i - 1] / NPP_M_min[, i - 1] * NPP_M_min[, i])
C_max[, idx] <- sapply(idx, function(i) C_max[, i - 1] / NPP_M_max[, i - 1] * NPP_M_max[, i])


saveRDS(C_rv,  "wu_C_rv.rds")
saveRDS(C_min, "wu_C_min.rds")
saveRDS(C_max, "wu_C_max.rds")


# rm(wu_df)
rm(PREC, PREC_sum, PREC_sum_min, PREC_sum_max)
rm(TEMP, TEMP_avg, TEMP_avg_min, TEMP_avg_max)
rm(NPP_M, NPP_M_min, NPP_M_max, NPP_rv, NPP_min, NPP_max)
rm(C_rv, C_min, C_max)



# Moisture effects per month ----
# if file sizes are too large run the loop, otherwise run the lines inside the loop
xi <- lapply(1:10, function(x) {
  
  fn <- paste0("wu_pts_sub_", x, ".rds")
  
  cat("computing ", fn , as.character(Sys.time()), "\n")
  
  wu_df <- readRDS(file = fn)
  
  pClay_r <- wu_df$CONUS_gnatsgo_fy20_1km_clay_wt
  TEMP <- wu_df[grepl("_tmmx$", names(wu_df))]
  PREC <- wu_df[grepl("_pr$", names(wu_df))]
  PET  <- wu_df[grepl("_pet$",  names(wu_df))]
  COV  <- wu_df[grepl("^COV_",  names(wu_df))]
  LU   <- wu_df$LU
  
  
  # loop over each year separately
  ncols  <- ncol(TEMP)
  nyrs    <- ncols / 12 
  df_idx <- data.frame(
    yr = rep(1:nyrs, each = 12),
    cn = 1:ncols
    )
  
  fW_r <- TEMP
  for (i in 1:nyrs) {
    idx <- df_idx[df_idx$yr == i, ]$cn
    fW_r[, idx]   <- fW(pClay_r      , PREC[, idx],        PET[, idx], COV, s_thk = 30, pE = 1)
  }
  
  fW_min <- TEMP
  for (i in 1:nyrs) {
    idx <- df_idx[df_idx$yr == i, ]$cn
    fW_min[, idx] <- fW(pClay_r * 0.9, PREC[, idx] * 0.95, PET[, idx], COV, s_thk = 30, pE = 1)
  }
  
  fW_max <- TEMP
  for (i in 1:nyrs) {
    idx <- df_idx[df_idx$yr == i, ]$cn
    fW_max[, idx] <- fW(pClay_r * 1.1, PREC[, idx] * 1.05, PET[, idx], COV, s_thk = 30, pE = 1)
  }
  
  
  
  # Temperature effects per month ----
  fT_r   <- as.data.frame(lapply(TEMP,        function(x) {
    temp <- fT.RothC(x)
   }))
  fT_min <- as.data.frame(lapply(TEMP * 1.02, function(x) {
    temp <- fT.RothC(x)
   }))
  fT_max <- as.data.frame(lapply(TEMP * 0.98, function(x) {
    temp <- fT.RothC(x)
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
  # saveRDS(xi_r,   "wu_effcts_r.rds")
  # saveRDS(xi_min, "wu_effcts_min.rds")
  # saveRDS(xi_max, "wu_effcts_max.rds")

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
# Load Model Inputs ----

aoi <- "CONUS"
nSim  <- nlayers(stack(paste0(aoi, "_Temp_Stack_228_01-19_TC.tif"))) / 12
years <- seq(1 / 12, 1, by = 1 / 12)

aoi <- "CONUS"
su_df <- readRDS("su_results_v3_analytical.rds")

# xi_r   <- as.data.frame(data.table::fread("wu_effcts_r.csv"))
# xi_min <- as.data.frame(data.table::fread("wu_effcts_min.csv"))
# xi_max <- as.data.frame(data.table::fread("wu_effcts_max.csv"))
xi_r   <- readRDS("wu_effcts_r.rds")
xi_min <- readRDS("wu_effcts_min.rds")
xi_max <- readRDS("wu_effcts_max.rds")

C_rv  <- readRDS("wu_C_rv.rds")
C_min <- readRDS("wu_C_min.rds")
C_max <- readRDS("wu_C_max.rds")



# set.seed(42)
# su_df2 <- su_df;
# idx <- sample(1:nrow(su_df), 20000)
# su_df <- su_df2[idx, ]


# xi_r2 <- xi_r; C_rv2 <- C_rv
# xi_r  <- xi_r2[idx, ]
# C_rv  <- C_rv2[idx, ]
# 
# xi_min2 <- xi_min; C_min2 <- C_min
# xi_min <- xi_min2[idx, ]
# C_min  <- C_min2[idx, ]
# 
# xi_max2 <- xi_max; C_max2 <- C_max
# xi_max <- xi_max2[idx, ]
# C_max  <- C_max2[idx, ]
su_df <- readRDS("su_results_v3_analytical.rds")
vv <- readRDS("Vector_points.rds")
su_df <- su_df[su_df$cell %in% vv$cell, ]

idx <- complete.cases(su_df)
su_df <- su_df[idx, ]
C_rv  <- C_rv[idx, ]
xi_r  <- xi_r[idx, ]


# Run RothC ----
library(parallel)

## C input equilibrium. (Ceq) ----
# su_df2 <- su_df[1:2000, ]; C_rv2 <- C_rv[1:2000, ]; xi_r2 <- xi_r[1:2000, ]

clus <- makeCluster(16)
clusterExport(clus, list("su_df", "C_rv", "xi_r", "years", "carbonTurnover", "rothC_wu"))

Sys.time()
rothC_r <- rothC_wu(
  time = years,
  su_df = su_df, pClay_var = "pClay.r", C0_vars = c("fract.dpm.r", "fract.rpm.r", "fract.bio.r", "fract.hum.r", "fract.iom.r"), DR_var = "DR", 
  C_m   = C_rv,
  xi_df = xi_r
)
Sys.time()
# saveRDS(rothC_r, file = "rothC_r_wu.rds")
stopCluster(clus)



### rerun on negative values using a different solver ----
rc_wu <- as.data.frame(do.call("rbind", readRDS(file = "rothC_r_wu.rds")))
idx <- which(apply(rc_wu, 1, function(x) any(x < 0)))

# su_df <- readRDS("su_results_v3_analytical.rds")[idx, ]
# C_rv  <- readRDS("wu_C_rv.rds")[idx, ]
# xi_r  <- as.data.frame(data.table::fread("wu_effcts_r.csv")[idx, ])
# xi_r   <- readRDS("wu_effcts_r.rds")[idx, ]

su_df2 <- su_df[idx, ]
C_rv2  <- C_rv[idx, ]
xi_r2  <- xi_r[idx, ]


clus <- makeCluster(15)
clusterExport(clus, list("idx", "su_df2", "C_rv2", "xi_r2", "years", "RothCModel", "getC", "rothC_wu"))

Sys.time()
rothC_r_nn <- rothC_wu_nn(
  time  = years,
  su_df = su_df2, pClay_var = "pClay.r", C0_vars = c("fract.dpm.r", "fract.rpm.r", "fract.bio.r", "fract.hum.r", "fract.iom.r"), DR_var = "DR", 
  C_m   = C_rv2,
  xi_df = xi_r2
)
Sys.time()
# saveRDS(rothC_r_nn, file = "rothC_r_wu_nonneg.rds")
stopCluster(clus)



## Cmin input equilibrium. (Ceq) ----
clus <- makeCluster(16)
clusterExport(clus, list("su_df", "C_min", "xi_min", "years", "carbonTurnover", "rothC_wu"))

Sys.time()
rothC_min <- rothC_wu(
  time = years,
  su_df = su_df, pClay_var = "pClay.min", C0_vars = c("fract.dpm.min", "fract.rpm.min", "fract.bio.min", "fract.hum.min", "fract.iom.min"), DR_var = "DR", 
  C_m = C_min,
  xi_df = xi_min
)
Sys.time()
# saveRDS(rothC_min, file = "rothC_min_wu.rds")
stopCluster(clus)


### rerun on negative values using a different solver ----
rc_wu_min <- as.data.frame(do.call("rbind", readRDS(file = "rothC_min_wu.rds")))
idx <- which(apply(rc_wu_min, 1, function(x) any(x < 0)))

# C_min  <- readRDS("wu_C_min.rds")[idx, ]
# xi_min <- as.data.frame(data.table::fread("wu_effcts_min.csv")[idx, ])
# xi_min <- readRDS("wu_effcts_min.rds")[idx, ]


su_df2 <- su_df[idx, ]
C_min2  <- C_min[idx, ]
xi_min2  <- xi_min[idx, ]


clus <- makeCluster(16)
clusterExport(clus, list("idx", "su_df2", "C_min2", "xi_min2", "years", "RothCModel", "getC", "rothC_wu_nn"))

Sys.time()
rothC_min_nn <- rothC_wu_nn(
  time = years,
  su_df = su_df2, pClay_var = "pClay.min", C0_vars = c("fract.dpm.min", "fract.rpm.min", "fract.bio.min", "fract.hum.min", "fract.iom.min"), DR_var = "DR", 
  C_m = C_min2,
  xi_df = xi_min2
)
Sys.time()
# saveRDS(rothC_min_nn, file = "rothC_min_wu_nonneg.rds")
stopCluster(clus)



## Cmax input equilibrium. (Ceq) ----

rn <- as.integer(1:nrow(su_df))
p  <- quantile(rn, probs = seq(0, 1, 0.1))
idx <- as.integer(cut(rn, breaks = p, include.lowest = TRUE))
df <- data.frame(rn = 1:length(rn), idx)

Sys.time()
lapply(1:10, function(i){
  
  # subset to fit into memory
  idx <- df[df$idx == i, ]$rn
  su_df2  <<- su_df[idx, ]
  xi_max2 <<- xi_max[idx, ]
  C_max2  <<- C_max[idx, ]
  
  clus <- makeCluster(16)
  clusterExport(clus, list("su_df", "C_max", "xi_max", "years", "carbonTurnover", "rothC_wu"))
  
  Sys.time()
  rothC_max <- rothC_wu(
    time = years,
    su_df = su_df, pClay_var = "pClay.max", C0_vars = c("fract.dpm.max", "fract.rpm.max", "fract.bio.max", "fract.hum.max", "fract.iom.max"), DR_var = "DR", 
    C_m = C_max,
    xi_df = xi_max
  )
  Sys.time()
  saveRDS(rothC_max, file = paste0("rothC_max_wu_", i, ".rds"))
  # saveRDS(rothC_max, file = "rothC_max_wu.rds")
  stopCluster(clus)
})
Sys.time()




### rerun on negative values using a different solver ----
rc_wu_max <- as.data.frame(do.call("rbind", readRDS(file = "rothC_max_wu.rds")))
idx <- which(apply(rc_wu_max, 1, function(x) any(x < 0)))

# su_df  <- readRDS("su_results_v3_analytical.rds")[idx, ]
# C_max  <- readRDS("wu_C_max.rds")[idx, ]
# xi_max <- as.data.frame(data.table::fread("wu_effcts_max.csv")[idx, ])
# xi_max <- readRDS("wu_effcts_max.rds")[idx, ]

su_df2  <- su_df[idx, ]
C_max2  <- C_max[idx, ]
xi_max2 <- xi_max[idx, ]


clus <- makeCluster(16)
clusterExport(clus, list("idx", "su_df2", "C_max2", "xi_max2", "years", "RothCModel", "getC", "rothC_wu_nn"))

Sys.time()
rothC_max_nn <- rothC_wu_nn(
  time = years,
  su_df = su_df2, pClay_var = "pClay.max", C0_vars = c("fract.dpm.max", "fract.rpm.max", "fract.bio.max", "fract.hum.max", "fract.iom.max"), DR_var = "DR", 
  C_m = C_max2,
  xi_df = xi_max2
)
Sys.time()
# saveRDS(rothC_max_nn, file = "rothC_max_wu_nonneg.rds")
stopCluster(clus)



# Combine RothC Outputs, C inputs and Spin-Up Results ----
rc_wu_r    <- as.data.frame(do.call("rbind", readRDS(file = "rothC_r_wu.rds")))
rc_wu_r_nn <- as.data.frame(do.call("rbind", readRDS(file = "rothC_r_wu_nonneg.rds")))

rc_wu_min     <- as.data.frame(do.call("rbind", readRDS(file = "rothC_min_wu.rds")))
rc_wu_min_nn  <- as.data.frame(do.call("rbind", readRDS(file = "rothC_min_wu_nonneg.rds")))

rc_wu_max     <- as.data.frame(do.call("rbind", readRDS(file = "rothC_max_wu.rds")))
rc_wu_max_nn  <- as.data.frame(do.call("rbind", readRDS(file = "rothC_max_wu_nonneg.rds")))

su_df <- as.data.frame(readRDS(file = "su_results_v3_analytical.rds"))
idx <- complete.cases(su_df)
su_df <- su_df[idx, ]

C_rv  <- readRDS("wu_C_rv.rds")
C_min <- readRDS("wu_C_min.rds")
C_max <- readRDS("wu_C_max.rds")



# replace
idx <- which(apply(rc_wu_r, 1,   function(x) any(x < 0)))
rc_wu_r[idx, -1] <- rc_wu_r_nn

idx <- which(apply(rc_wu_min, 1, function(x) any(x < 0)))
rc_wu_min[idx, -1] <- rc_wu_min_nn

idx <- which(apply(rc_wu_max, 1, function(x) any(x < 0)))
rc_wu_max[idx, -1] <- rc_wu_max_nn


# check
format(summary(rc_wu_r), big.mark = ",", scientific = FALSE)
sum(apply(rc_wu_r, 1, function(x) any(x < 0)), na.rm = TRUE)


# combine
vars <- c("cell", "x", "y", "SOC.r", "Cin.r")

rc_wu_all <- rbind(
  cbind(source = "r",   cell = su_df$cell, Cinput = C_rv[19],  CinputFORWARD = rowMeans(C_rv),  rc_wu_r),
  cbind(source = "min", cell = su_df$cell, Cinput = C_min[19], CinputFORWARD = rowMeans(C_min), rc_wu_min),
  cbind(source = "max", cell = su_df$cell, Cinput = C_max[19], CinputFORWARD = rowMeans(C_max), rc_wu_max)
)
rc_wu_all$SOC_t0 <- rowSums(rc_wu_all[5:10])
names(rc_wu_all)[5:10] <- c("time", "DPM_wu", "RPM_wu", "BIO_wu", "HUM_wu", "IOM_wu")


rc_wu_all <- reshape(rc_wu_all, direction = "wide",
                     idvar   = "cell",
                     timevar = "source",
                     v.names = names(rc_wu_all[c(3:4, 6:11)])
)
all(su_df$cell == rc_wu_all$cell)
rc_wu_all$cell <- NULL
rc_wu_all <- cbind(su_df[vars], rc_wu_all)


# save final results
saveRDS(rc_wu_all, file = "rothC_r_wu_final.rds")

