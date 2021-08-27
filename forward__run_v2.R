
# Setup ----
library(SoilR)
library(raster)
library(rgdal)
library(soilassessment)
library(sf)

source('C:/Users/stephen.roecker/OneDrive - USDA/projects/gsp-gsocseq/code/functions.R')

setwd("D:/geodata/project_data/gsp-gsocseq/CONUS")


# Load inputs ----
# forward stack
fr_df <- readRDS(file = "fr_df.RDS")
fr_df <- fr_df[order(fr_df$cell), ]


# warm up
wu_df <- readRDS(file = "rothC_r_wu_final.rds")


# combine
fr_df <- merge(fr_df, wu_df, by = "cell", all.y = TRUE)
fr_df <-fr_df[complete.cases(fr_df), ]


# Set the increase in Carbon input for each land use and each scenario
# Crops and Crop trees
# LU == 2 | LU == 12
Low_Crops  <- 1.05
Med_Crops  <- 1.10
High_Crops <- 1.2

# Shrublands, Grasslands , Herbaceous vegetation flooded & Sparse Vegetation
# LU == 3 | LU == 5 | LU == 6 | LU == 8
Low_Grass  <- 1.05
Med_Grass  <- 1.10
High_Grass <- 1.2

# Paddy Fields
# LU == 13
Low_PaddyFields  <- 1.05
Med_PaddyFields  <- 1.10
High_PaddyFields <- 1.2


# apply increase ----
lu_cl <- c(2, 12, 3, 5, 6, 8, 13)

fr_df <- within(fr_df, {
  CinputFORWARD_low     <- ifelse(LU %in% lu_cl, CinputFORWARD.r   * 1.05,        CinputFORWARD.r)
  CinputFORWARD_med     <- ifelse(LU %in% lu_cl, CinputFORWARD.r   * 1.10,        CinputFORWARD.r)
  CinputFORWARD_high    <- ifelse(LU %in% lu_cl, CinputFORWARD.r   * 1.20,        CinputFORWARD.r)
  
  CinputFORWARD_med_min <- ifelse(LU %in% lu_cl, CinputFORWARD.min * 1.05 - 0.15, CinputFORWARD.min)
  CinputFORWARD_med_max <- ifelse(LU %in% lu_cl, CinputFORWARD.max * 1.20 + 0.15, CinputFORWARD.max)
})

saveRDS(fr_df, file = "fr_df_Cinputs.rds")


# Extract variables
pClay_r <- fr_df$CLAY
TEMP    <- fr_df[grepl("TEMP_", names(fr_df))]
PREC    <- fr_df[grepl("PREC_", names(fr_df))]
PET     <- fr_df[grepl("PET_",  names(fr_df))]
COV     <- fr_df[grepl("COV_",  names(fr_df))]
LU      <- fr_df$LU


# Moisture effects per month ----
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
fC <- COV


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


saveRDS(xi_r,   "fr_effcts_r.rds")
saveRDS(xi_min, "fr_effcts_min.rds")
saveRDS(xi_max, "fr_effcts_max.rds")



# RUN THE MODEL from soilassessment ----
##  load inputs ----
years <- seq(1 / 12, 20, by = 1 / 12)

fr_df <- readRDS("fr_df_Cinputs.rds")
fr_df <- fr_df[!grepl("TEMP_|PREC_|PET_|COV_", names(fr_df))]
fr_df$CLAY_min <- fr_df$CLAY * 0.9
fr_df$CLAY_max <- fr_df$CLAY * 1.1

xi_r   <- readRDS("fr_effcts_r.rds")
xi_min <- readRDS("fr_effcts_min.rds")
xi_max <- readRDS("fr_effcts_max.rds")



# Roth C soilassesment in parallel
library(parallel)


## Cinput bau ----
clus <- makeCluster(16)
clusterExport(clus, list("fr_df", "xi_r", "years", "carbonTurnover"))

Sys.time()
rothC_bau <- parLapply(clus, 1:nrow(fr_df), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(fr_df$DPM_wu.r[i], fr_df$RPM_wu.r[i], fr_df$BIO_wu.r[i], fr_df$HUM_wu.r[i], fr_df$IOM_wu.r[i]),
    In   = fr_df$CinputFORWARD.r[i],
    Dr   = fr_df$DR[i],
    clay = fr_df$CLAY[i],
    effcts = data.frame(years, rep(unlist(xi_r[i, 1:12]), length.out = length(years))),
    solver = "euler"
  )
  fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_bau, file = "rothC_fr_bau.rds")
stopCluster(clus)



## Cinput bau min ----
clus <- makeCluster(16)
clusterExport(clus, list("fr_df", "xi_min", "years", "carbonTurnover"))

Sys.time()
rothC_fr_bau_min <- parLapply(clus, 1:nrow(fr_df), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(fr_df$DPM_wu.min[i], fr_df$RPM_wu.min[i], fr_df$BIO_wu.min[i], fr_df$HUM_wu.min[i], fr_df$IOM_wu.min[i]),
    In   = fr_df$CinputFORWARD.min[i],
    Dr   = fr_df$DR[i],
    clay = fr_df$CLAY_min[i],
    effcts = data.frame(years, rep(unlist(xi_min[i, 1:12]), length.out = length(years))),
    solver = "euler"
  )
  fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_bau_min, file = "rothC_fr_bau_min.rds")
stopCluster(clus)



## Cinput bau max ----
clus <- makeCluster(16)
clusterExport(clus, list("fr_df", "xi_max", "years", "carbonTurnover"))

Sys.time()
rothC_fr_bau_max <- parLapply(clus, 1:nrow(fr_df), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(fr_df$DPM_wu.max[i], fr_df$RPM_wu.max[i], fr_df$BIO_wu.max[i], fr_df$HUM_wu.max[i], fr_df$IOM_wu.max[i]),
    In   = fr_df$CinputFORWARD.max[i],
    Dr   = fr_df$DR[i],
    clay = fr_df$CLAY_max[i],
    effcts = data.frame(years, rep(unlist(xi_max[i, 1:12]), length.out = length(years))),
    solver = "euler"
  )
  fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_bau_max, file = "rothC_fr_bau_max.rds")
stopCluster(clus)



## Cinput low ----
clus <- makeCluster(16)
clusterExport(clus, list("fr_df", "xi_r", "years", "carbonTurnover"))

Sys.time()
rothC_fr_low <- parLapply(clus, 1:nrow(fr_df), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(fr_df$DPM_wu.r[i], fr_df$RPM_wu.r[i], fr_df$BIO_wu.r[i], fr_df$HUM_wu.r[i], fr_df$IOM_wu.r[i]),
    In   = fr_df$CinputFORWARD_low[i],
    Dr   = fr_df$DR[i],
    clay = fr_df$CLAY[i],
    effcts = data.frame(years, rep(unlist(xi_r[i, 1:12]), length.out = length(years))),
    solver = "euler"
  )
  fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_low, file = "rothC_fr_low.rds")
stopCluster(clus)



## Cinput med ----
clus <- makeCluster(16)
clusterExport(clus, list("fr_df", "xi_r", "years", "carbonTurnover"))

Sys.time()
rothC_fr_med <- parLapply(clus, 1:nrow(fr_df), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(fr_df$DPM_wu.r[i], fr_df$RPM_wu.r[i], fr_df$BIO_wu.r[i], fr_df$HUM_wu.r[i], fr_df$IOM_wu.r[i]),
    In   = fr_df$CinputFORWARD_med[i],
    Dr   = fr_df$DR[i],
    clay = fr_df$CLAY[i],
    effcts = data.frame(years, rep(unlist(xi_r[i, 1:12]), length.out = length(years))),
    solver = "euler"
  )
  fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_med, file = "rothC_fr_med.rds")
stopCluster(clus)



## Cinput high ----
clus <- makeCluster(16)
clusterExport(clus, list("fr_df", "xi_r", "years", "carbonTurnover"))

Sys.time()
rothC_fr_high <- parLapply(clus, 1:nrow(fr_df), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(fr_df$DPM_wu.r[i], fr_df$RPM_wu.r[i], fr_df$BIO_wu.r[i], fr_df$HUM_wu.r[i], fr_df$IOM_wu.r[i]),
    In   = fr_df$CinputFORWARD_high[i],
    Dr   = fr_df$DR[i],
    clay = fr_df$CLAY[i],
    effcts = data.frame(years, rep(unlist(xi_r[i, 1:12]), length.out = length(years))),
    solver = "euler"
  )
  fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_high, file = "rothC_fr_high.rds")
stopCluster(clus)



## Cinput med min ----
clus <- makeCluster(16)
clusterExport(clus, list("fr_df", "xi_min", "years", "carbonTurnover"))

Sys.time()
rothC_fr_medmin <- parLapply(clus, 1:nrow(fr_df), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(fr_df$DPM_wu.min[i], fr_df$RPM_wu.min[i], fr_df$BIO_wu.min[i], fr_df$HUM_wu.min[i], fr_df$IOM_wu.min[i]),
    In   = fr_df$CinputFORWARD_med_min[i],
    Dr   = fr_df$DR[i],
    clay = fr_df$CLAY_min[i],
    effcts = data.frame(years, rep(unlist(xi_min[i, 1:12]), length.out = length(years))),
    solver = "euler"
  )
  fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_medmin, file = "rothC_fr_medmin.rds")
stopCluster(clus)



## Cinput med max ----
clus <- makeCluster(16)
clusterExport(clus, list("fr_df", "xi_max", "years", "carbonTurnover"))

Sys.time()
rothC_fr_medmax <- parLapply(clus, 1:nrow(fr_df), function(i) {
  
  temp <- carbonTurnover(
    tt   = years,
    C0   = c(fr_df$DPM_wu.max[i], fr_df$RPM_wu.max[i], fr_df$BIO_wu.max[i], fr_df$HUM_wu.max[i], fr_df$IOM_wu.max[i]),
    In   = fr_df$CinputFORWARD_med_max[i],
    Dr   = fr_df$DR[i],
    clay = fr_df$CLAY_max[i],
    effcts = data.frame(years, rep(unlist(xi_max[i, 1:12]), length.out = length(years))),
    solver = "euler"
  )
  fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_medmax, file = "rothC_fr_medmax.rds")
stopCluster(clus)



# Load rothC outputs ----
fr_df <- readRDS("fr_df_Cinputs.rds")

vars <- c(
  bau    = "rothC_fr_bau.rds",
  baumin = "rothC_fr_bau_min.rds",
  baumax = "rothC_fr_bau_max.rds",
  low    = "rothC_fr_low.rds",
  med    = "rothC_fr_med.rds",
  high   = "rothC_fr_high.rds",
  medmin = "rothC_fr_medmin.rds",
  medmax = "rothC_fr_medmax.rds"
)
roth_l <- lapply(1:length(vars), function(x) {
  cat("reading", vars[x], "\n")
  temp <- cbind(
    source = names(vars)[x],
    cell   = fr_df$cell,
    as.data.frame(do.call("rbind", readRDS(vars[x])))
  )
  return(temp)
})
rc_fr_all     <- do.call("rbind", roth_l)
idx <- 4:8
rc_fr_all$f_t <- rowSums(rc_fr_all[idx])
names(rc_fr_all)[idx] <- c("DPM_fr", "RPM_fr", "BIO_fr", "HUM_fr", "IOM_fr")


# rc_fr_all <- reshape(rc_fr_all, direction = "wide",
#                      idvar   = "cell",
#                      timevar = "source",
#                      v.names = names(rc_fr_all[c(3:8)])
# )


library(data.table)

rc_fr_all <- dcast(
  as.data.table(rc_fr_all),
  cell + time ~ source,
  value.var = c("DPM_fr", "RPM_fr", "BIO_fr", "HUM_fr", "IOM_fr", "f_t"),
  sep = "."
)
rc_fr_all <- as.data.frame(rc_fr_all)

vars <- c("aoi", "x", "y", "cell", "SOC", "CLAY", "LU", "SOC_t0.r", 
          "SOC_t0.min", "SOC_t0.max",
          "CinputFORWARD.r", 
          "CinputFORWARD.min", "CinputFORWARD.max"
          )
all(fr_df$cell == rc_fr_all$cell)
rc_fr_all$cell <- NULL
rc_fr_final <- cbind(fr_df[vars], rc_fr_all)


saveRDS(rc_fr_final, file = "rothC_fr_final.rds")


