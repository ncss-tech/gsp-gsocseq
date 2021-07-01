library(SoilR)
library(raster)
library(rgdal)
library(soilassessment)
library(sf)


# Set working directory
setwd("D:/geodata/project_data/gsp-gsocseq")


# load forward stack
fr_df <- as.data.frame(readRDS(file = "fr_sdf.RDS"))
fr_df <- fr_df[complete.cases(fr_df), ]
fr_sf <- st_as_sf(
    fr_df,
    coords = c("x", "y"),
    crs    = 4326
    ) %>%
    st_transform(crs = 5070)
fr_df <- cbind(id = 1:nrow(fr_sf), round(st_coordinates(fr_sf)), st_drop_geometry(fr_sf))
fr_df$xy <- paste(fr_df$X, fr_df$Y, sep = "_")


# load warm up
wu_df <- readRDS(file = "rothC_r_wu_final.rds")
wu_df$xy <- paste(wu_df$X, wu_df$Y, sep = "_")


# combine
fr_df <- merge(fr_df, wu_df, by = "xy", all.x = TRUE)
fr_df <- fr_df[complete.cases(fr_df), ]


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
fr_df <- within(fr_df, {
    CinputFORWARD_low     <- ifelse(LU %in% c(2, 12, 3, 5, 6, 8, 13), CinputFORWARD.r    * 1.05,        CinputFORWARD.r)
    CinputFORWARD_med     <- ifelse(LU %in% c(2, 12, 3, 5, 6, 8, 13), CinputFORWARD.r   * 1.10,        CinputFORWARD.r)
    CinputFORWARD_high    <- ifelse(LU %in% c(2, 12, 3, 5, 6, 8, 13), CinputFORWARD.r   * 1.20,        CinputFORWARD.r)
    
    CinputFORWARD_med_min <- ifelse(LU %in% c(2, 12, 3, 5, 6, 8, 13), CinputFORWARD.min * 1.05 - 0.15, CinputFORWARD.min) 
    CinputFORWARD_med_max <- ifelse(LU %in% c(2, 12, 3, 5, 6, 8, 13), CinputFORWARD.max * 1.20 + 0.15, CinputFORWARD.max)
})

saveRDS(fr_df, file = "fr_df.rds")


# Extract variables
pClay_r <- fr_df$CLAY
TEMP    <- fr_df[grepl("TEMP_", names(fr_df))]
PREC    <- fr_df[grepl("PREC_", names(fr_df))]
PET     <- fr_df[grepl("PET_",  names(fr_df))]
COV     <- fr_df[grepl("COV_",  names(fr_df))]
LU      <- fr_df$LU


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


xi_r   <- cbind(id = 1:nrow(xi_r),   xi_r)
xi_min <- cbind(id = 1:nrow(xi_min), xi_min)
xi_max <- cbind(id = 1:nrow(xi_max), xi_max)


saveRDS(xi_r,   "fr_effcts_r.rds")
saveRDS(xi_min, "fr_effcts_min.rds")
saveRDS(xi_max, "fr_effcts_max.rds")


# RUN THE MODEL from soilassessment ----
# Roth C soilassesment in parallel
years <- seq(1 / 12, 20, by = 1 / 12)

fr_df <- readRDS("fr_df.rds")
fr_df <- fr_df[!grepl("TEMP_|PREC_|PET_|COV_", names(fr_df))]
fr_df$CLAY_min <- fr_df$CLAY * 0.9
fr_df$CLAY_max <- fr_df$CLAY * 1.1

xi_r   <- readRDS("fr_effcts_r.rds")
xi_min <- readRDS("fr_effcts_min.rds")
xi_max <- readRDS("fr_effcts_max.rds")


library(parallel)


# Cinput bau ----
clus <- makeCluster(10)
clusterExport(clus, list("fr_df", "xi_r", "years", "carbonTurnover"))

Sys.time()
rothC_bau <- parLapply(clus, 1:nrow(fr_df), function(i) {
    
    temp <- carbonTurnover(
        tt   = years,
        C0   = c(fr_df$DPM_wu.r[i], fr_df$RPM_wu.r[i], fr_df$BIO_wu.r[i], fr_df$HUM_wu.r[i], fr_df$IOM_wu.r[i]),
        In   = fr_df$CinputFORWARD.r[i],
        Dr   = fr_df$DR[i],
        clay = fr_df$CLAY[i],
        effcts = data.frame(years, rep(unlist(xi_r[i, 2:13]), length.out = length(years))),
        solver = "euler"
    )
    fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_r, file = "rothC_fr_bau.rds")
stopCluster(clus)



# Cinput bau min ----
clus <- makeCluster(10)
clusterExport(clus, list("fr_df", "xi_min", "years", "carbonTurnover"))

Sys.time()
rothC_fr_bau_min <- parLapply(clus, 1:nrow(fr_df), function(i) {
    
    temp <- carbonTurnover(
        tt   = years,
        C0   = c(fr_df$DPM_wu.min[i], fr_df$RPM_wu.min[i], fr_df$BIO_wu.min[i], fr_df$HUM_wu.min[i], fr_df$IOM_wu.min[i]),
        In   = fr_df$CinputFORWARD.min[i],
        Dr   = fr_df$DR[i],
        clay = fr_df$CLAY_min[i],
        effcts = data.frame(years, rep(unlist(xi_min[i, 2:13]), length.out = length(years))),
        solver = "euler"
    )
    fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_bau_min, file = "rothC_fr_bau_min.rds")
stopCluster(clus)



# Cinput bau max ----
clus <- makeCluster(10)
clusterExport(clus, list("fr_df", "xi_max", "years", "carbonTurnover"))

Sys.time()
rothC_fr_bau_max <- parLapply(clus, 1:nrow(fr_df), function(i) {
    
    temp <- carbonTurnover(
        tt   = years,
        C0   = c(fr_df$DPM_wu.max[i], fr_df$RPM_wu.max[i], fr_df$BIO_wu.max[i], fr_df$HUM_wu.max[i], fr_df$IOM_wu.max[i]),
        In   = fr_df$CinputFORWARD.max[i],
        Dr   = fr_df$DR[i],
        clay = fr_df$CLAY_max[i],
        effcts = data.frame(years, rep(unlist(xi_max[i, 2:13]), length.out = length(years))),
        solver = "euler"
    )
    fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_bau_max, file = "rothC_fr_bau_max.rds")
stopCluster(clus)



# Cinput low ----
clus <- makeCluster(10)
clusterExport(clus, list("fr_df", "xi_r", "years", "carbonTurnover"))

Sys.time()
rothC_fr_low <- parLapply(clus, 1:nrow(fr_df), function(i) {
    
    temp <- carbonTurnover(
        tt   = years,
        C0   = c(fr_df$DPM_wu.r[i], fr_df$RPM_wu.r[i], fr_df$BIO_wu.r[i], fr_df$HUM_wu.r[i], fr_df$IOM_wu.r[i]),
        In   = fr_df$CinputFORWARD_low[i],
        Dr   = fr_df$DR[i],
        clay = fr_df$CLAY[i],
        effcts = data.frame(years, rep(unlist(xi_r[i, 2:13]), length.out = length(years))),
        solver = "euler"
    )
    fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_low, file = "rothC_fr_low.rds")
stopCluster(clus)



# Cinput med ----
clus <- makeCluster(10)
clusterExport(clus, list("fr_df", "xi_r", "years", "carbonTurnover"))

Sys.time()
rothC_fr_med <- parLapply(clus, 1:nrow(fr_df), function(i) {
    
    temp <- carbonTurnover(
        tt   = years,
        C0   = c(fr_df$DPM_wu.r[i], fr_df$RPM_wu.r[i], fr_df$BIO_wu.r[i], fr_df$HUM_wu.r[i], fr_df$IOM_wu.r[i]),
        In   = fr_df$CinputFORWARD_med[i],
        Dr   = fr_df$DR[i],
        clay = fr_df$CLAY[i],
        effcts = data.frame(years, rep(unlist(xi_r[i, 2:13]), length.out = length(years))),
        solver = "euler"
    )
    fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_med, file = "rothC_fr_med.rds")
stopCluster(clus)



# Cinput high ----
clus <- makeCluster(10)
clusterExport(clus, list("fr_df", "xi_r", "years", "carbonTurnover"))

Sys.time()
rothC_fr_high <- parLapply(clus, 1:nrow(fr_df), function(i) {
    
    temp <- carbonTurnover(
        tt   = years,
        C0   = c(fr_df$DPM_wu.r[i], fr_df$RPM_wu.r[i], fr_df$BIO_wu.r[i], fr_df$HUM_wu.r[i], fr_df$IOM_wu.r[i]),
        In   = fr_df$CinputFORWARD_high[i],
        Dr   = fr_df$DR[i],
        clay = fr_df$CLAY[i],
        effcts = data.frame(years, rep(unlist(xi_r[i, 2:13]), length.out = length(years))),
        solver = "euler"
    )
    fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_high, file = "rothC_fr_high.rds")
stopCluster(clus)



# Cinput med min ----
clus <- makeCluster(10)
clusterExport(clus, list("fr_df", "xi_min", "years", "carbonTurnover"))

Sys.time()
rothC_fr_medmin <- parLapply(clus, 1:nrow(fr_df), function(i) {
    
    temp <- carbonTurnover(
        tt   = years,
        C0   = c(fr_df$DPM_wu.min[i], fr_df$RPM_wu.min[i], fr_df$BIO_wu.min[i], fr_df$HUM_wu.min[i], fr_df$IOM_wu.min[i]),
        In   = fr_df$CinputFORWARD_med_min[i],
        Dr   = fr_df$DR[i],
        clay = fr_df$CLAY_min[i],
        effcts = data.frame(years, rep(unlist(xi_min[i, 2:13]), length.out = length(years))),
        solver = "euler"
    )
    fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_medmin, file = "rothC_fr_medmin.rds")
stopCluster(clus)



# Cinput med max ----
clus <- makeCluster(10)
clusterExport(clus, list("fr_df", "xi_max", "years", "carbonTurnover"))

Sys.time()
rothC_fr_medmax <- parLapply(clus, 1:nrow(fr_df), function(i) {
    
    temp <- carbonTurnover(
        tt   = years,
        C0   = c(fr_df$DPM_wu.max[i], fr_df$RPM_wu.max[i], fr_df$BIO_wu.max[i], fr_df$HUM_wu.max[i], fr_df$IOM_wu.max[i]),
        In   = fr_df$CinputFORWARD_med_max[i],
        Dr   = fr_df$DR[i],
        clay = fr_df$CLAY_max[i],
        effcts = data.frame(years, rep(unlist(xi_max[i, 2:13]), length.out = length(years))),
        solver = "euler"
    )
    fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_medmax, file = "rothC_fr_medmax.rds")
stopCluster(clus)



# load rothC outputs ----
fr_df <- readRDS("fr_df.rds")

rc_fr_bau    <- as.data.frame(do.call("rbind", readRDS(file = "rothC_fr_bau.rds")))
rc_fr_baumin <- as.data.frame(do.call("rbind", readRDS(file = "rothC_fr_bau_min.rds")))
rc_fr_baumax <- as.data.frame(do.call("rbind", readRDS(file = "rothC_fr_bau_max.rds")))
rc_fr_low    <- as.data.frame(do.call("rbind", readRDS(file = "rothC_fr_low.rds")))
rc_fr_med    <- as.data.frame(do.call("rbind", readRDS(file = "rothC_fr_med.rds")))
rc_fr_high   <- as.data.frame(do.call("rbind", readRDS(file = "rothC_fr_high.rds")))
rc_fr_medmin <- as.data.frame(do.call("rbind", readRDS(file = "rothC_fr_medmin.rds")))
rc_fr_medmax <- as.data.frame(do.call("rbind", readRDS(file = "rothC_fr_medmax.rds")))


ids <- 1:nrow(rc_fr_bau)
rc_fr_all <- rbind(
    cbind(source = "bau",    id = ids, rc_fr_bau),
    cbind(source = "baumin", id = ids, rc_fr_baumin),
    cbind(source = "baumax", id = ids, rc_fr_baumax),
    cbind(source = "low",    id = ids, rc_fr_low),
    cbind(source = "med",    id = ids, rc_fr_med),
    cbind(source = "high",   id = ids, rc_fr_high),
    cbind(source = "medmin", id = ids, rc_fr_medmin),
    cbind(source = "medmax", id = ids, rc_fr_medmax)
)
rc_fr_all$f_t <- rowSums(rc_fr_all[4:8])
names(rc_fr_all)[4:8] <- c("DPM_fr", "RPM_fr", "BIO_fr", "HUM_fr", "IOM_fr")

rc_fr_all <- reshape(rc_fr_all, direction = "wide",
                     idvar = c("id"),
                     timevar = "source", 
                     v.names = names(rc_fr_all[c(4:9)])
)

vars <- c("x", "y", "SOC", "CLAY", "LU", "SOC_t0.r", "SOC_t0.min", "SOC_t0.max", "CinputFORWARD.r", "CinputFORWARD.min", "CinputFORWARD.max")
names(fr_df)[3:6] <- vars[1:4]

rc_fr_all2 <- cbind(fr_df[vars], rc_fr_all)



# Uncertainties
rc_fr_all <- within(rc_fr_all, {
    unc_soc = (f_t.baumax - f_t.baumin) / (2 * f_t.bau)  * 100
    unc_t0  = (SOC_t0.max - SOC_t0.min) / (2 * SOC_t0.r) * 100
    unc_ssm = (f_t.medmax - f_t.medmin) / (2 * f_t.med)  * 100
})



# Convert to points
names(rc_fr_all) <- gsub("\\.", "_", names(rc_fr_all))

rc_fr_final <- st_as_sf(
    rc_fr_all,
    coords = c("x", "y"),
    crs    = 5070
    ) %>%
    st_transform(4326)

saveRDS(rc_fr_final, file = "rothC_fr_final.rds")
write_sf(rc_fr_final, dsn = "rothC_fr_final.gpkg", driver = "GPKG", overwrite = TRUE) 



