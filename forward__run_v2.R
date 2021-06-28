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


# apply increase
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
fr_df$Clay_min <- fr_df$CLAY * 0.9
fr_df$Clay_max <- fr_df$CLAY * 1.1

xi_r   <- readRDS("fr_effcts_r.rds")
xi_min <- readRDS("fr_effcts_min.rds")
xi_max <- readRDS("fr_effcts_max.rds")


library(parallel)


# C input equilibrium. (Ceq) ----
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



# bau min ----
clus <- makeCluster(10)
clusterExport(clus, list("fr_df", "xi_min", "years", "carbonTurnover"))

Sys.time()
rothC_fr_bau_min <- parLapply(clus, 1:nrow(fr_df), function(i) {
    
    temp <- carbonTurnover(
        tt   = years,
        C0   = c(fr_df$DPM_wu.min[i], fr_df$RPM_wu.min[i], fr_df$BIO_wu.min[i], fr_df$HUM_wu.min[i], fr_df$IOM_wu.min[i]),
        In   = fr_df$CinputFORWARD.min[i],
        Dr   = fr_df$DR[i],
        clay = fr_df$Clay_min[i],
        effcts = data.frame(years, rep(unlist(xi_min[i, 2:13]), length.out = length(years))),
        solver = "euler"
    )
    fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_bau_min, file = "rothC_fr_bau_min.rds")
stopCluster(clus)


# bau min ----
clus <- makeCluster(10)
clusterExport(clus, list("fr_df", "xi_max", "years", "carbonTurnover"))

Sys.time()
rothC_fr_bau_max <- parLapply(clus, 1:nrow(fr_df), function(i) {
    
    temp <- carbonTurnover(
        tt   = years,
        C0   = c(fr_df$DPM_wu.max[i], fr_df$RPM_wu.max[i], fr_df$BIO_wu.max[i], fr_df$HUM_wu.max[i], fr_df$IOM_wu.max[i]),
        In   = fr_df$CinputFORWARD.max[i],
        Dr   = fr_df$DR[i],
        clay = fr_df$Clay_max[i],
        effcts = data.frame(years, rep(unlist(xi_max[i, 2:13]), length.out = length(years))),
        solver = "euler"
    )
    fp <- tail(temp, 1)
})
Sys.time()
# saveRDS(rothC_fr_bau_max, file = "rothC_fr_bau_max.rds")
stopCluster(clus)


f_bau_min<-Roth_C(Cinputs=Cinputs_min,years=years,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
f_bau_t_min<-f_bau_min[1]+f_bau_min[2]+f_bau_min[3]+f_bau_min[4]+f_bau_min[5]

#Unc BAU maximum

f_bau_max<-Roth_C(Cinputs=Cinputs_max,years=years,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
f_bau_t_max<-f_bau_max[1]+f_bau_max[2]+f_bau_max[3]+f_bau_max[4]+f_bau_max[5]


