library(SoilR)
library(raster)
library(rgdal)
library(soilassessment)


# Set working directory
setwd("D:/geodata/project_data/gsp-gsocseq")


# Stack_Set_1 is a stack that contains the spatial variables 
wu_df <- readRDS(file = "wu_pts.rds")
wu_df  <- wu_df[complete.cases(wu_df), ]
wu_df  <- cbind(ID = 1:nrow(wu_df), wu_df)


# Extract the layers from the Vector
SOC_r    <- wu_df$CONUS_GSOCmap1.5.0
SOC_min  <- wu_df$SOC * 0.8
SOC_max  <- wu_df$SOC * 1.2

pClay_r   <- wu_df$CONUS_gnatsgo_fy20_1km_clay_wt
pClay_min <- wu_df$CLAY * 0.9
pClay_max <- wu_df$CLAY * 1.1

DR   <- wu_df$CONUS_glc_shv10_DOM_DR
LU   <- wu_df$CONUS_glc_shv10_DOM
TEMP <- wu_df[grepl("_tmmx$", names(wu_df))]
PREC <- wu_df[grepl("_pr$", names(wu_df))]
PET  <- wu_df[grepl("_pet$",  names(wu_df))]
COV  <- wu_df[grepl("^CON_",  names(wu_df))]
COV  <- do.call("cbind", rep(list(COV), 19))


# Define Years for Cinputs calculations ----
year <- seq(1 / 12, 1, by = 1 / 12)
nSim <- ncol(TEMP) / 12

# Paddy fields coefficent  ----
# fPR = 0.4 if the target point is class = 13 , else fPR=1
# From Shirato and Yukozawa 2004
fPR <- (LU == 13) * 0.4 + (LU != 13) * 1


# Temperature effects per month ----
fT <- as.data.frame(lapply(TEMP, function(x) fT.RothC(x)))
