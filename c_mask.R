rm(list = ls())

library(raster)


setwd("D:/geodata/masks")
ev <- c()
rs <- list.files("D:\\geodata\\masks", pattern=".tif$", full.names = TRUE)

names(rs) <- 1:length(rs)
names(rs)

rs

# remove conus
rs <- rs[-9]

names(rs) <- 1:length(rs)
rs
mast <- raster(rs[9])
fips <- "USVI850"

plot(mast)



setwd("D:\\geodata\\project_data\\gsp-gsocseqp")

dir.create("RESULTS")
# fips <- "__GU316"



socrs <- list.files(".", pattern = "^PRVI")

for(soc in socrs){
  r <- raster(soc)
  rn <- gsub("PRVI", fips, soc)
  # # print(nn)
  # writeRaster(msk, filename =  paste0("./RESULTS/", rn), progress = "text", overwrite = TRUE)
  
  
  gdalUtilities::gdalwarp(
    srcfile = soc,
    dstfile = paste0("./RESULTS/", rn),
    te      = c(bbox(mast)),
    tr      = res(r),
    r       = "bilinear",
    ot      = "Float64",
    co      = "COMPRESS=DEFLATE",
    overwrite = TRUE
  )
  
}
traceback()
  