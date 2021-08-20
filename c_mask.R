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
mast <- raster(rs[4])
fips <- "USHI840"

plot(mast)



setwd("D:\\geodata\\fin_project_data4\\fin_project_data4\\gsp-gsocseqwp")

dir.create("RESULTS")
# fips <- "__GU316"



socrs <- list.files(".", pattern = "^HIAS")

for(soc in socrs){
  r <- raster(soc)
  rn <- gsub("HIAS", fips, soc)
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
  