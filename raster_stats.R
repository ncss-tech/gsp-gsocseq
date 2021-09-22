library(raster)
library(xlsx)

setwd("D:/geodata/project_data/gsp-gsocseqp")
# ml <- list("file", "min", "max")

rasters <- list.files(".", pattern = "^PRVI")

ml <- list()
for (rast in rasters){
  r <- raster(rast)
  name <- basename(r@file@name)
  # mn <- r@data@min
  # mx <- r@data@max
  vals <- raster::values(r)
  vals <- na.omit(vals)
  vals <- vals[!vals %in% c(0, -32768)]
  mn <- min(vals)
  med <- median(vals)
  mea <- mean(vals)
  mx <- max(vals)
  # print(paste0(name, as.character(mn), as.character(mx), sep = "\t"))
  data <-  c(name,mn,med, mea,mx)
  ml <- append(ml, list(data)) 
  # (name, "\t", as.character(mn), "\t", as.character(mx), "\n")
}

df <- data.frame(matrix(unlist(ml), nrow=length(ml), byrow = TRUE))
names(df) <- c("file", "min", "median", "mean", "max")
write.xlsx(df, file = 'raster_stats.xlsx', sheetName = 'Sheet1', col.names = TRUE, row.names = FALSE)

saveRDS(df, "D:\\geodata\\new_raster_stats.rds")

inDF <- readRDS("D:\\geodata\\new_raster_stats.rds")
oldDF <- readRDS("D:\\geodata\\old_raster_stats.rds")
