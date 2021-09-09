# Spin up (equilibrium run)

test1 <- readRDS("su_results_v3_er.rds")
test2 <- readOGR(dsn = "SPIN_UP_County_AOI_er.shp")

test1 <- test1[test1$cell %in% test2$cell, ]
test2 <- test2[test2$cell %in% test1$cell, ]

test1 <- test1[order(test1$cell), ]
test2 <- test2[order(test2$cell), ]

all(test1$cell == test2$cell)

test3 <- data.frame(
  cell = test1$cell,
  Cin1 = round(test1$Cin.r, 5), 
  Cin2 = round(test2$Cnpt_EQ, 5)
)
test3$dif <- test3$Cin1 - test3$Cin2

all(test3$Cin1 == test3$Cin2)


# Spin up

test1 <- readRDS("su_results_v3_analytical.rds")
test2 <- readOGR(dsn = "SPIN_UP_AOI.shp")

test1 <- test1[test1$cell %in% test2$cell, ]
test2 <- test2[test2$cell %in% test1$cell, ]

test1 <- test1[order(test1$cell), ]
test2 <- test2[order(test2$cell), ]

all(test1$cell == test2$cell)

test3 <- data.frame(
  cell = test1$cell,
  Cin1 = round(test1$Cin.r, 5), 
  Cin2 = round(test2$Cnpt_EQ, 5)
  )
test3$dif <- test3$Cin1 - test3$Cin2

all(test3$Cin1 == test3$Cin2)


# Warm up
test1 <- readRDS("rothC_r_wu_final.rds")
test2 <- readOGR(dsn = "WARM_UP_County_AOI.shp")

test1 <- test1[test1$cell %in% test2$cell, ]
test2 <- test2[test2$cell %in% test1$cell, ]

test1 <- test1[order(test1$cell), ]
test2 <- test2[order(test2$cell), ]

all(test1$cell == test2$cell)

test3 <- data.frame(
  cell = test1$cell,
  SOC1 = round(test1$SOC_t0.r, 5), 
  SOC2 = round(test2$SOC_t0, 5)
)
test3$dif <- test3$SOC1 - test3$SOC2

all(test3$SOC1 == test3$SOC2)


# Forward run
test1 <- readRDS("rothC_fr_final.rds")
test2 <- readOGR(dsn = "FOWARD_County_AOI.shp")

test1 <- test1[test1$cell %in% test2$cell, ]
test2 <- test2[test2$cell %in% test1$cell, ]

test1 <- test1[order(test1$cell), ]
test2 <- test2[order(test2$cell), ]

all(test1$cell == test2$cell)

test3 <- data.frame(
  cell = test1$cell,
  SOC1 = round(test1$f_t, 5), 
  SOC2 = round(test2$SOC_BAU_20, 5)
)
test3$dif <- test3$SOC1 - test3$SOC2

all(test3$SOC1 == test3$SOC2)



