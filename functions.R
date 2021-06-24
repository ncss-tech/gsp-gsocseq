
NPPmodel <- function (rain, temperature, model = "miami") 
{
  if (model == "miami") {
    npprain = 3000 * (1 - exp(-0.000664 * rain))
    npptemp = 3000/(1 + exp(1.315 - 0.119 * temperature))
    npp = pmin(npprain, npptemp)
  }
  else if (model == "schuur") {
    npprain = (0.005215 * (rain)^1.12363)/exp(0.000459532 * 
                                                rain)
    npptemp = 17.6243/(1 + exp(1.3496 - 0.071514 * temperature))
    npp = pmin(npprain, npptemp) * 200
  }
  else if (model == "NCEAS") {
    npprain = (0.551 * (rain)^1.055)/exp(0.000306 * rain)
    npptemp = 2540/(1 + exp(1.584 - 0.0622 * temperature))
    npp = pmin(npprain, npptemp) * 2
  }
  return(npp)
}



rothC_wu <- function(
  
  time = seq(1 / 12, 1, by = 1 / 12), 
  nSim = 19, 
  su_df, pClay_var, C0_vars, DR_var, 
  C_m, 
  xi_df, 
  model = "euler" 
  # n_clus = 7
  ) {
  
  # clus <- makeCluster(n_clus)
  # 
  # clusterExport(clus, list("su_df", "C_m", "xi_df", "time", "nSim", "carbonTurnover"))
  
  wu_l <- parLapply(clus, 1:nrow(su_df), function(i) {
    
    temp <- carbonTurnover(
      tt   = time,
      C0   = as.numeric(su_df[i, C0_vars]),
      In   = C_m[i, 1],
      Dr   = su_df[i, DR_var],
      clay = su_df[i, pClay_var],
      effcts = data.frame(
        time, 
        rep(unlist(xi_df[i, 2:ncol(xi_df)]), length.out = length(time))
        ),
      solver = model
      )
    fp <- list(tail(temp, 1))
    
    for (j in 2:nSim) {
      temp <- carbonTurnover(
        tt   = time,
        C0   = c(fp[[1]][2:6]),
        In   = C_m[i, j],
        Dr   = su_df[i, DR_var],
        clay = su_df[i, pClay_var],
        effcts = data.frame(
          time, rep(unlist(xi_df[i, 2:ncol(xi_df)]), length.out = length(time))
        ),
        solver = model
        )
      fp[[1]] <- tail(temp, 1)
    }
    fp <- unlist(fp)
    fp[1] <- nSim
    
    return(fp)
    })
  
  # stopCluster(clus)
}


rothC_wu_nn <- function(
  time = seq(1 / 12, 1, by = 1 / 12), 
  nSim = 19, 
  su_df, pClay_var, C0_vars, DR_var, 
  C_m, 
  xi_df
  ) {
    
  
  wu_l <- parLapply(clus, 1:nrow(su_df), function(i) {
    
    temp <- RothCModel(
      t    = time,
      C0   = as.numeric(su_df[i, C0_vars]),
      In   = C_m[i, 1],
      DR   = su_df[i, DR_var],
      clay = su_df[i, pClay_var],
      xi   = data.frame(
        time, 
        rep(unlist(xi_df[i, 2:ncol(xi_df)]), length.out = length(time))
        ),
      pass = TRUE
      )
    fp <- list(tail(getC(temp), 1))
    
    for (j in 2:19) {
      temp <- RothCModel(
        t    = time,
        C0   = as.numeric(fp[[1]][1:5]),
        In   = C_m[i, j],
        DR   = su_df[i, DR_var],
        clay = su_df[i, pClay_var],
        xi   = data.frame(
          years, rep(unlist(xi_df[i, 2:ncol(xi_df)]), length.out = length(time))
          ),
        pass = TRUE
        )
      fp[[1]] <- tail(getC(temp), 1)
    }
    
    fp <- unlist(fp)
    
    return(fp)
  })
}
