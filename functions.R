# TerraClimate ----
# function to average months
avgRS <- function(rs, var, interval) {
  
  idx_r  <- class(rs) == "RasterStack"
  int_nm <- names(interval)
  
  rs_l <- lapply(interval, function(x) {
    vars <- paste0(x, "_", var)
    idx  <- which(names(rs) %in% vars)
    cat("averaging", paste(names(rs[[idx]]), collapse = ", "), "\n", "\n")
    mean(rs[[idx]], na.rm = TRUE)
  })
  
  if (idx_r) {
    final_rs <- stack(rs_l)
  } else final_rs <- rast(rs_l)
  
  if (!is.null(int_nm)) names(final_rs) <- int_nm
  
  return(final_rs)
}


# function to sum layers
sumRS <- function(rs, var, interval) {
  
  idx_r  <- class(rs) == "RasterStack"
  int_nm <- names(interval)
  
  rs_l <- lapply(interval, function(x) {
    vars <- paste0(x, "_", var)
    idx  <- which(names(rs) %in% vars)
    cat("summing", paste(names(rs[[idx]]), collapse = ", "), "\n", "\n")
    sum(rs[[idx]], na.rm = TRUE)
  })
  
  if (idx_r) {
    final_rs <- stack(rs_l)
  } else final_rs <- rast(rs_l)
  
  if (!is.null(int_nm)) names(final_rs) <- int_nm
  
  return(final_rs)
}


# function to average years
avgRSyr <- function(rs, var, yrs) { 
  
  idx_r <- class(rs) == "RasterStack"
  
  rs_l <- lapply(yrs, function(x) {
    vars <- paste0(x, "_", var)
    idx  <- which(names(rs) %in% vars)
    cat("averaging", paste(names(rs[[idx]]), collapse = ", "), "\n", "\n")
    mean(rs[[idx]], na.rm = TRUE)
  })
  
  if (idx_r) {
    final_rs <- stack(rs_l)
  } else final_rs <- rast(rs_l)
  
  # names(final_rs) <- paste0("X", substr(x[[1]][1], 2, 5))
  return(final_rs)
}


fW <- function(pClay, PREC, PET, COV, s_thk = 30, pE = 1) {
  
  M     <- PREC - PET * pE
  
  n   <- ncol(M)  / 12
  idx <- rep(1:12, n)
  B   <- as.data.frame(lapply(COV[, idx], function(x) {
    ifelse(x < 0.8, 1, 1.8)
    }))
  
  Max.TSMD <- -(20 + 1.3 * pClay - 0.01 * (pClay^2)) * (s_thk/23) * (1/B)
  
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
  
  n <- ncol(Max.TSMD)
  for (i in 1:ncol(Acc.TSMD)) {
    
    Acc.TSMD[, i] <- ifelse(
      # Acc.TSMD[, i] > 0.444 * Max.TSMD[, i], 
      Acc.TSMD[, i] > 0.444 * Max.TSMD[, n], 
      1, 
      (0.2 + 0.8 * ((Max.TSMD[, n] - Acc.TSMD[, i])/(Max.TSMD[, n] - 0.444 * Max.TSMD[, n])))
      # (0.2 + 0.8 * ((Max.TSMD[, i] - Acc.TSMD[, i])/(Max.TSMD[, i] - 0.444 * Max.TSMD[, i])))
      )
    Acc.TSMD[, i] <- raster::clamp(Acc.TSMD[, i], lower = 0.2)
  }
  
  return(Acc.TSMD)
}



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
) {

  wu_l <- parLapply(clus, 1:nrow(su_df), function(i) {
    
    ncols  <- ncol(xi_df)
    nyrs    <- ncols / 12 
    df_idx <- data.frame(
      yr = rep(1:nyrs, each = 12),
      cn = 1:ncols
    )
    idx <- df_idx[df_idx$yr == 1, ]$cn
    
    temp <- carbonTurnover(
      tt   = time,
      C0   = as.numeric(su_df[i, C0_vars]),
      In   = C_m[i, 1],
      Dr   = su_df[i, DR_var],
      clay = su_df[i, pClay_var],
      effcts = data.frame(
        time, 
        xi = unlist(xi_df[i, idx])
      ),
      solver = model
    )
    fp <- list(tail(temp, 1))
    
    
    for (j in 2:nSim) {

      idx2 <- df_idx[df_idx$yr == j, ]$cn

      temp <- carbonTurnover(
        tt   = time,
        C0   = fp[[1]][2:6],
        In   = C_m[i, j],
        Dr   = su_df[i, DR_var],
        clay = su_df[i, pClay_var],
        effcts = data.frame(
          time,
          xi = unlist(xi_df[i, idx2])
          ),
        solver = model
      )
      fp[[1]] <- tail(temp, 1)
    }
    fp <- unlist(fp)
    fp[1] <- nSim
    
    return(fp)
  })
}


rothC_wu_nn <- function(
  time = seq(1 / 12, 1, by = 1 / 12), 
  nSim = 19, 
  su_df, pClay_var, C0_vars, DR_var, 
  C_m, 
  xi_df
) {
  
  
  wu_l <- parLapply(clus, 1:nrow(su_df), function(i) {
    
    ncols  <- ncol(xi_df)
    nyrs    <- ncols / 12 
    df_idx <- data.frame(
      yr = rep(1:nyrs, each = 12),
      cn = 1:ncols
    )
    idx <- df_idx[df_idx$yr == 1, ]$cn
    
    temp <- RothCModel(
      t    = time,
      C0   = as.numeric(su_df[i, C0_vars]),
      In   = C_m[i, 1],
      DR   = su_df[i, DR_var],
      clay = su_df[i, pClay_var],
      xi   = data.frame(
        time, 
        xi = unlist(xi_df[i, idx])
      ),
      pass = TRUE
    )
    fp <- list(tail(getC(temp), 1))
    
    for (j in 2:nSim) {
      
      idx2 <- df_idx[df_idx$yr == j, ]$cn
      
      temp <- RothCModel(
        t    = time,
        C0   = as.numeric(fp[[1]][1:5]),
        In   = C_m[i, j],
        DR   = su_df[i, DR_var],
        clay = su_df[i, pClay_var],
        xi   = data.frame(
          years, 
          xi = unlist(xi_df[i, idx2])
        ),
        pass = TRUE
      )
      fp[[1]] <- tail(getC(temp), 1)
    }
    
    fp <- unlist(fp)
    
    return(fp)
  })
}



fget_equilibrium_fractions.RothC_input=function(xi=1,C.tot,clay, fractI)
{   
  rmf=xi
  IOM= fIOM.Falloon.RothC(SOC = C.tot)
  C.active=C.tot-IOM 
  
  ########################################################################
  #The analytical solution of RothC
  ########################################################################
  
  ########################################################################
  # Parameter
  ########################################################################
  fract.rooted.to.bio = 0.46
  fract.rooted.to.hum = 0.54
  ks = c(k.DPM = 10, k.RPM = 0.3, k.BIO = 0.66, k.HUM = 0.02, 
         k.IOM = 0)
  ks=as.numeric(ks)
  k.dpm=ks[1]
  k.rpm=ks[2]
  k.bio=ks[3]
  k.hum=ks[4]
  ########################################################################
  # the carbon use efficiency
  ########################################################################
  cue=  1/(1+ 1.67 * (1.85 + 1.6 * exp(-0.0786 * clay)))
  
  ########################################################################
  # All the coefficients alpha.1 und alpha.2
  ########################################################################
  alpha.1=cue*fract.rooted.to.bio
  alpha.2=cue*fract.rooted.to.hum
  
  ########################################################################
  # All the coefficients a.1.1, a.1.2, a.2.1, a2.2
  ########################################################################
  a.1.1=k.bio*rmf*(alpha.1-1)
  a.1.2=alpha.1*k.hum*rmf
  a.2.1=alpha.2*k.bio*rmf
  a.2.2=k.hum*rmf*(alpha.2-1)
  
  #########################################################################
  #########################################################################
  # The Eigenvalues lambda 1 and lambda 2
  #########################################################################
  lambda.1= (a.1.1+a.2.2)/2-sqrt(((a.1.1+a.2.2)/2)*((a.1.1+a.2.2)/2)+a.1.2*a.2.1-a.1.1*a.2.2)
  lambda.2= (a.1.1+a.2.2)/2+sqrt(((a.1.1+a.2.2)/2)*((a.1.1+a.2.2)/2)+a.1.2*a.2.1-a.1.1*a.2.2)
  #########################################################################
  # The c.0.1; c.0.2; c.0.3 values
  #########################################################################
  c.0.1= (alpha.2 * a.1.2 - alpha.1 * a.2.2)/(a.1.1*a.2.2-a.1.2*a.2.1)
  c.0.2= (alpha.2 * a.1.2 - alpha.1 * a.2.2)/(a.1.1*a.2.2-a.1.2*a.2.1)
  c.0.3= (a.1.2)/(a.1.1*a.2.2-a.1.2*a.2.1)
  
  ######################################################################################################
  # BIO pool quantification
  ######################################################################################################
  u.bio.dpm=(c.0.2) #65
  u.bio.rpm=(c.0.1) #66
  u.bio.hum=(c.0.3) #67
  
  
  ######################################################################################################
  # HUM pool quantification ( is all C.78)
  ######################################################################################################
  u.hum.dpm= 1/a.1.2*((-c.0.2*a.1.1-alpha.1))
  u.hum.rpm= 1/a.1.2*(-c.0.2*a.1.1-alpha.1)
  u.hum.hum= 1/a.1.2*(-c.0.3*a.1.1)
  
  
  ######################################################################################################
  # DPM C ( is all C.79)
  ######################################################################################################
  u.dpm.dpm=1/k.dpm/rmf 
  
  #C.dpm=i.dpm * u.dpm.dpm + C0 * s.dpm
  
  
  ######################################################################################################
  # RPM C ( is all C.80)
  ######################################################################################################
  u.rpm.rpm=1/k.rpm/rmf
  
  #C.rpm=i.rpm * u.rpm.rpm + C0 *s.rpm
  
  ######################################################################################################
  # Total C ( is all C.78)
  ######################################################################################################
  u.dpm=u.dpm.dpm+u.bio.dpm+u.hum.dpm
  u.rpm=u.rpm.rpm+u.bio.rpm+u.hum.rpm
  u.hum=u.bio.hum+u.hum.hum
  
  denominator= fractI[1]*u.dpm+fractI[2]*u.rpm+fractI[3]*u.hum
  
  fract.dpm= fractI[1]*u.dpm.dpm/denominator
  fract.rpm= fractI[2]*u.rpm.rpm/denominator
  fract.bio= (fractI[1]*u.bio.dpm+fractI[2]*u.bio.rpm+fractI[3]*u.bio.hum)/denominator
  fract.hum= (fractI[1]*u.hum.dpm+fractI[2]*u.hum.rpm+fractI[3]*u.hum.hum)/denominator   
  
  fract.all=c(fract.dpm,fract.rpm,fract.bio,fract.hum)
  
  ###################################################
  # so unfortunately we have the IOM
  ###################################################
  fract.all_stock=(fract.all*C.active)
  fract.all=fract.all_stock/C.tot
  fract.all=append(fract.all,IOM/C.tot)
  pools=fract.all*C.tot
  Cin=(C.tot-pools[5])/denominator
  # list(pools,Cin)
  names(pools) <- c("fract.dpm", "fract.rpm", "fract.bio", "fract.hum", "fract.iom")
  return(c(pools, Cin = Cin))
}



fIOM.Falloon.RothC <- function(SOC, par1 = -1.31, par2 = 1.139) {
  # IOM=10^(par1+par2*log10(c))
  IOM = 0.049 * SOC^1.139 
  IOM
}
