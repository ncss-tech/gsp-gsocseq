rm(list=ls()) 

library(SoilR)
library(raster)
library(rgdal)
library(soilassessment)
library(sf)
library(terra)


# Set working directory 

# Set working directory 

setwd("D:/geodata/project_data/gsp-gsocseq/CONUS")
su_sf <- readRDS(file = "su_sf.RDS")
su_sf <- su_sf[order(su_sf$cell), ]
Vector <- as(su_sf, "Spatial")[1]


# Vector<-readOGR("INPUTS/TARGET_POINTS/Target_Points_sub.shp")


# Stack_Set_1 is a stack that contains the spatial variables 

Stack_Set_1<- stack("Stack_Set_SPIN_UP_AOI.tif")

# extract variables to points

Vector_variables<-extract(Stack_Set_1,Vector,df=TRUE)
Vector_variables$ID <- Vector$cell
Vector_variables <- Vector_variables[order(Vector_variables$ID), ]
idx <- complete.cases(Vector_variables)
Vector_variables <- Vector_variables[idx, ]
saveRDS(Vector_variables, "vv.rds")


# Create A vector to save the results

C_INPUT_EQ<-Vector[idx, 1]


# use this only for backup

# C_INPUT_EQ<-readOGR("OUTPUTS/1_SPIN_UP/SPIN_UP_BSAS_27-03-2020_332376.shp")


# Extract the layers from the Vector

SOC_im<-Vector_variables[[2]] # primera banda del stack

clay_im<-Vector_variables[[3]] # segunda banda del stack 

DR_im<-Vector_variables[[40]]

LU_im<-Vector_variables[[41]]

# Define Years for Cinputs calculations

years=seq(1/12,500,by=1/12)


########################################################
# calculates some iom in t / ha
# input
# 1. c total carbpn stock in t /ha
#####################################################
fIOM.Falloon.RothC =function(c, par1=-1.31, par2=1.139)
{
  
  # IOM=10^(par1+par2*log10(c))
  IOM=0.049*SOC^(1.139) 
  IOM
}


#################################################################################
# fget_equilibrium_fractions.RothC_input 
# brief: quantifies pool distribution and C input for RothC at equilibrium
#Input
# xi= scalar representing an averaged modifying factor
# C.tot = initial C stock (and C stock in equilibrium)
# clay = clay content
# fractI = vector of Cinput fractions that enter the DPM, RPM, HUM  
#          with a DR of 1.44 fractI becomes [1] 0.5901639 0.4098361 0.0000000
#          by fractI=c((DR)/(DR+1),1-(DR)/(DR+1),0)
#Output
# list with pools at equilibrium and C input at equilibrium
################################################################################

fget_equilibrium_fractions.RothC_input=function(xi=1,C.tot,clay, fractI)
{   
  rmf=xi
  IOM= fIOM.Falloon.RothC(c = C.tot)
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
  # IOM
  ###################################################
  fract.all_stock=(fract.all*C.active)
  fract.all=fract.all_stock/C.tot
  fract.all=append(fract.all,IOM/C.tot)
  pools=fract.all*C.tot
  
  Cin=(C.tot-pools[5])/denominator
  list(pools,Cin)
}

# ROTH C MODEL FUNCTION . 

########## function set up starts###############
Roth_C_equi_analy<-function(Cinputs,Temp,Precip,Evp,Cov2,soil.thick,SOC,clay,DR,bare1,LU)
  
{
  
  #------------------------
  # c input distribution
  #-------------------------
  fractI=c((DR)/(DR+1),1-(DR)/(DR+1),0)
  
  # Paddy fields coefficent fPR = 0.4 if the target point is class = 13 , else fPR=1
  # From Shirato and Yukozawa 2004
  
  
  fPR=(LU == 13)*0.4 + (LU!=13)*1
  
  #Temperature effects per month
  fT=fT.RothC(Temp[,2]) 
  
  #Moisture effects per month . 
  
  fw1func<-function(P, E, S.Thick = 30, pClay = 32.0213, pE = 1, bare) 
  {
    
    M = P - E * pE
    Acc.TSMD = NULL
    for (i in 2:length(M)) {
      B = ifelse(bare[i] == FALSE, 1, 1.8)
      Max.TSMD = -(20 + 1.3 * pClay - 0.01 * (pClay^2)) * (S.Thick/23) * (1/B)
      Acc.TSMD[1] = ifelse(M[1] > 0, 0, M[1])
      if (Acc.TSMD[i - 1] + M[i] < 0) {
        Acc.TSMD[i] = Acc.TSMD[i - 1] + M[i]
      }
      else (Acc.TSMD[i] = 0)
      if (Acc.TSMD[i] <= Max.TSMD) {
        Acc.TSMD[i] = Max.TSMD
      }
    }
    b = ifelse(Acc.TSMD > 0.444 * Max.TSMD, 1, (0.2 + 0.8 * ((Max.TSMD - 
                                                                Acc.TSMD)/(Max.TSMD - 0.444 * Max.TSMD))))
    b<-clamp(b,lower=0.2)
    return(data.frame(b))   
  }
  
  fW_2<- fw1func(P=(Precip[,2]), E=(Evp[,2]), S.Thick = soil.thick, pClay = clay, pE = 1, bare=bare1)$b 
  
  #Vegetation Cover effects 
  
  fC<-Cov2[,2]
  
  # Set the factors frame for Model calculations
  
  xi=mean(fT*fW_2*fC*fPR)
  
  # RUN THE MODEL 
  
  result=fget_equilibrium_fractions.RothC_input(xi=xi,C.tot=SOC,clay=clay, fractI)
  
  # RUN THE MODEL FROM SOILR
  #Model3_spin=RothCModel(t=years,C0=c(DPMptf, RPMptf, BIOptf, HUMptf, FallIOM),In=Cinputs,DR=DR,clay=clay,xi=xi.frame, pass=TRUE) 
  #Ct3_spin=getC(Model3_spin)
  
  # Get the final pools of the time series
  
  
  return(result)
}

######### function set up ends###############

# Iterates over the area of interest
#source("D:/projecte/Rlibs/dataframe_ops.R")

########for loop starts###############3
for (i in 1:dim(Vector_variables)[1]) {
  
  
  # Extract the variables 
  
  Vect<-as.data.frame(Vector_variables[i,])
  
  Temp<-as.data.frame(t(Vect[4:15]))
  Temp<-data.frame(Month=1:12, Temp=Temp[,1])
  
  Precip<-as.data.frame(t(Vect[16:27]))
  Precip<-data.frame(Month=1:12, Precip=Precip[,1])
  
  Evp<-as.data.frame(t(Vect[28:39]))
  Evp<-data.frame(Month=1:12, Evp=Evp[,1])
  
  Cov<-as.data.frame(t(Vect[42:53]))
  Cov1<-data.frame(Cov=Cov[,1])
  Cov2<-data.frame(Month=1:12, Cov=Cov[,1])
  
  #Avoid calculus over Na values 
  
  if (any(is.na(Evp[,2])) | any(is.na(Temp[,2])) | any(is.na(SOC_im[i])) | any(is.na(clay_im[i])) | 
      any(is.na(Precip[,2]))  |  any(is.na(Cov2[,2]))  |  any(is.na(Cov1[,1]))  | any(is.na(DR_im[i])) |  
      (SOC_im[i]<0) | (clay_im[i]<0) ) {C_INPUT_EQ[i,2]<-0
  }else{
    
    # Set the variables from the images
    
    soil.thick=30  #Soil thickness (organic layer topsoil), in cm
    SOC<-SOC_im[i]      #Soil organic carbon in Mg/ha 
    clay<-clay_im[i]        #Percent clay %
    
    DR<-DR_im[i]              # DPM/RPM (decomplosable vs resistant plant material.)
    bare1<-(Cov1>0.8)           # If the surface is bare or vegetated
    LU<-LU_im[i]
    
    #IOM using Falloon method
    FallIOM=0.049*SOC^(1.139) 
    
    # If you use a SOC uncertainty layer turn on this. First open the layer SOC_UNC 
    #(it must have the same extent and resolution of the SOC layer)
    
    #SOC_min<-(1-(SOC_UNC/100))*SOC
    #SOC_max<-(1+(SOC_UNC/100))*SOC
    
    # Define SOC min, max Clay min and max. 
    SOC_min<-SOC*0.8
    SOC_max<-SOC*1.2
    clay_min<-clay*0.9
    clay_max<-clay*1.1
    
    ##############################################################################  
    # C input equilibrium. (Ceq) + Ceq_MIN + Ceq_MAX are quantified here
    ##############################################################################  
    
    #fb<-Roth_C(Cinputs=b,years=years,DPMptf=0, RPMptf=0, BIOptf=0, HUMptf=0, FallIOM=FallIOM,Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    #fb_t<-fb[1]+fb[2]+fb[3]+fb[4]+fb[5]
    #pool.equi.goodi[i,]=fb
    #m<-(fb_t-FallIOM)/(b)
    
    #Ceq<-(SOC-FallIOM)/m
    
    #Cin.equi$spinup[i]=Ceq;
    
    result=Roth_C_equi_analy(Cinputs=b,Temp=Temp,Precip=Precip,Evp=Evp,Cov2=Cov2,soil.thick,SOC,clay,DR,bare1,LU)
    Ceq = result[[2]]
    pool.equi.mean = result[[1]]
    
    # UNCERTAINTIES C input equilibrium (MINIMUM)
    
    result=Roth_C_equi_analy(Cinputs=b,Temp=Temp,Precip=Precip,Evp=Evp,Cov2=Cov2,soil.thick,SOC_min,clay_min,DR,bare1,LU)
    Ceq_MIN = result[[2]]
    pool.equi.min = result[[1]]
    
    # UNCERTAINTIES C input equilibrium (MAXIMUM)
    
    result=Roth_C_equi_analy(Cinputs=b,Temp=Temp,Precip=Precip,Evp=Evp,Cov2=Cov2,soil.thick,SOC_max,clay_max,DR,bare1,LU)
    Ceq_MAX = result[[2]]
    pool.equi.max = result[[1]]
    
    # SOC POOLS AFTER 500 YEARS RUN WITH C INPUT EQUILIBRIUM
    good_landuse_classes=c(2,12,13,4,3,5,6,8)
    if (LU %in% good_landuse_classes){
      
      C_INPUT_EQ[i,2]<-SOC
      C_INPUT_EQ[i,3]<-Ceq
      C_INPUT_EQ[i,4]<-sum(pool.equi.mean)
      C_INPUT_EQ[i,5]<-pool.equi.mean[1] #DPM
      C_INPUT_EQ[i,6]<-pool.equi.mean[2] #RPM
      C_INPUT_EQ[i,7]<-pool.equi.mean[3] #BIO
      C_INPUT_EQ[i,8]<-pool.equi.mean[4] #HUM
      C_INPUT_EQ[i,9]<-pool.equi.mean[5] #IOM
      C_INPUT_EQ[i,10]<-Ceq_MIN
      C_INPUT_EQ[i,11]<-Ceq_MAX
      C_INPUT_EQ[i,12]<-sum(pool.equi.min)
      C_INPUT_EQ[i,13]<-pool.equi.min[1] #DPM
      C_INPUT_EQ[i,14]<-pool.equi.min[2] #RPM
      C_INPUT_EQ[i,15]<-pool.equi.min[3] #BIO
      C_INPUT_EQ[i,16]<-pool.equi.min[4] #HUM
      C_INPUT_EQ[i,17]<-pool.equi.min[5] #IOM
      C_INPUT_EQ[i,18]<-sum(pool.equi.max)
      C_INPUT_EQ[i,19]<-pool.equi.max[1] #DPM 
      C_INPUT_EQ[i,20]<-pool.equi.max[2] #RPM 
      C_INPUT_EQ[i,21]<-pool.equi.max[3] #BIO 
      C_INPUT_EQ[i,22]<-pool.equi.max[4] #HUM 
      C_INPUT_EQ[i,23]<-pool.equi.max[5] #IOM 
      
    }else {
      C_INPUT_EQ[i,2]<-SOC
      C_INPUT_EQ[i,3]<-Ceq
      C_INPUT_EQ[i,4]<-0
      C_INPUT_EQ[i,5]<-0
      C_INPUT_EQ[i,6]<-0
      C_INPUT_EQ[i,7]<-0
      C_INPUT_EQ[i,8]<-0
      C_INPUT_EQ[i,9]<-0
      C_INPUT_EQ[i,10]<-0
      C_INPUT_EQ[i,11]<-0
      C_INPUT_EQ[i,12]<-0
      C_INPUT_EQ[i,13]<-0
      C_INPUT_EQ[i,14]<-0
      C_INPUT_EQ[i,15]<-0
      C_INPUT_EQ[i,16]<-0
      C_INPUT_EQ[i,17]<-0
      C_INPUT_EQ[i,18]<-0
      C_INPUT_EQ[i,19]<-0
      C_INPUT_EQ[i,20]<-0
      C_INPUT_EQ[i,21]<-0
      C_INPUT_EQ[i,22]<-0
      C_INPUT_EQ[i,23]<-0
      
    }
    print(c(i,SOC,Ceq))
    
  } # NA problems
}
###############for loop ends##############

#rename de columns

colnames(C_INPUT_EQ@data)[2]="SOC_FAO"
colnames(C_INPUT_EQ@data)[3]="Cinput_EQ"
colnames(C_INPUT_EQ@data)[4]="SOC_pedotransfer"
colnames(C_INPUT_EQ@data)[5]="DPM_pedotransfer"
colnames(C_INPUT_EQ@data)[6]="RPM_pedotransfer"
colnames(C_INPUT_EQ@data)[7]="BIO_pedotransfer"
colnames(C_INPUT_EQ@data)[8]="HUM_pedotransfer"
colnames(C_INPUT_EQ@data)[9]="IOM_pedotransfer"
colnames(C_INPUT_EQ@data)[10]="CIneq_min"
colnames(C_INPUT_EQ@data)[11]="CIneq_max"
colnames(C_INPUT_EQ@data)[12]="SOC_min"
colnames(C_INPUT_EQ@data)[13]="DPM_min"
colnames(C_INPUT_EQ@data)[14]="RPM_min"
colnames(C_INPUT_EQ@data)[15]="BIO_min"
colnames(C_INPUT_EQ@data)[16]="HUM_min"
colnames(C_INPUT_EQ@data)[17]="IOM_min"
colnames(C_INPUT_EQ@data)[18]="SOC_max"
colnames(C_INPUT_EQ@data)[19]="DPM_max"
colnames(C_INPUT_EQ@data)[20]="RPM_max"
colnames(C_INPUT_EQ@data)[21]="BIO_max"
colnames(C_INPUT_EQ@data)[22]="HUM_max"
colnames(C_INPUT_EQ@data)[23]="IOM_max"

# SAVE the Points (shapefile)

# setwd("D:/TRAINING_MATERIALS_GSOCseq_MAPS_12-11-2020/OUTPUTS/1_SPIN_UP")
writeOGR(C_INPUT_EQ, ".", "SPIN_UP_Country_AOI", driver="ESRI Shapefile",overwrite=TRUE)

