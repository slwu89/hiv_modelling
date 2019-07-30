
#Uganda SHARE-trial model
rm(list=ls());gc()
library(FME)
#library(ggplot2)
pars=list(T_f=100000, T_m=100000, beta_FM=0.002, beta_MF=0.004, theta=25, circum=0.60, phi_0=0.12, phi=0.0,phi_v=0,treat=0.015,
          f_FGP=0.02, f_FBB=0.64, f_FTS=0.2,f_mf=0.31, f_mc=0.16, gamma_V=7, gamma=0.125,gamma_t=0.031,
          mu = 0.017,p0_FGP=0.08, p0_MGP=0.04, p0_FTS=0.2, p0_CTS=0.25,p0_CBB=0.2, p0_FBB=0.39,P_FBB=0.005,P_FTS=0.03,
          pmf=0.12,pmc=0.38,df_BB=2,df_TS=2, dm_BB=2,dm_TS=2,Z=0.6,CFBB_CBB1=48,CFTS_CTS1=2,CCBB_FBB=295,CMGP_FGP=1,CMGP_FGP1=1,
          CCTS_FTS=6,CCBB_FTS=6,eta=0.04,d=0.02, tau=0.14,ep=0.8,psi_v1=0.8,psi_v2=1.11)
FSW.ode<- function(pars,times=seq(1979,2011,by=1))
{derivs <- function(t,state,pars)
  
{ # returns rate of change
  with(as.list(c(state, pars)), 
{
  # browser()
#circumcision
  psi=1-circum
  if(t<2002){
    phi=0
  }else if(t>=2002&t<=2006){
    phi=0.04
  }else {
    phi=0.01
  }
  
#  #Treatment rates
   
  if(t<2004){
    tau =0
  } else if(t>=2004&t<=2012){
    tau=(t-2004)*treat
   }
# else {
#      tau=0.14
#    }
#violence terms
if(t<2005){
  phi_v=0
}else if(t>=2005&t<=2006){
  phi_v=0.68
}else {
  phi_v=0.04
}

if(t<2006){
  pci =1
} else if(t>=2006&t<=2009){
  pci=psi_v1
}else {
  pci=psi_v2
}


  #Duration 
  #t_gamma=1/t_du #treatment rate
  alpha_BB=1/df_BB
  alpha_TS=1/df_TS
  delta_BB=1/dm_BB
  delta_TS=1/dm_TS
  
  #Population sizes
  #initial sizes
  N0_FGP=T_f-(P_FBB*T_f+P_FTS*T_f) 
  N0_MGP=abs(T_m-(CCBB_FBB*P_FBB*T_f/CFBB_CBB1+(1-Z)*CCTS_FTS*P_FTS*T_f/CFTS_CTS1))
  N0_MGPint=T_m-(CCBB_FBB*P_FBB*T_f/CFBB_CBB1+(1-Z)*CCTS_FTS*P_FTS*T_f/CFTS_CTS1)
  N0_MGP=(T_f-(P_FBB*T_f+P_FTS*T_f))
  N0_CBB=CCBB_FBB*P_FBB*T_f/CFBB_CBB1
  N0_CTS=(1-Z)*CCTS_FTS*P_FTS*T_f/CFTS_CTS1
  N0_FBB=P_FBB*T_f
  N0_FTS=P_FTS*T_f
  
  #variable population sizes
  N_FGP=S_FGP+S_FGI+I_FGV+I_FGP+T_FGP
  N_MGP=S_MGP+S_MGC+I_MGV+I_MGP+T_MGP
  
  N_FBB=S_FBB+S_FBI+I_FBV+I_FBB+T_FBB
  N_FTS=S_FTS+S_FTI+I_FTV+I_FTS+T_FTS
  
  N_CBB=S_CBB+S_CBC+I_CBV+I_CBB+T_CBB
  N_CTS=S_CTS+S_CTC+I_CTV+I_CTS+T_CTS
  
  N_FSW=N_FBB+N_FTS
  N_C=N_CBB+N_CTS
  N_F=N_FGP+N_FSW
  N_M=N_MGP+N_C
  
 
  P_MGPFGP=N_FGP/(pmf*N_FBB+N_FGP)
  P_MGPFBB=pmf*N_FBB/(pmf*N_FBB+N_FGP)
  P_FGPMGP=N_MGP/(pmc*N_CBB+N_MGP)
  P_FGPCBB=pmc*N_CBB/(pmc*N_CBB+N_MGP)
 
  #contact constraints
  CFBB_CBB=CCBB_FBB*N_FBB/N_CBB
  CFTS_CTS=(1-Z)*(CCTS_FTS*N_FTS)/N_CTS
  CFTS_CBB=Z*CCTS_FTS*N_FTS/N_CBB
  CFGP_MGP=P_FGPMGP*(CMGP_FGP*N_FGP)/(P_MGPFGP*N_MGP)
  CMGP_FBB=P_MGPFBB*(CFGP_MGP*N_MGP)/(pmf*N_FBB)
  CFGP_CBB=P_FGPCBB*(CMGP_FGP*N_FGP)/(pmc*N_CBB)
 
  #General population
  lambda_MGP=beta_FM*CFGP_MGP*54*((1-ep*f_FGP)*P_MGPFBB*(theta*I_FBV+I_FBB+eta*T_FBB)/N_FBB+(1-ep*f_FGP)*P_MGPFGP*(theta*I_FGV+I_FGP+eta*T_FGP)/N_FGP)
  lambda_FGP=beta_MF*CMGP_FGP*54*((1-ep*f_FGP)*P_FGPCBB*(theta*I_CBV+I_CBB+eta*T_CBB)/N_CBB+(1-ep*f_FGP)*P_FGPMGP*(theta*I_MGV+I_MGP+eta*T_MGP)/N_MGP)
 
  #clients
  lambda_CBB=beta_FM*(CFBB_CBB*15*(1-ep*f_FBB)*(theta*I_FBV+I_FBB+eta*T_FBB)/N_FBB + CFTS_CBB*20*(1-ep*f_FTS)*(theta*I_FTV+I_FTS+eta*T_FTS)/N_FTS+CFGP_CBB*54*(1-ep*f_FGP)*(theta*I_FGV+I_FGP+eta*T_FGP)/N_FGP)
  lambda_CTS=beta_FM*CFTS_CTS*20*(1-ep*f_FTS)*(theta*I_FTV+I_FTS+eta*T_FTS)/N_FTS
  
  #sex workers
  lambda_FBB=beta_MF*(CCBB_FBB*15*(1-ep*f_FBB)*(theta*I_CBV+I_CBB+eta*T_CBB)/N_CBB+CMGP_FBB*54*(1-ep*f_FGP)*(theta*I_MGV+I_MGP+eta*T_MGP)/N_MGP)
  lambda_FTS=beta_MF*CCTS_FTS*20*(1-ep*f_FTS)*(Z*(theta*I_CBV+I_CBB+eta*T_CBB)/N_CBB+(1-Z)*(theta*I_CTV+I_CTS+eta*T_CTS)/N_CTS)
    
  #Tinci=(lambda_FGP*S_FGP+lambda_MGP*S_MGP+psi*lambda_MGP*S_MGC+lambda_CTS*S_CTS+psi*lambda_CTS*S_CTC+lambda_CBB*S_CBB+psi*lambda_CBB*S_CBC+lambda_FTS*S_FTS+lambda_FBB*S_FBB)*100
  #General population
  
  #Females
  dS_FGP = p0_FGP*N0_FGP+alpha_BB*S_FBB+alpha_TS*S_FTS-lambda_FGP*S_FGP-(mu+phi_v)*S_FGP #HIV -ve individuals 
  dS_FGI = phi_v*S_FGP+alpha_BB*S_FBI+alpha_TS*S_FTI-pci*lambda_FGP*S_FGI-mu*S_FGI
  dI_FGV= lambda_FGP*S_FGP+pci*lambda_FGP*S_FGI-(gamma_V+mu)*I_FGV                       #HIV +ve individuals
  dI_FGP=gamma_V*I_FGV+alpha_BB*I_FBB+alpha_TS*I_FTS-(mu+gamma+tau)*I_FGP+d*T_FGP
  dT_FGP=tau*I_FGP+delta_BB*T_FBB+delta_TS*T_FTS-(mu+d+gamma_t)*T_FGP #HIV treatment class
  #Males
  dS_MGP =(1-phi_0)*p0_MGP*N0_MGP+delta_BB*S_CBB+delta_TS*S_CTS-lambda_MGP*S_MGP-(phi+mu)*S_MGP      #HIV -ve individuals
  dS_MGC = phi_0*p0_MGP*N0_MGP+phi*S_MGP+delta_BB*S_CBC+delta_TS*S_CTC-psi*lambda_MGP*S_MGC-mu*S_MGC
  dI_MGV = lambda_MGP*S_MGP+psi*lambda_MGP*S_MGC-(gamma_V+mu)*I_MGV                                  #HIV +ve individuals
  dI_MGP = gamma_V*I_MGV+delta_BB*I_CBB+delta_TS*I_CTS-(mu+gamma+tau)*I_MGP+d*T_MGP
  dT_MGP = tau*I_MGP+delta_BB*T_CBB+delta_TS*T_CTS-(mu+d+gamma_t)*T_MGP #HIV treatment class
  #High risk males
    
  #Female 2+ partners
  dS_CTS=(1-phi_0)*p0_CTS*N0_CTS-lambda_CTS*S_CTS-(delta_TS+phi+mu)*S_CTS
  dS_CTC=phi_0*p0_CTS*N0_CTS+phi*S_CTS-psi*lambda_CTS*S_CTC-(delta_TS+mu)*S_CTC
  dI_CTV=lambda_CTS*S_CTS+psi*lambda_CTS*S_CTC-(gamma_V+mu)*I_CTV
  dI_CTS=gamma_V*I_CTV-(mu+gamma+tau+delta_TS)*I_CTS+d*T_CTS
  dT_CTS = tau*I_CTS-(mu+d+gamma_t+delta_TS)*T_CTS #HIV treatment class
  #FSW clients
  dS_CBB=(1-phi_0)*p0_CBB*N0_CBB-lambda_CBB*S_CBB-(delta_BB+phi+mu)*S_CBB
  dS_CBC=phi_0*p0_CBB*N0_CBB+ phi*S_CBB-psi*lambda_CBB*S_CBC-(delta_BB+mu)*S_CBC
  dI_CBV=lambda_CBB*S_CBB+psi*lambda_CBB*S_CBC-(gamma_V+mu)*I_CBV
  dI_CBB=gamma_V*I_CBV-(mu+gamma+tau+delta_BB)*I_CBB+d*T_CBB
  dT_CBB = tau*I_CBB-(mu+d+gamma_t+delta_BB)*T_CBB #HIV treatment class
 
  #High risk females
  
  #Female 2+
  dS_FTS=p0_FTS*N0_FTS-lambda_FTS*S_FTS-(alpha_TS+mu+phi_v)*S_FTS
  dS_FTI=phi_v*S_FTS-pci*lambda_FTS*S_FTI-(alpha_TS+mu)*S_FTI
  dI_FTV=lambda_FTS*S_FTS+pci*lambda_FTS*S_FTI-(gamma_V+mu)*I_FTV
  dI_FTS=gamma_V*I_FTV-(mu+gamma+tau+alpha_TS)*I_FTS+d*T_FTS
  dT_FTS = tau*I_FTS-(mu+d+gamma_t+alpha_TS)*T_FTS
  #FSWs
  dS_FBB=p0_FBB*N0_FBB-lambda_FBB*S_FBB-(alpha_BB+mu+phi_v)*S_FBB
  dS_FBI=phi_v*S_FBB-pci*lambda_FBB*S_FBI-(alpha_BB+mu)*S_FBI
  dI_FBV=lambda_FBB*S_FBB+pci*lambda_FBB*S_FBI-(gamma_V+mu)*I_FBV
  dI_FBB=gamma_V*I_FBV-(mu+gamma+tau+alpha_BB)*I_FBB+d*T_FBB
  dT_FBB = tau*I_FBB-(mu+d+gamma_t+alpha_BB)*T_FBB #HIV treatment class
  
  
  return(list(c(dS_FGP,dS_FGI,dI_FGV,dI_FGP,dT_FGP,dS_MGP,dS_MGC,dI_MGV,dI_MGP,dT_MGP, dS_CBB,dS_CBC,dI_CBV,dI_CBB,dT_CBB, dS_CTS,dS_CTC,dI_CTV,dI_CTS,dT_CTS,dS_FBB,dS_FBI,dI_FBV,dI_FBB,dT_FBB,dS_FTS,dS_FTI,dI_FTV,dI_FTS,dT_FTS), 
              N_F=N_FGP+N_FSW,N_M=N_MGP+N_C,
              N0_MGPint=T_m-(CCBB_FBB*P_FBB*T_f/CFBB_CBB1+(1-Z)*CCTS_FTS*P_FTS*T_f/CFTS_CTS1),
              N0_CBB=CCBB_FBB*P_FBB*T_f/CFBB_CBB1,
              N0_CTS=(1-Z)*CCTS_FTS*P_FTS*T_f/CFTS_CTS1,
              Prev_F=((I_FGP+I_FGV+T_FGP+I_FBB+I_FBV+T_FBB+I_FTS+I_FTV+T_FTS)/N_F)*100,
              Prev_M=((I_MGP+I_MGV+T_MGP+I_CBB+I_CBV+T_CBB+I_CTS+I_CTV+T_CTS)/N_M)*100,
              Prev_CBB=((I_CBB+I_CBV)/N_CBB)*100,
              Prev_CTS=((I_CTS+I_CTV)/N_CTS)*100,Prev_FBB=((I_FBB+I_FBV)/N_FBB)*100,
              CFBB_CBB=CCBB_FBB*N_FBB/N_CBB,
              CFTS_CTS=(1-Z)*(CCTS_FTS*N_FTS)/N_CTS, 
              CFTS_CBB=Z*CCTS_FTS*N_FTS/N_CBB,
              CFGP_MGP=P_FGPMGP*(CMGP_FGP*N_FGP)/(N_MGP),
              CMGP_FBB=P_MGPFBB*(CFGP_MGP*N_MGP)/(pmf*N_FBB),
              CFGP_CBB=P_FGPCBB*(CMGP_FGP*N_FGP)/(pmc*N_CBB),
              P_MGPFGP=N_FGP/(pmf*N_FBB+N_FGP),
              P_MGPFBB=pmf*N_FBB/(pmf*N_FBB+N_FGP),
              P_FGPMGP=N_MGP/(pmc*N_CBB+N_MGP),
              P_FGPCBB=pmc*N_CBB/(pmc*N_CBB+N_MGP),
              N_FGP=S_FGP+S_FGI+I_FGV+I_FGP+T_FGP,
              N_MGP=S_MGP+S_MGC+I_MGV+I_MGP+T_MGP,
              N_FBB=S_FBB+S_FBI+I_FBV+I_FBB+T_FBB,
              N_FTS=S_FTS+S_FTI+I_FTV+I_FTS+T_FTS,
              N_CBB=S_CBB+S_CBC+I_CBV+I_CBB+T_CBB,
              N_CTS=S_CTS+S_CTC+I_CTV+I_CTS+T_CTS,
              P_VF=(S_FGI+S_FBI+S_FTI)*100/(S_FGP+S_FGI+S_FBB+S_FBI+S_FTS+S_FTI),
              T_MC=(S_MGC+S_CBC+S_CTC)*100/(S_MGP+S_MGC+S_CBB+S_CBC+S_CTS+S_CTC),
              T_TR=(T_FGP+T_FBB+T_FTS+T_MGP+T_CBB+T_CTS)*100/(I_FGP+I_FBB+I_FTS+I_MGP+I_CBB+I_CTS+T_FGP+T_FBB+T_FTS+T_MGP+T_CBB+T_CTS),
              T_TOT=(T_FGP+T_FBB+T_FTS+T_MGP+T_CBB+T_CTS),
              Tinci=(lambda_FGP*S_FGP+pci*lambda_FGP*S_FGI+lambda_MGP*S_MGP+
                       psi*lambda_MGP*S_MGC+lambda_CTS*S_CTS+psi*lambda_CTS*S_CTC+
                       lambda_CBB*S_CBB+psi*lambda_CBB*S_CBC+lambda_FTS*S_FTS+pci*
                       lambda_FTS*S_FTI+lambda_FBB*S_FBB+pci*lambda_FBB*S_FBI)*100/
                (S_FGP+S_FGI+S_MGP+S_MGC+S_CTS+S_CTC+S_CBB+S_CBC+S_FTS+S_FTI+S_FBB+S_FBI)))})
}
 #initial conditions
  state<-c(S_FGP=pars$T_f-(pars$P_FBB*pars$T_f+pars$P_FTS*pars$T_f),S_FGI=0, 
           I_FGV=0,I_FGP=2000,T_FGP=0, S_MGP=(1-pars$phi_0)*abs(pars$T_m-(pars$CCBB_FBB*pars$P_FBB*pars$T_f/pars$CFBB_CBB1+
         (1-pars$Z)*pars$CCTS_FTS*pars$P_FTS*pars$T_f/pars$CFTS_CTS1)),S_MGC=pars$phi_0*abs(pars$T_m-(pars$CCBB_FBB*pars$P_FBB*pars$T_f/pars$CFBB_CBB1+(1-pars$Z)*pars$CCTS_FTS*pars$P_FTS*pars$T_f/pars$CFTS_CTS1)),
          I_MGV=0,I_MGP=1500,T_MGP=0, S_CBB=(1-pars$phi_0)*pars$CCBB_FBB*pars$P_FBB*pars$T_f/pars$CFBB_CBB1,S_CBC=pars$phi_0*pars$CCBB_FBB*pars$P_FBB*pars$T_f/pars$CFBB_CBB1,I_CBV=0,I_CBB=100,T_CBB=0, S_CTS=(1-pars$phi_0)*
           (1-pars$Z)*pars$CCTS_FTS*pars$P_FTS*pars$T_f/pars$CFTS_CTS1,S_CTC=pars$phi_0*(1-pars$Z)*pars$CCTS_FTS*pars$P_FTS*pars$T_f/pars$CFTS_CTS1,I_CTV=0,I_CTS=100,T_CTS=0,S_FBB=pars$P_FBB*pars$T_f,S_FBI=0,I_FBV=0,I_FBB=1000,T_FBB=0, S_FTS=pars$P_FTS*pars$T_f,S_FTI=0,I_FTV=0,I_FTS=100, T_FTS=0) 
# ode solves the model by integration
 return(as.data.frame(ode(y = state,times = times, func = derivs, parms = pars)))
}

out=FSW.ode(pars)

#### up til here ######

#Time series plots
R0<-cbind(out$N_F,out$N_M)
R1<-cbind(out$S_FGP,out$I_FGP,out$S_MGP,out$I_MGP)
R2<-cbind(out$S_CBB,out$I_CBB,out$S_CTS,out$I_CTS)
R3<-cbind(out$S_FBB,out$I_FBB,out$S_FTS,out$I_FTS)


#Prevalences
R4<-cbind(out$Prev_F,out$Prev_M)
R5<-cbind(out$Prev_CBB,out$Prev_CTS)
R6<-cbind(out$Prev_FBB,out$Prev_FTS)

#High Viraemia
R7<-cbind(out$S_FGP,out$I_FGV,out$S_MGP,out$I_MGV)
R8<-cbind(out$S_CBB,out$I_CBV,out$S_CTS,out$I_CTV)
R9<-cbind(out$S_FBB,out$I_FBV,out$S_FTS,out$I_FTV)


t<-out[,1]

par(mfrow=c(1,1))
matplot(t,R0,type="l",lwd = c(2,2), col =c("blue","red"), xlab = "Time (years)", ylab = "Population size")
legend("topright", c("N_F","N_M"),lty = 1:2,bty="n",col =c("blue","red"), lwd = c(2,2),cex = 0.6)


par(mfrow=c(1,1))
matplot(t,R1,type="l",lwd = c(2,2,2,2), col =c("blue","red","black","brown"), xlab = "Time (years)", ylab = "Population size")
legend("center", c("S_FGP","I_FGP", "S_MGP","I_MGP"),lty = 1:4,bty="n",col =c("blue","red","black","yellow"), lwd = c(2,2,2,2),cex = 0.6)

par(mfrow=c(1,2))
matplot(t,R2,type="l",lwd = c(2,2,2,2), col =c("forestgreen","green","blue","red"), xlab = "Time (years)", ylab = "Populations size")
legend("topright", c("S_CBB","I_CBB","S_CTS", "I_CTS"),lty = 1:6,bty="n",col=c("forestgreen","green","blue","red"), lwd =c(2,2,2,2),cex = 0.6)

matplot(t,R3,type="l",lwd = c(2,2,2,2), col =c("forestgreen","green","blue","red"), xlab = "Time (years)", ylab = "Population size")
legend("topright", c("S_FBB","I_FBB","S_FTS", "I_FTS"),lty = 1:6,bty="n",col=c("forestgreen","green","blue","red"), lwd = c(2,2,2,2),cex = 0.6)
#Prevtes
par(mfrow=c(1,1))
matplot(t,R4,type="l",lwd = c(2,2), col =c("blue","red"), xlab = "Time (years)", ylab = "Prevalence")
legend("topright", c("Prev_F","Prev_M"),lty = 1:2,bty="n",col =c("blue","red"), lwd = c(2,2),cex = 0.6)

par(mfrow=c(1,2))
matplot(t,R5, type="l",lwd = c(2,2), col =c("forestgreen","green"), xlab = "Time (years)", ylab = "Prevalence")
legend("topright", c("Prev_CBB","Prev_CTS"),lty = 1:3,bty="n",col=c("forestgreen","green"), lwd =c(2,2),cex = 0.6)

matplot(t,R6,type="l",lwd = c(2,2,2), col =c("forestgreen","blue"), xlab = "Time (years)", ylab = "Prevalence")
legend("topright", c("Prev_FBB","Prev_FTS"),lty = 1:3,bty="n", col=c("forestgreen","green"), lwd = c(2,2),cex = 0.6)
#write.table(out,"out.txt",sep="\t")

# parameter ranges
#v     "beta_FM","beta_MF","P_FBB", "P_FTS",  "pmf"    "pmc"  "f_mf"  "f_mc"   "df_BB", "df_TS",  "dm_BB", "dm_TS"  "Z" , CCBB_FBB,"CFBB_CBB1" "CCTS_FTS", f_FGP, "f_FBB" "f_FTS" "circum" "treat",  "gamma_V","theta","eta",   "d",      ep, psi_v1,psi_v2)
min<-c(0.0006,   0.0006,   0.001,    0.001,   0.1,     0.28,   0.26,    0.8  ,   1,       5,        5,        1 ,     0,   200   ,   12,         3,       0.05,   0.6,     0.10,    0.5,     0.01,     4.8,     20,    0,    0.01,     0.8, 0.42,0.79)
max<-c(0.0015,   0.004,    0.006,    0.30,    0.15,    0.50,   0.47,   0.32 ,    6,       15,       10,       10 ,    1,    1200  ,    48,         15,      0.07,  0.7,     0.39,    0.6,    0.02,     7.3,    25 ,   0.16,   0.07,    0.95,0.97,1.55)
parRanges<-cbind(min,max)
rownames(parRanges)= c("beta_FM","beta_MF","P_FBB","P_FTS","pmf","pmc","f_mf","f_mc","df_BB","df_TS","dm_BB","dm_TS","Z","CCBB_FBB","CFBB_CBB1","CCTS_FTS","f_FGP","f_FBB","f_FTS","circum","treat","gamma_V", "theta","eta","d","ep","psi_v1","psi_v2")
parRanges

# Sensitivity analysis of all parameters (LHS)
SA0=sensRange(func =FSW.ode, parms = pars, dist = "latin",sensvar = c("N_F","N_M","Prev_F","Prev_M","T_MC","T_TR","N_FGP","N_MGP","N_FBB","N_CBB","N_CTS","N_FTS","Tinci","P_MGPFGP","P_MGPFBB","P_FGPCBB","S_MGC","T_TOT","P_VF"),parRange = parRanges, num=80000)
write.table(SA0,"SA04.txt",sep="\t")
SA1 <- summary(SA0)
par(mfrow=c(1,1))
plot(SA1, xlab = "time (years)",ylab = "HIV prevalence", mfrow = NULL,quant = F, col = c("lightblue","darkblue"), legpos = "topright")
mtext(outer = TRUE, line = -1.5, side = 1, "",cex = 1.25)

#selecting output parameters

d1 = read.table("SA04.txt", sep="\t",header = T)

#d3=subset(d1,(Prev_F2005>=6&Prev_F2005<=14)&(Prev_M2005>=6&Prev_M2005<=14))
#d3=subset(d1,(Prev_F2005>=6&Prev_F2005<=14)&(Prev_M2005>=6&Prev_M2005<=14)&(Tinci2005>=1.02&Tinci2005<=1.68)&(Tinci2006>=0.82&Tinci2006<=1.33)&(Tinci2008>=1.02&Tinci2008<=1.54)&(Tinci2009>=1.03&Tinci2009<=1.5)&(Tinci2011>=0.67&Tinci2011<=1.07))
d3=subset(d1,(Prev_F1992>=11&Prev_F1992<=14)&(Prev_M1992>=8&Prev_M1992<=12)&(Prev_F2001>=8&Prev_F2001<=11)&(Prev_M2001>=6&Prev_M2001<=9)&(Prev_F2005>=6&Prev_F2005<=14)&(Prev_M2005>=6&Prev_M2005<=14)&(Tinci2005>=0.62&Tinci2005<=1.22)&(Tinci2006>=0.79&Tinci2006<=1.37)&(Tinci2008>=0.57&Tinci2008<=0.98)&(Tinci2009>=0.65&Tinci2009<=1.04)&(Tinci2011>=0.72&Tinci2011<=1.16))
#d3=subset(d,(Prev_F2005>=0&Prev_F2005<=14)&(Prev_M2005>=0&Prev_M2005<=14)&(Tinci2005>=1.02&Tinci2005<=1.68)&(Tinci2006>=0.82&Tinci2006<=1.68)&(Tinci2008>=0.92&Tinci2008<=1.54)&(Tinci2009>=0.9&Tinci2009<=1.68)&(Tinci2011>=0.67&Tinci2011<=1.68))
#d3=subset(d,(Prev_F2005>=0&Prev_F2005<=140)&(Prev_M2005>=0&Prev_M2005<=140)&(Tinci2005>=0&Tinci2005<=150)&(Tinci2006>=0&Tinci2006<=150)&(Tinci2008>=0 &Tinci2008<=150)&(Tinci2009>=0 &Tinci2009<=150)&(Tinci2011>=0 &Tinci2011<=150))

write.table(d3,"Fits.txt",sep="\t")

SS<-(d3$Tinci2005-1.32)^2+(d3$Tinci2006-1.05)^2+(d3$Tinci2008-0.75)^2+(d3$Tinci2009-0.83)^2+(d3$Tinci2011-0.93)^2
attach(d3) 
d3$sums<- SS 
detach(d3) 
d4=d3[order(d3$sums), ]
d44=subset(d4,select=c(Tinci1979,Tinci1980,Tinci1981,Tinci1982,Tinci1983,Tinci1984,Tinci1985,Tinci1986,Tinci1987,Tinci1988,Tinci1989,Tinci1990, Tinci1991,Tinci1992,Tinci1993,Tinci1994,Tinci1995,Tinci1996,Tinci1997,Tinci1998,Tinci1999,Tinci2000,Tinci2001,Tinci2002,Tinci2003,Tinci2004,Tinci2005,Tinci2006,Tinci2007,Tinci2008,Tinci2009,Tinci2010,Tinci2011))
d4p1=subset(d4,select=c(Prev_F1979,Prev_F1980,Prev_F1981,Prev_F1982,Prev_F1983,Prev_F1984,Prev_F1985,Prev_F1986,Prev_F1987,Prev_F1988,Prev_F1989,Prev_F1990, Prev_F1991,Prev_F1992,Prev_F1993,Prev_F1994,Prev_F1995,Prev_F1996,Prev_F1997,Prev_F1998,Prev_F1999,Prev_F2000,Prev_F2001,Prev_F2002,Prev_F2003,Prev_F2004,Prev_F2005,Prev_F2006,Prev_F2007,Prev_F2008,Prev_F2009,Prev_F2010,Prev_F2011))
d4p2=subset(d4,select=c(Prev_M1979,Prev_M1980,Prev_M1981,Prev_M1982,Prev_M1983,Prev_M1984,Prev_M1985,Prev_M1986,Prev_M1987,Prev_M1988,Prev_M1989,Prev_M1990, Prev_M1991,Prev_M1992,Prev_M1993,Prev_M1994,Prev_M1995,Prev_M1996,Prev_M1997,Prev_M1998,Prev_M1999,Prev_M2000,Prev_M2001,Prev_M2002,Prev_M2003,Prev_M2004,Prev_M2005,Prev_M2006,Prev_M2007,Prev_M2008,Prev_M2009,Prev_M2010,Prev_M2011))


write.table(d4,"Fits2.txt",sep="\t")

d4 = read.table("Fits2.txt", sep="\t",header = T)

d45=subset(d4,select=c(beta_FM,beta_MF,P_FBB,P_FTS,pmf,pmc,f_mf,f_mc,df_BB,df_TS,dm_BB,dm_TS,Z,CCBB_FBB,CFBB_CBB1,CCTS_FTS,f_FGP,f_FBB,f_FTS,circum,treat,gamma_V, theta,eta,d,ep,psi_v1,psi_v2))
write.table(d45,"parameterz.txt",sep="\t")
years <- c(1979,1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,1990, 1991,1992,1993,1994,1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011)
d5 <- rbind(years, d4[1:nrow(d4), ])
d6=d44
d7=d4p1
d8=d4p2
 
par(mfrow = c(1,2)) 
plot(years,d6[1,], type = "l", ylim = c(0,8),  xlab = "Time (years)", ylab = "HIV incidence per 100 person yrs",col = 1) 
for (i in 1:137){
  lines(years,d6[i,], col = topo.colors(137)[i]) 
} 

segments(2005,0.62,2005,1.22)
arrows(2005,0.62,2005,1.22,length=0.03,angle=90,code=3)
points(2005,0.88,pch=0,cex=0.6)

segments(2006,0.79,2006,1.37)
arrows(2006,0.79,2006,1.37,length=0.03,angle=90,code=3)
points(2006,1.05,pch=0,cex=0.6)

segments(2008,0.57,2008,0.98)
arrows(2008,0.57,2008,0.98,length=0.03,angle=90,code=3)
points(2008,0.75,pch=0,cex=0.6)

segments(2009,0.65,2009,1.04)
arrows(2009,0.65,2009,1.04,length=0.03,angle=90,code=3)
points(2009,0.83,pch=0,cex=0.6)

segments(2011,0.72,2011,1.16)
arrows(2011,0.72,2011,1.16,length=0.03,angle=90,code=3)
points(2011,0.93,pch=0,cex=0.6)
plot(years,sapply(d6[-1, ], quantile, 0.50, names=FALSE), type = "l", ylim = c(0, 15),  xlab = "Time (years)", ylab = "HIV incidence per 100 person yrs",col = "blue") 
lines(years,sapply(d6[-1, ], quantile, 0.10, names=FALSE), col = "green") 
lines(years,sapply(d6[-1, ], quantile, 0.25, names=FALSE), col = "brown") 
lines(years,sapply(d6[-1, ], quantile, 0.75, names=FALSE), col = "yellow") 
lines(years,sapply(d6[-1, ], quantile, 0.90, names=FALSE), col = "red") 
legend("topright", c("10 percentile","25 percentile","50 percentile", "75 percentile", "90 percentile"),lty = 1:5,bty="n", col=c("green","brown","blue","yellow","red"), lwd = c(2,2,2,2,2),cex = 0.6)


segments(2005,0.62,2005,1.22)
arrows(2005,0.62,2005,1.22,length=0.03,angle=90,code=3)
points(2005,0.88,pch=0,cex=0.6)

segments(2006,0.79,2006,1.37)
arrows(2006,0.79,2006,1.37,length=0.03,angle=90,code=3)
points(2006,1.05,pch=0,cex=0.6)

segments(2008,0.57,2008,0.98)
arrows(2008,0.57,2008,0.98,length=0.03,angle=90,code=3)
points(2008,0.75,pch=0,cex=0.6)

segments(2009,0.65,2009,1.04)
arrows(2009,0.65,2009,1.04,length=0.03,angle=90,code=3)
points(2009,0.83,pch=0,cex=0.6)

segments(2011,0.72,2011,1.16)
arrows(2011,0.72,2011,1.16,length=0.03,angle=90,code=3)
points(2011,0.93,pch=0,cex=0.6)

#########################################################

plot(years,d7[1,], type = "l", ylim = c(0, 20),  xlab = "Time (years)", ylab = "HIV female prevalence",col = 1) 
for (i in 1:137){
  lines(years,d7[i,], col = topo.colors(137)[i]) 
} 

segments(1992,11,1992,14)
arrows(1992,11,1992,14,length=0.03,angle=90,code=3)
points(1992,13,pch=0,cex=0.6)

segments(2001,8,2001,11)
arrows(2001,8,2001,11,length=0.03,angle=90,code=3)
points(2001,9,pch=0,cex=0.6)

segments(2005,6,2005,14)
arrows(2005,6,2005,14,length=0.03,angle=90,code=3)
points(2005,10,pch=0,cex=0.6)



plot(years,d8[1,], type = "l", ylim = c(0, 20),  xlab = "Time (years)", ylab = "HIV male prevalence",col = 1) 
for (i in 1:137){
  lines(years,d8[i,], col = topo.colors(137)[i]) 
}
segments(1992,8,1992,12)
arrows(1992,8,1992,12,length=0.03,angle=90,code=3)
points(1992,10,pch=0,cex=0.6)

segments(2001,6,2001,9)
arrows(2001,6,2001,9,length=0.03,angle=90,code=3)
points(2001,7,pch=0,cex=0.6)

segments(2005,6,2005,14)
arrows(2005,6,2005,14,length=0.03,angle=90,code=3)
points(2005,10,pch=0,cex=0.6)  


