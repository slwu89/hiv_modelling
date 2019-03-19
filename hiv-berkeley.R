################################################################################
#       ___    ____               __  ___          __     __
#      /   |  / __/_  ______ _   /  |/  /___  ____/ /__  / /
#     / /| | / /_/ / / / __ `/  / /|_/ / __ \/ __  / _ \/ /
#    / ___ |/ __/ /_/ / /_/ /  / /  / / /_/ / /_/ /  __/ /
#   /_/  |_/_/  \__, /\__,_/  /_/  /_/\____/\__,_/\___/_/
#              /____/
#
#   Partow Imani & Sean Wu
#   March 2019
#
################################################################################

rm(list=ls());gc()
library(here)
library(Rcpp)
library(deSolve)

sourceCpp(here("hiv-berkeley.cpp"))

# parameters
pars <- list(NF_0 = 12864738, NM_0 = 12594866, P_transmission = 0.01, HighV_factor = 10,
            T_factor = 0.08, C_FGMG = 1, C_FSWC = 100, C_F2M2 = 3,
            Prop_FSW = 0.08, Prop_MC = 0.088, Prop_IFG = 0.05, Prop_IMG = 0.039, Prop_IFSW = 0.266, Prop_IMC = 0.06,
            Prop_F2 = 0.125, Prop_M2 = 0.125, Prop_IF2 = 0.06, Prop_IM2 = 0.06,
            Tau = 0.242, mu = 0.001,
            mu_I.2 = 0.06, mu_T.2 = 0.01, mu_A.2 = 0.0935, omega = 0.1, rho = 0.08,
            mu_I.1 = 0.1, mu_T.1 = 0.01, mu_A.1 = 0.3,
            con_F2.2 = 0.274, con_FSW.2 = 0.632, con_FG.2 = 0.576, con_MG.2 = 0.587, con_M2.2 = 0.213, con_MC.2 = 0.45,con_eff = 0.875,
            con_F2.1 = 0.1, con_FSW.1 = 0.05, con_FG.1 = 0.05, con_MG.1 = 0.05, con_M2.1 = 0.05, con_MC.1 = 0.05,
            circum.2 = 0.54, circum.1 = 0.8, circum_eff = 0.6, sup.2.rate = 0.08, perc.2nd.line = 0.1, unsup.1.rate = 0.05,
            unsup.2.rate = 0.05, fall.off.2 = 0.06, birth = 0.015, sigma.1 = 0.05, sigma.2 = 0.4,
            sup.1.rate.pre2013 = 0.554, sup.1.rate.pre2015 = 0.634, sup.1.rate.post2015 = 0.554,
            sup.1.rate.pre2015int = 0.85, sup.1.rate.post2015int = 0.749,
            fall.off.1.pre2013 = 0.173, fall.off.1.pre2015 = 0.109, fall.off.1.post2015 = 0.173,
            fall.off.1.pre2015int = 0.009, fall.off.1.post2015int = 0.067)

time <- seq(1990, 2033,by=1)

# initial conditions
NF_0 <- 12864738
NM_0 <- 12594866
state <- c(S_FG=NF_0*(1-pars$Prop_IFG), IV_FG=0 ,I_FG=NF_0*pars$Prop_IFG,T1_FG=0, T1S_FG=0,T2_FG=0, T2S_FG=0,A_FG = 0,
           S_MG=NM_0*(1-pars$Prop_IMG), IV_MG=0, I_MG=NM_0*pars$Prop_IMG, T1_MG=0,T1S_MG=0, T2_MG=0,T2S_MG=0,A_MG = 0,
           S_FSW=NF_0*pars$Prop_FSW, IV_FSW=0, I_FSW=NF_0*pars$Prop_FSW*pars$Prop_IFSW, T1_FSW=0,T1S_FSW=0,T2_FSW=0,T2S_FSW=0, A_FSW=0,
           S_MC=NM_0*pars$Prop_MC, IV_MC=0, I_MC=NM_0*pars$Prop_MC*pars$Prop_IMC, T1_MC=0,T1S_MC=0,T2_MC=0,T2S_MC=0, A_MC=0,
           S_F2=NF_0*pars$Prop_F2, IV_F2=0 ,I_F2=NF_0*pars$Prop_F2*pars$Prop_IF2,T1_F2=0,T1S_F2=0,T2_F2=0,T2S_F2=0, A_F2 = 0,
           S_M2=NM_0*pars$Prop_M2, IV_M2=0, I_M2=NM_0*pars$Prop_IM2*pars$Prop_M2, T1_M2=0,T1S_M2=0,T2_M2=0,T2S_M2=0, A_M2 = 0,D=0, D_HIV = 0)

out <- ode(y = state,times = time,func = hiv_fsw,parms = pars)
plot(out)