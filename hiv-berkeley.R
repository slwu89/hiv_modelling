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
library(ggplot2)

sourceCpp(here("src/hiv-berkeley.cpp"))

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


data_prev <- data.frame(
  male = c(5.1, 5.9, 6.5, 7.1, 7.5, 7.9, 8.1, 8.1, 8.1, 7.9, 7.7, 7.4, 7.1, 6.7,
           6.3, 6, 5.6, 5.3, 5.1, 4.8, 4.6, 4.4, 4.2, 4.1, 3.9, 3.7, 3.6, 3.1),
  female = c(5.7, 6.7, 7.7, 8.6, 9.3, 9.9, 10.3, 10.5, 10.5, 10.4, 10.2, 9.9, 9.5,
             9.1, 8.6, 8.2, 7.7, 7.3, 7.1, 6.8, 6.6, 6.4, 6.2, 6.1, 6, 5.9, 5.8, 6.4),
  year = 1990:2017
)

data_pop <- data.frame(
  male = c(12594866,14818600,16914713,19441222,22749655,26627475,27474896,28339804,29232512),
  female = c(12864738,15142176,17263329,19969323,23348936,27252482,28097305,28970215,29858880),
  year = c(1990, 1995, 2000, 2005, 2010, 2015, 2016, 2017, 2018)
)

ggplot(data = data_prev) +
  geom_line(aes(x=year,y=male),col="firebrick3") +
  geom_line(aes(x=year,y=female),col="steelblue") +
  theme_bw() +
  xlab("Year") + ylab("Prevalence (blue = F, red = M)")

ggplot(data = data_pop) +
  geom_line(aes(x=year,y=male),col="firebrick3") +
  geom_line(aes(x=year,y=female),col="steelblue") +
  theme_bw() +
  xlab("Year") + ylab("Population Size (blue = F, red = M)")

# parameters to fit
par_ranges <- matrix(data = ,
                     c(0.0060,0.0600, # P_transmission
                       5.5000,24.0000, # HighV_factor
                       0.0000,0.1600, # T_factor
                       0.5000,1.5000, # C_FGMG
                       70.0000,150.0000, # C_FSWC
                       2.0000,7.0000, # C_F2M2
                       0.1030,0.3810, # Tau
                       0.0005,0.0015, # mu
                       0.0300,0.0700, # mu_I_2
                       0.0070,0.0150, # mu_T_2
                       0.0700,0.1800, # mu_A_2
                       0.0300,0.1500, # mu_I_1
                       0.0070,0.1000, # mu_T_1
                       0.0700,0.4000, # mu_A_1
                       0.0500,0.1500, # omega
                       0.0500,0.1500, # rho
                       0.4000,0.6500, # con_FG_2
                       0.0200,0.5000, # con_FG_1
                       0.5000,0.7500, # con_FSW_2
                       0.0200,0.5000, # con_FSW_1
                       0.2000,0.4000, # con_F2_2
                       0.0200,0.5000, # con_F2_1
                       0.8000,0.9500, # con_eff
                       0.4000,0.9500, # circum_2
                       0.4000,0.8500, # circum_1
                       0.5000,0.6000, # circum_eff
                       0.0300,0.1500, # sup2rate
                       0.0400,0.1500, # perc2ndline
                       0.0200,0.1000, # unsup2rate
                       0.0200,0.1000, # falloff2
                       0.0200,0.1000, # unsup1rate
                       0.4600,0.6640, # sup1rate_pre2013
                       0.1116,0.2587, # falloff1_pre2013
                       0.5400,0.7180, # sup1rate_pre2015
                       0.0628,0.1830, # falloff1_pre2015
                       0.4600,0.7180, # sup1rate_post2015
                       0.1116,0.2587, # falloff1_post2015
                       0.0100,0.0300, # birth
                       0.0500,0.9000, # sigma_2
                       0.0000,0.5000, # sigma_1
                       0.0600,0.1000, # Prop_FSW
                       0.0700,0.1000, # Prop_MC
                       0.0400,0.0600, # Prop_IFG
                       0.0300,0.0500, # Prop_IMG
                       0.1800,0.3500, # Prop_IFSW
                       0.0400,0.1000, # Prop_IMC
                       0.1000,0.1500, # Prop_F2
                       0.1000,0.1500, # Prop_M2
                       0.0300,0.1000, # Prop_IF2
                       0.0100,0.2000 # Prop_IM2
                       ),
                     nrow = 50,ncol = 2,byrow = TRUE,dimnames = list(
                       c(
                         "P_transmission",
                           "HighV_factor",
                           "T_factor",
                           "C_FGMG",
                           "C_FSWC",
                           "C_F2M2",
                           "Tau",
                           "mu",
                           "mu_I_2",
                           "mu_T_2",
                           "mu_A_2",
                           "mu_I_1",
                           "mu_T_1",
                           "mu_A_1",
                           "omega",
                           "rho",
                           "con_FG_2",
                           "con_FG_1",
                           "con_FSW_2",
                           "con_FSW_1",
                           "con_F2_2",
                           "con_F2_1",
                           "con_eff",
                           "circum_2",
                           "circum_1",
                           "circum_eff",
                           "sup2rate",
                           "perc2ndline",
                           "unsup2rate",
                           "falloff2",
                           "unsup1rate",
                           "sup1rate_pre2013",
                           "falloff1_pre2013",
                           "sup1rate_pre2015",
                           "falloff1_pre2015",
                           "sup1rate_post2015",
                           "falloff1_post2015",
                           "birth",
                           "sigma_2",
                           "sigma_1",
                           "Prop_FSW",
                           "Prop_MC",
                           "Prop_IFG",
                           "Prop_IMG",
                           "Prop_IFSW",
                           "Prop_IMC",
                           "Prop_F2",
                           "Prop_M2",
                           "Prop_IF2",
                           "Prop_IM2"
                       ),
                       NULL
                     ))

