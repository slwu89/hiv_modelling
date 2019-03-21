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

sse <- function(x,data_prev,data_pop){
  x
}


### with ten parameter sets
#v     "P_transmission*, HighV_factor*, T_factor*, C_FGMG, C_FSWC, C_F2M2, Prop_FSW, Prop_MC, Prop_IFG, Prop_IMG, Prop_IFSW, Prop_IMC, Prop_F2, Prop_M2, Prop_IF2, Prop_IM2,  Tau,     mu,     mu_I.2,  mu_T.2,  mu_A.2, mu_I.1,  mu_T.1,  mu_A.1, omega, rho, con_F2.2, con_FSW.2, con_FG.2, con_MG.2, con_M2.2, con_MC.2,con_F2.1, con_FSW.1, con_FG.1, con_MG.1, con_M2.1, con_MC.1,con_eff*,circum.2 , circum.1, circum_eff*, sup.2.rate, perc.2nd.line, unsup.1.rate, unsup.2.rate, fall.off.2, birth, sigma.1, sigma.2, sup.1.rate.pre2013, sup.1.rate.pre2015, sup.1.rate.post2015, sup.1.rate.pre2015int, sup.1.rate.post2015.int, fall.1.off.pre2013, fall.1.off.pre2015, fall.1.off.post2015, fall.1.off.pre2015int, fall.1.off.post2015int, phiM, phiF"
minNew = c(0.006,           5.5,          0,            0.5,   70   ,  2     ,    0.06,   0.07,      0.04 ,    0.03  ,    0.18   ,   0.04,    0.1  ,  0.1 ,    0.03,     0.01     ,0.103,  0.0005,   0.03,    0.007,   0.07,   0.03,   0.007,    0.07,   0.05 , 0.05, 0.2  ,    0.5  ,     0.4  ,    0.4  ,    0.15 ,    0.2  ,     0.02,    0.02,      0.02,     0.02,     0.02,      0.02,    0.8  ,  0.4   ,   0.4,         0.5      ,  0.03     ,    0.04      ,   0.02      ,      0.02   ,     0.02  , 0.01,  0,       0.05,       0.46,                0.540,                0.46,                 0.808,                 0.7,                    0.1116,             0.0628,             0.1116,              0.0028,                 0.0438, 0, 0)
maxNew = c(0.06,            24,           0.16,         1.5,   150  ,  7     ,    0.1,     0.1,      0.06 ,    0.05  ,    0.35   ,   0.1,     0.15 , 0.15,     0.1,      0.2     ,0.381,    0.0015,   0.07,   0.015,   0.18,   0.15,   0.1,      0.4,    0.15 , 0.15, 0.4 ,    0.75  ,     0.65 ,    0.65 ,    0.4  ,    0.6  ,     0.5,     0.5,       0.5,      0.5,      0.5,       0.5,     0.95 ,  0.95  ,   0.85,        0.6      ,  0.15     ,    0.15      ,   0.1       ,      0.1    ,     0.1   , 0.03,  0.5,     0.9 ,       0.664,               0.718,                0.718,                0.884,                 0.792,                  0.2587,             0.183,              0.2587,              0.0268,                 0.1, 1, 1)

parRangesNew = cbind(minNew, maxNew)
rownames(parRangesNew) = c("P_transmission", "HighV_factor", "T_factor", "C_FGMG", "C_FSWC", "C_F2M2", "Prop_FSW", "Prop_MC", "Prop_IFG", "Prop_IMG", "Prop_IFSW", "Prop_IMC", "Prop_F2", "Prop_M2", "Prop_IF2", "Prop_IM2",  "Tau",     "mu",     "mu_I.2",  "mu_T.2",  "mu_A.2", "mu_I.1",  "mu_T.1",  "mu_A.1", "omega", "rho", "con_F2.2", "con_FSW.2", "con_FG.2", "con_MG.2", "con_M2.2", "con_MC.2","con_F2.1", "con_FSW.1", "con_FG.1", "con_MG.1", "con_M2.1", "con_MC.1","con_eff" ,"circum.2" , "circum.1", "circum_eff", "sup.2.rate", "perc.2nd.line", "unsup.1.rate", "unsup.2.rate", "fall.off.2", "birth", "sigma.1", "sigma.2", "sup.1.rate.pre2013", "sup.1.rate.pre2015", "sup.1.rate.post2015", "sup.1.rate.pre2015int", "sup.1.rate.post2015int", "fall.off.1.pre2013", "fall.off.1.pre2015", "fall.off.1.post2015", "fall.off.1.pre2015int", "fall.off.1.post2015int", "phiM", "phiF")
