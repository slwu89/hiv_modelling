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

library(Rcpp)
library(lhs)
library(deSolve)

library(foreach)
library(parallel)
library(doParallel)

sourceCpp("hiv-berkeley.cpp")

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
prop_FG <- (1-(pars$Prop_FSW+pars$Prop_F2))
prop_MG <- (1-(pars$Prop_M2+pars$Prop_MC))
state <- c(S_FG=NF_0*prop_FG*(1-pars$Prop_IFG), IV_FG=0 ,I_FG=NF_0*prop_FG*pars$Prop_IFG,T1_FG=0, T1S_FG=0,T2_FG=0, T2S_FG=0,A_FG = 0,
           S_MG=NM_0*prop_MG*(1-pars$Prop_IMG), IV_MG=0, I_MG=NM_0*prop_MG*pars$Prop_IMG, T1_MG=0,T1S_MG=0, T2_MG=0,T2S_MG=0,A_MG = 0,
           S_FSW=NF_0*pars$Prop_FSW*(1-pars$Prop_IFSW), IV_FSW=0, I_FSW=NF_0*pars$Prop_FSW*pars$Prop_IFSW, T1_FSW=0,T1S_FSW=0,T2_FSW=0,T2S_FSW=0, A_FSW=0,
           S_MC=NM_0*pars$Prop_MC*(1-pars$Prop_IMC), IV_MC=0, I_MC=NM_0*pars$Prop_MC*pars$Prop_IMC, T1_MC=0,T1S_MC=0,T2_MC=0,T2S_MC=0, A_MC=0,
           S_F2=NF_0*pars$Prop_F2*(1-pars$Prop_IF2), IV_F2=0 ,I_F2=NF_0*pars$Prop_F2*pars$Prop_IF2,T1_F2=0,T1S_F2=0,T2_F2=0,T2S_F2=0, A_F2 = 0,
           S_M2=NM_0*pars$Prop_M2*(1-pars$Prop_IM2), IV_M2=0, I_M2=NM_0*pars$Prop_IM2*pars$Prop_M2, T1_M2=0,T1S_M2=0,T2_M2=0,T2S_M2=0, A_M2 = 0,D=0, D_HIV = 0)

data_prev <- data.frame(
  # male = c(5.1, 5.9, 6.5, 7.1, 7.5, 7.9, 8.1, 8.1, 8.1, 7.9, 7.7, 7.4, 7.1, 6.7,
  #          6.3, 6, 5.6, 5.3, 5.1, 4.8, 4.6, 4.4, 4.2, 4.1, 3.9, 3.7, 3.6, 3.1),
  # female = c(5.7, 6.7, 7.7, 8.6, 9.3, 9.9, 10.3, 10.5, 10.5, 10.4, 10.2, 9.9, 9.5,
  #            9.1, 8.6, 8.2, 7.7, 7.3, 7.1, 6.8, 6.6, 6.4, 6.2, 6.1, 6, 5.9, 5.8, 6.4),
  male_min = c(4.4, 5.1, 5.7, 6.2, 6.6, 6.9, 7, 7.1, 7, 6.9, 6.7, 6.5, 6.1, 5.7, 5.3, 5, 4.7, 4.4, 4.2, 4, 3.9, 3.7, 3.6, 3.5, 3.3, 3.2, 3.1, 2.7),
  male_max = c(5.9, 6.8, 7.5, 8.1, 8.7, 9, 9.2, 9.3, 9.2, 8.9, 8.6, 8.3, 7.9, 7.5, 7.1, 6.7, 6.4, 6, 5.8, 5.5, 5.3, 5.1, 4.8, 4.6, 4.5, 4.3, 4.1, 3.6),
  female_min = c(5, 5.8, 6.7, 7.5, 8.1, 8.6, 8.9, 9.1, 9.2, 9.1, 8.9, 8.6, 8.2, 7.7, 7.3, 6.8, 6.4, 6.1, 5.9, 5.7, 5.5, 5.4, 5.3, 5.2, 5.1, 5, 4.9, 5.8),
  female_max = c(6.6, 7.8, 8.9, 9.8, 10.7, 11.3, 11.7, 11.9, 12, 11.8, 11.4, 11.1, 10.7, 10.2, 9.7, 9.2, 8.8, 8.4, 8.1, 7.8, 7.5, 7.3, 7.1, 7, 6.9, 6.7, 6.6, 6.9),
  year = 1990:2017
)
data_prev$male <- (data_prev$male_min + data_prev$male_max)/2
data_prev$female <- (data_prev$female_min + data_prev$female_max)/2

# calc SDs such that 99% of prob is within the min and max
sd_fn <- function(sd,min,max,mean){
  (integrate(f = dnorm,lower = min,upper = max,mean=mean,sd=sd)$value - 0.99)^2
}

for(i in 1:nrow(data_prev)){
  data_prev$sd_m[i] <- optimize(f = sd_fn,interval = c(0,10),min=data_prev$male_min[i],max=data_prev$male_max[i],mean=data_prev$male[i])$minimum
  data_prev$sd_f[i] <- optimize(f = sd_fn,interval = c(0,10),min=data_prev$female_min[i],max=data_prev$female_max[i],mean=data_prev$female[i])$minimum
}


data_pop <- data.frame(
  male = c(12594866,14818600,16914713,19441222,22749655,26627475,27474896,28339804,29232512),
  female = c(12864738,15142176,17263329,19969323,23348936,27252482,28097305,28970215,29858880),
  year = c(1990, 1995, 2000, 2005, 2010, 2015, 2016, 2017, 2018)
)

# parameters to fit
par_ranges <- matrix(data =
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
                       # 0,0.01 # sd_p
                       # 0,1e4 # sd_pop
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
                           # "sd_p"
                       ),
                       c("min","max")
                     ))

pnames <- rownames(par_ranges)
snames <- c("S_FG","IV_FG","I_FG","T1_FG","T1S_FG","T2_FG","T2S_FG",
  "A_FG","S_MG","IV_MG","I_MG","T1_MG","T1S_MG","T2_MG",
  "T2S_MG","A_MG","S_FSW","IV_FSW","I_FSW","T1_FSW","T1S_FSW",
  "T2_FSW","T2S_FSW","A_FSW","S_MC","IV_MC","I_MC","T1_MC",
  "T1S_MC","T2_MC","T2S_MC","A_MC","S_F2","IV_F2","I_F2",
  "T1_F2","T1S_F2","T2_F2","T2S_F2","A_F2","S_M2","IV_M2",
  "I_M2","T1_M2","T1S_M2","T2_M2","T2S_M2","A_M2"
)

# useful parameter transformations
expit <- function(x){
  1 / (1 + exp(-x))
}

expit_ab <- function(x,a,b){
  a + ((b-a)*expit(x))
}

logit <- function(x){
  log(x / (1-x))
}

logit_ab <- function(x,a,b){
  logit((x-a)/(b-a))
}

# LHS sampling for optimally uniform points over transformed hypercube
nstart <- 1e3
p <- nrow(par_ranges)
# LHS samples on unit hypercube
starts <- as.list(data.frame(t(lhs::randomLHS(n = nstart,k = p))))
starts <- lapply(starts,function(x){
  # transform to (a,b)
  x <- setNames(qunif(p = x,min = par_ranges[,"min"],max = par_ranges[,"max"]),pnames)
  # transform to (-inf,inf)
  logit_ab(x = x,a = par_ranges[,"min"],b = par_ranges[,"max"])
})

# x begins at one of the "starts"
obj <- function(x,ranges){

  # parameters
  theta <- x[1:40] # parameters

  init_state <- setNames(initstate(x = x,ranges = ranges),snames)

  sim <- deSolve::ode(y = init_state,times = 1990:2020,func = hiv_fsw_fit,parms = theta,ranges = ranges)

  # minimize negative log-likelihood
  negll <- sum(dnorm(x = data_prev$male, mean = sim[1:28,"Prev_M"],
                 sd = data_prev$sd_m*0.01,
                 log = TRUE))
  negll <- negll + sum(dnorm(x = data_prev$female, mean = sim[1:28,"Prev_F"],
                             sd = data_prev$sd_f*0.01,
                             log = TRUE))

  negll <- negll + sum(dpois(x = data_pop$female, lambda = sim[sim[,"time"] %in% data_pop$year,"N_F"],log = TRUE))
  negll <- negll + sum(dpois(x = data_pop$male, lambda = sim[sim[,"time"] %in% data_pop$year,"N_M"],log = TRUE))

  return(-negll)
}


cl <- parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)

opt <- foreach(it = starts,.combine = "c",.inorder = FALSE, .packages = c("deSolve","stats")) %do% {
    optim(par = it,fn = obj,
          method = "Nelder-Mead", control=list(trace=0,maxit=5e4,reltol=1e-8),
          ranges = par_ranges)
}

parallel::stopCluster(cl)
rm(cl);gc()