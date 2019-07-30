################################################################################
#       ___    ____               __  ___          __     __
#      /   |  / __/_  ______ _   /  |/  /___  ____/ /__  / /
#     / /| | / /_/ / / / __ `/  / /|_/ / __ \/ __  / _ \/ /
#    / ___ |/ __/ /_/ / /_/ /  / /  / / /_/ / /_/ /  __/ /
#   /_/  |_/_/  \__, /\__,_/  /_/  /_/\____/\__,_/\___/_/
#              /____/
#
#   Partow Imani & Sean Wu
#   August 2019
#
################################################################################


rm(list=ls());gc()

library(here)
library(Rcpp)
library(deSolve)

library(lhs)
library(pbmcapply)
library(ggplot2)
library(gridExtra)

sourceCpp(here("src/hiv-berkeley-int.cpp"))

# grab the MLE fits and take the best
mle <- readRDS(here("fits.rds"))

opt_negll <- unname(sapply(X = mle,FUN = function(x){x$value},USE.NAMES = FALSE))
opt_vals <- unname(lapply(X = mle,FUN = function(x){x$par}))

# take the 10 best and follow up with BFGS
mle_par <- opt_vals[order(opt_negll,decreasing = F)[1]][[1]]

# transform from (-inf,inf) to (a,b) ... back to natural range of parameter
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

# going from unconstrained to constrained space
expit <- function(x){
  1 / (1 + exp(-x))
}

expit_ab <- function(x,a,b){
  a + ((b-a)*expit(x))
}

# the MLE parameters
mle_par_trans <- mapply(FUN = function(x,name){
  a <- par_ranges[name,1]
  b <- par_ranges[name,2]
  return(expit_ab(x,a,b))
},x = mle_par,name = names(mle_par),SIMPLIFY = TRUE,USE.NAMES = TRUE)

samples_names <- names(mle_par_trans)
samples_names[which(samples_names == "sup1rate_pre2015")] <- "sup1rate_pre2015_int"
samples_names[which(samples_names == "sup1rate_post2015")] <- "sup1rate_post2015_int"
samples_names[which(samples_names == "falloff1_pre2015")] <- "falloff1_pre2015_int"
samples_names[which(samples_names == "falloff1_post2015")] <- "falloff1_post2015_int"

samples_mle <- mle_par_trans
names(samples_mle) <- samples_names

# samples for the intervention parameters
int_ranges <- matrix(data = c(
  0.8080,0.8840, # sup1rate_pre2015_int
  0.7000,0.7920, # sup1rate_post2015_int
  0.0028,0.0268, # falloff1_pre2015_int
  0.0438,0.1000 # falloff1_post2015_int
),nrow = 4,ncol = 2,byrow = TRUE,
dimnames = list(c("sup1rate_pre2015_int",
                  "sup1rate_post2015_int",
                  "falloff1_pre2015_int",
                  "falloff1_post2015_int"),
                c("min","max"))
)

# number of samples
nsamp <- 1e4
samples <- as.list(data.frame(t(lhs::randomLHS(n = nsamp,k = 4))))
samples <- lapply(samples,function(x){
  # transform to uniform density on (a,b)
  x <- setNames(object = qunif(p = x,min = int_ranges[,"min"],max = int_ranges[,"max"]),nm = rownames(int_ranges))
  samples_mle[names(x)] <- x
  return(samples_mle)
})

# initial state
snames <- c("S_FG","IV_FG","I_FG","T1_FG","T1S_FG","T2_FG","T2S_FG",
            "A_FG","S_MG","IV_MG","I_MG","T1_MG","T1S_MG","T2_MG",
            "T2S_MG","A_MG","S_FSW","IV_FSW","I_FSW","T1_FSW","T1S_FSW",
            "T2_FSW","T2S_FSW","A_FSW","S_MC","IV_MC","I_MC","T1_MC",
            "T1S_MC","T2_MC","T2S_MC","A_MC","S_F2","IV_F2","I_F2",
            "T1_F2","T1S_F2","T2_F2","T2S_F2","A_F2","S_M2","IV_M2",
            "I_M2","T1_M2","T1S_M2","T2_M2","T2S_M2","A_M2"
)
init_state <- setNames(object = rep(0,length(snames)),nm = snames)

NF_0 = 12864738
NM_0 = 12594866

prop_FG = (1.0 - (mle_par_trans["Prop_FSW"] + mle_par_trans["Prop_F2"]));
prop_MG = (1.0 - (mle_par_trans["Prop_M2"] + mle_par_trans["Prop_MC"]));

init_state[1] = NF_0 * prop_FG * (1.0 - mle_par_trans["Prop_IFG"]) # S_FG
init_state[3] = NF_0 * prop_FG * mle_par_trans["Prop_IFG"] # I_FG
init_state[9] = NM_0 * prop_MG * (1.0 - mle_par_trans["Prop_IMG"]) # S_MG
init_state[11] = NM_0 * prop_MG * mle_par_trans["Prop_IMG"] # I_MG
init_state[17] = NF_0 * mle_par_trans["Prop_FSW"] * (1.0 - mle_par_trans["Prop_IFSW"]) # S_FSW
init_state[19] = NF_0 * mle_par_trans["Prop_FSW"] * mle_par_trans["Prop_IFSW"] # I_FSW
init_state[25] = NM_0 * mle_par_trans["Prop_MC"] * (1.0 - mle_par_trans["Prop_IMC"]) # S_MC
init_state[27] = NM_0 * mle_par_trans["Prop_MC"] * mle_par_trans["Prop_IMC"] # I_MC
init_state[33] = NF_0 * mle_par_trans["Prop_F2"] * (1.0 - mle_par_trans["Prop_IF2"]) # S_F2
init_state[35] = NF_0 * mle_par_trans["Prop_F2"] * mle_par_trans["Prop_IF2"] # I_F2
init_state[41] = NM_0 * mle_par_trans["Prop_M2"] * (1.0 - mle_par_trans["Prop_IM2"]) # S_M2
init_state[43] = NM_0 * mle_par_trans["Prop_IM2"] * mle_par_trans["Prop_M2"] # I_M2

# run simulations
int_sims <- pbmcapply::pbmclapply(X = samples,FUN = function(samp){
  deSolve::ode(y = init_state,times = 1990:2020,func = hiv_fsw_int,parms = samp)
},mc.cores = 6)
