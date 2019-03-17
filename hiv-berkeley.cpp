/* ################################################################################
#
#   Translation of the basic model (by Partow Imani) into C++
#   Sean Wu
#   March 2019
#
################################################################################ */

/* libraries */
#include <Rcpp.h>

/* get things out of Rcpp::List */
#define RLIST(name) Rcpp::as<double>(pars[name]);

/* model */
Rcpp::List hiv_fsw(const double t, const Rcpp::NumericVector& state, const Rcpp::List& pars){

  /* out */
  Rcpp::List out(10);

  /* dx/dt */
  Rcpp::NumericVector dx(state.size());
  out[0] = dx;

  /* time dependent parameters */
  double sigma; /* rate infecteds get 1st line treatment */

  double con_FG; /* P(condom use | MG encounters with FG) */
  double con_MG; /* P(condom use | FG encounters with MG) */
  double con_MC; /* P(condom use | FSW encounters with MC) */
  double con_FSW; /* P(condom use | MC encounters with FSW) */
  double con_F2; /* P(condom use | M2 encounters with F2) */
  double con_M2; /* P(condom use | F2 encounters with M2) */

  double circum; /* prop. circumcision */

  double mu_A; /* additional mortality due to AIDS */
  double mu_T; /* additional mortality when on treatment */
  double mu_I; /* additional mortality in chronic infection period */

  double sup1rate; /* rate at which individuals become suppressed on T1 within 6 months of starting */
  double falloff1; /* rate at which individuals fall off treatment on 1st line */

  if(t < 2015.0){

    sigma = RLIST("sigma");

    con_FG = RLIST("con_FG.1");
    con_MG = RLIST("con_MG.1");
    con_MC = RLIST("con_MC.1");
    con_FSW = RLIST("con_FSW.1");
    con_F2 = RLIST("con_F2.1");
    con_M2 = RLIST("con_M2.1");
    circum = RLIST("circum.1");
    mu_A = RLIST("mu_A.1");
    mu_T = RLIST("mu_T.1");
    mu_I = RLIST("mu_I.1");
  } else if(t >= 2015.0){

    sigma = RLIST(pars["sigma.2"]);

    con_FG = RLIST("con_FG.2");
    con_MG = RLIST("con_MG.2");
    con_MC = RLIST("con_MC.2");
    con_FSW = RLIST("con_FSW.2");
    con_F2 = RLIST("con_F2.2");
    con_M2 = RLIST("con_M2.2");
    circum = RLIST("circum.2");
    mu_A = RLIST("mu_A.2");
    mu_T = RLIST("mu_T.2");
    mu_I = RLIST("mu_I.2");
  }

  if(t < 2013.0){
    sup1rate = RLIST("sup.1.rate.pre2013");
    falloff1 = RLIST("fall.off.1.pre2013");
  } else if(t >= 2013.0 && t < 2015.0){
    sup1rate = RLIST("sup.1.rate.pre2015");
    falloff1 = RLIST("fall.off.1.pre2015");
  } else if(t >= 2015.0){
    sup1rate = RLIST("sup.1.rate.post2015");
    falloff1 = RLIST("fall.off.1.post2015");
  }

  /* constants */
  double Tau = RLIST("Tau");
  double mu = RLIST("mu");
  double omega = RLIST("omega");
  double falloff2 = RLIST("fall.off.2");
  double sup2rate = RLIST("sup.2.rate");
  double unsup2rate = RLIST("unsup.2.rate");
  double unsup1rate = RLIST("unsup.1.rate");
  double perc2ndline = RLIST("perc.2nd.line");

  /* general female */
  double S_FG = state[0];
  double IV_FG = state[1];
  double I_FG = state[2];
  double T1_FG = state[3];
  double T1S_FG = state[4];
  double T2_FG = state[5];
  double T2S_FG = state[6];
  double A_FG = state[7];

  /* general male */
  double S_MG = state[8];
  double IV_MG = state[9];
  double I_MG = state[10];
  double T1_MG = state[11];
  double T1S_MG = state[12];
  double T2_MG = state[13];
  double T2S_MG = state[14];
  double A_MG = state[15];

  /* female sex workers */
  double S_FSW = state[16];
  double IV_FSW = state[17];
  double I_FSW = state[18];
  double T1_FSW = state[19];
  double T1S_FSW = state[20];
  double T2_FSW = state[21];
  double T2S_FSW = state[22];
  double A_FSW = state[23];

  /* male clients */
  double S_MC = state[24];
  double IV_MC = state[25];
  double I_MC = state[26];
  double T1_MC = state[27];
  double T1S_MC = state[28];
  double T2_MC = state[29];
  double T2S_MC = state[30];
  double A_MC = state[31];

  /* female 2+ */
  double S_F2 = state[32];
  double IV_F2 = state[33];
  double I_F2 = state[34];
  double T1_F2 = state[35];
  double T1S_F2 = state[36];
  double T2_F2 = state[37];
  double T2S_F2 = state[38];
  double A_F2 = state[39];

  /* male 2+ */
  double S_M2 = state[40];
  double IV_M2 = state[41];
  double I_M2 = state[42];
  double T1_M2 = state[43];
  double T1S_M2 = state[44];
  double T2_M2 = state[45];
  double T2S_M2 = state[46];
  double A_M2 = state[47];

  /* what are these */
  double D = state[48];
  double D_HIV = state[49];

  /* aggregated pop sizes */
  double N_FG = S_FG+IV_FG+I_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG;
  double N_MG = S_MG+IV_MG+I_MG+T1_MG+T2_FG+T1S_MG+T2S_MG+A_MG;

  double N_FSW = S_FSW+IV_FSW+I_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW;
  double N_MC = S_MC+IV_MC+I_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC;

  doube N_F2 = S_F2+IV_F2+I_F2+T1_F2+T2_F2+T1S_F2+T2S_F2+A_F2;
  doube N_M2 = S_M2+IV_M2+I_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2;

  double N_F = N_FG + N_FSW + N_F2;
  double N_M = N_MG + N_MC + N_M2;

  /* contact constraints */
  double C_FGMG = RLIST("C_FGMG")
  double C_MGFG = C_FGMG * N_FG/N_MG;
  double C_FSWC = RLIST("C_FSWC");
  double C_CFSW = C_FSWC * N_FSW/N_MC;
  double C_F2M2 = RLIST("C_F2M2");
  double C_M2F2 = C_F2M2 * N_F2/N_M2;

  /* transmission rates */
  double P_transmission = RLIST("P_transmission");
  double HighV_factor = RLIST("HighV_factor");
  double T_factor = RLIST("T_factor");

  double beta_FGI = P_transmission*C_FGMG;
  double beta_FGV = P_transmission*HighV_factor*C_FGMG;
  double beta_FGT = P_transmission*C_FGMG;
  double beta_FGTS = P_transmission*T_factor*C_FGMG;
  double beta_FGA = P_transmission*C_FGMG;
  double beta_MGI = P_transmission*C_MGFG;
  double beta_MGV = P_transmission*HighV_factor*C_MGFG;
  double beta_MGT = P_transmission*C_MGFG;
  double beta_MGTS = P_transmission*T_factor*C_MGFG;
  double beta_MGA = P_transmission*C_MGFG;

  /* FSW and C */
  double beta_FSWI = P_transmission*C_FSWC;
  double beta_FSWV = P_transmission*HighV_factor*C_FSWC;
  double beta_FSWT = P_transmission*C_FSWC;
  double beta_FSWTS = P_transmission*T_factor*C_FSWC;
  double beta_FSWA = P_transmission*C_FSWC;
  double beta_MCI = P_transmission*C_CFSW;
  double beta_MCV = P_transmission*HighV_factor*C_CFSW;
  double beta_MCT = P_transmission*C_CFSW;
  double beta_MCTS = P_transmission*T_factor*C_CFSW;
  double beta_MCA = P_transmission*C_CFSW;

  /* F2 and M2 */
  double beta_F2I = P_transmission*C_F2M2;
  double beta_F2V = P_transmission*HighV_factor*C_F2M2;
  double beta_F2T = P_transmission*C_F2M2;
  double beta_F2TS = P_transmission*T_factor*C_F2M2;
  double beta_F2A = P_transmission*C_F2M2;
  double beta_M2I = P_transmission*C_M2F2;
  double beta_M2V = P_transmission*HighV_factor*C_M2F2;
  double beta_M2T = P_transmission*C_M2F2;
  double beta_M2TS = P_transmission*T_factor*C_M2F2;
  double beta_M2A = P_transmission*C_M2F2;

  /* force of infection */
  double lambda_MG;

  double lambda_FG;

  /* FSW and C */
  double lambda_MC;

  double lambda_FSW;

  /* F2 and M2 */
  double lambda_M2;

  double lambda_F2;

  /* general females */
  dx[0] = (birth * N_FG) - ((lambda_FG + mu) * S_FG); /* S_FG */
  dx[1] = (lambda_FG * S_FG) - ((mu + Tau) * IV_FG); /* IV_FG */
  dx[2] = (Tau * IV_FG) + (rho * A_FG) + (falloff1 * T1_FG) + (falloff2 * T2_FG) - ((mu_I + mu + sigma + omega) * I_FG); /* I_FG */
  dx[3] = (sigma * I_FG) + (unsup1rate * T1S_FG) - ((mu + mu_I + perc2ndline + sup1rate) * T1_FG); /* T1_FG */
  dx[4] = (sup1rate * T1_FG) - ((mu + mu_T + unsup1rate) * T1S_FG); /* T1S_FG */
  dx[5] = (perc2ndline * T1_FG) - ((mu + mu_T + sup2rate + falloff2) * T2_FG); /* T2_FG */
  dx[6] = (sup2rate * T2_FG) - ((mu + mu_T + unsup2rate) * T2S_FG); /* T2S_FG */
  dx[7] = (omega * I_FG) - ((mu + mu_A + rho) * A_FG) /* A_FG */

  /* general males */

  /* female sex workers */

  /* male clients */

  /* female 2+ */

  /* male 2+ */
};
