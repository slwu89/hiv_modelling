/* ################################################################################
#
#   Translation of the basic model (by Partow Imani) into C++
#   ugly code, optimized for speed
#   Sean Wu
#   March 2019
#
################################################################################ */

/* libraries */
#include <Rcpp.h>
#include <math.h>
#include <unordered_map>
#include <string>

// [[Rcpp::plugins(cpp14)]]


/* ################################################################################
#   Helper functions
################################################################################ */

/* get things out of Rcpp::List */
#define RLIST(name) Rcpp::as<double>(pars[name]);

/* get things out of Rcpp::NumericVector */
#define RVEC(name) Rcpp::as<double>(pars[name]);

/* names for output */
static const Rcpp::CharacterVector outnames = Rcpp::CharacterVector::create("dx","N_F","N_M","Prev_F","Prev_M","Prev_F");

/* if x in (0,1) logit(x) -> (-inf,inf) */
inline double logit(const double x){
  return log(x / (1.0 - x));
};

/* if x in (a,b) logit_ab(x) -> (-inf,inf) */
inline double logit_ab(const double x, const double a, const double b){
  return logit((x-a)/(b-a));
};

/* if x in (-inf,inf) expit(x) -> (0,1) */
inline double expit(const double x){
  return 1.0 / (1.0 + exp(-x));
};

/* if x in (-inf,inf) expit_ab(x) -> (a,b) */
inline double expit_ab(const double x, const double a, const double b){
  return a + ((b-a)*expit(x));
};

/* lookup table for named parameters (getting a,b bounds from par ranges): returns what row it's in */
static const std::unordered_map<std::string,size_t> par_ranges = {
  {"P_transmission",0},
  {"HighV_factor",1},
  {"T_factor",2},
  {"C_FGMG",3},
  {"C_FSWC",4},
  {"C_F2M2",5},
  {"Tau",6},
  {"mu",7},
  {"mu_I_2",8},
  {"mu_T_2",9},
  {"mu_A_2",10},
  {"mu_I_1",11},
  {"mu_T_1",12},
  {"mu_A_1",13},
  {"omega",14},
  {"rho",15},
  {"con_FG_2",16},
  {"con_FG_1",17},
  {"con_FSW_2",18},
  {"con_FSW_1",19},
  {"con_F2_2",20},
  {"con_F2_1",21},
  {"con_eff",22},
  {"circum_2",23},
  {"circum_1",24},
  {"circum_eff",25},
  {"sup2rate",26},
  {"perc2ndline",27},
  {"unsup2rate",28},
  {"falloff2",29},
  {"unsup1rate",30},
  {"sup1rate_pre2013",31},
  {"falloff1_pre2013",32},
  {"sup1rate_pre2015",33},
  {"falloff1_pre2015",34},
  {"sup1rate_post2015",35},
  {"falloff1_post2015",36},
  {"birth",37},
  {"sigma_2",38},
  {"sigma_1",39},
  {"Prop_FSW",40},
  {"Prop_MC",41},
  {"Prop_IFG",42},
  {"Prop_IMG",43},
  {"Prop_IFSW",44},
  {"Prop_IMC",45},
  {"Prop_F2",46},
  {"Prop_M2",47},
  {"Prop_IF2",48},
  {"Prop_IM2",49}
};

/* helper function to extract parameters and transform from (-Inf,Inf) -> (a,b) */
inline double get_par(const std::string par, const Rcpp::NumericVector& pars, const Rcpp::NumericMatrix& ranges){
  double par_a = ranges(par_ranges.at(par),0);
  double par_b = ranges(par_ranges.at(par),1);
  return expit_ab(Rcpp::as<double>(pars[par]),par_a,par_b);
};


/* ################################################################################
#   state variables
################################################################################ */

// [[Rcpp::export]]
Rcpp::NumericVector initstate(const Rcpp::NumericVector& x, const Rcpp::NumericMatrix& ranges){

  double Prop_FSW = get_par("Prop_FSW",x,ranges);
  double Prop_MC = get_par("Prop_MC",x,ranges);
  double Prop_IFG = get_par("Prop_IFG",x,ranges);
  double Prop_IMG = get_par("Prop_IMG",x,ranges);
  double Prop_IFSW = get_par("Prop_IFSW",x,ranges);
  double Prop_IMC = get_par("Prop_IMC",x,ranges);
  double Prop_F2 = get_par("Prop_F2",x,ranges);
  double Prop_M2 = get_par("Prop_M2",x,ranges);
  double Prop_IF2 = get_par("Prop_IF2",x,ranges);
  double Prop_IM2 = get_par("Prop_IM2",x,ranges);

  static const double NF_0 = 12864738;
  static const double NM_0 = 12594866;

  double prop_FG = (1.0 - (Prop_FSW + Prop_F2));
  double prop_MG = (1.0 - (Prop_M2 + Prop_MC));

  Rcpp::NumericVector state(48);

  state[0] = NF_0 * prop_FG * (1.0 - Prop_IFG); /* S_FG */
  state[2] = NF_0 * prop_FG * Prop_IFG; /* I_FG */
  state[8] = NM_0 * prop_MG * (1.0 - Prop_IMG); /* S_MG */
  state[10] = NM_0 * prop_MG * Prop_IMG; /* I_MG */
  state[16] = NF_0 * Prop_FSW * (1.0 - Prop_IFSW); /* S_FSW */
  state[18] = NF_0 * Prop_FSW * Prop_IFSW; /* I_FSW */
  state[24] = NM_0 * Prop_MC * (1.0 - Prop_IMC); /* S_MC */
  state[26] = NM_0 * Prop_MC * Prop_IMC; /* I_MC */
  state[32] = NF_0 * Prop_F2 * (1.0 - Prop_IF2); /* S_F2 */
  state[34] = NF_0 * Prop_F2 * Prop_IF2; /* I_F2 */
  state[40] = NM_0 * Prop_M2 * (1.0 - Prop_IM2); /* S_M2 */
  state[42] = NM_0 * Prop_IM2 * Prop_M2; /* I_M2 */

  return state;
};


/* ################################################################################
#   model for fitting:
#   all parameters are optimized on an unconstrained space
#   and then transformed back to the (a,b) interval they lie within
################################################################################ */

// [[Rcpp::export]]
Rcpp::List hiv_fsw_fit(const double t, const Rcpp::NumericVector& state, const Rcpp::NumericVector& pars, const Rcpp::NumericMatrix& ranges){

  /* out */
  Rcpp::List out(6);

  /* dx/dt */
  Rcpp::NumericVector dx(state.size());
  out[0] = dx;

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

  /* aggregated pop sizes */
  double N_FG = S_FG+IV_FG+I_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG;
  double N_MG = S_MG+IV_MG+I_MG+T1_MG+T2_FG+T1S_MG+T2S_MG+A_MG;

  double N_FSW = S_FSW+IV_FSW+I_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW;
  double N_MC = S_MC+IV_MC+I_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC;

  double N_F2 = S_F2+IV_F2+I_F2+T1_F2+T2_F2+T1S_F2+T2S_F2+A_F2;
  double N_M2 = S_M2+IV_M2+I_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2;

  double N_F = N_FG + N_FSW + N_F2;
  double N_M = N_MG + N_MC + N_M2;

  /* contact constants */

  /* general pop */
  double C_FGMG = get_par("C_FGMG",pars,ranges);
  double C_MGFG = C_FGMG * N_FG/N_MG;

  /* female sex workers & clients */
  double C_FSWC = get_par("C_FSWC",pars,ranges);
  double C_CFSW = C_FSWC * N_FSW/N_MC;

  /* 2+ pop */
  double C_F2M2 = get_par("C_F2M2",pars,ranges);
  double C_M2F2 = C_F2M2 * N_F2/N_M2;

  /* transmission parameters */
  double P_transmission = get_par("P_transmission",pars,ranges);
  double HighV_factor = get_par("HighV_factor",pars,ranges);
  double T_factor = get_par("T_factor",pars,ranges);

  /* other constants */
  double Tau = get_par("Tau",pars,ranges);
  double mu = get_par("mu",pars,ranges);
  double omega = get_par("omega",pars,ranges);
  double falloff2 = get_par("falloff2",pars,ranges);
  double sup2rate = get_par("sup2rate",pars,ranges);
  double unsup2rate = get_par("unsup2rate",pars,ranges);
  double unsup1rate = get_par("unsup1rate",pars,ranges);
  double perc2ndline = get_par("perc2ndline",pars,ranges);

  double con_eff = get_par("con_eff",pars,ranges);
  double circum_eff = get_par("circum_eff",pars,ranges);

  double birth = get_par("birth",pars,ranges);
  double rho = get_par("rho",pars,ranges);

  /* time dependent parameters */
  double sigma; /* rate infecteds get 1st line treatment */

  double con_FG; /* P(condom use | MG encounters with FG) */
  double con_FSW; /* P(condom use | MC encounters with FSW) */
  double con_F2; /* P(condom use | M2 encounters with F2) */

  double circum; /* prop. circumcision */

  double mu_A; /* additional mortality due to AIDS */
  double mu_T; /* additional mortality when on treatment */
  double mu_I; /* additional mortality in chronic infection period */

  double sup1rate; /* rate at which individuals become suppressed on T1 within 6 months of starting */
  double falloff1; /* rate at which individuals fall off treatment on 1st line */

  if(t < 2015.0){

    sigma = get_par("sigma_1",pars,ranges);

    con_FG = get_par("con_FG_1",pars,ranges);
    con_FSW = get_par("con_FSW_1",pars,ranges);
    con_F2 = get_par("con_F2_1",pars,ranges);
    circum = get_par("circum_1",pars,ranges);
    mu_A = get_par("mu_A_1",pars,ranges);
    mu_T = get_par("mu_T_1",pars,ranges);
    mu_I = get_par("mu_I_1",pars,ranges);
  } else if(t >= 2015.0){

    sigma = get_par("sigma_2",pars,ranges);

    con_FG = get_par("con_FG_2",pars,ranges);
    con_FSW = get_par("con_FSW_2",pars,ranges);
    con_F2 = get_par("con_F2_2",pars,ranges);
    circum = get_par("circum_2",pars,ranges);
    mu_A = get_par("mu_A_2",pars,ranges);
    mu_T = get_par("mu_T_2",pars,ranges);
    mu_I = get_par("mu_I_2",pars,ranges);
  } else {

    sigma = NAN;

    con_FG = NAN;
    con_FSW = NAN;
    con_F2 = NAN;
    circum = NAN;
    mu_A = NAN;
    mu_T = NAN;
    mu_I = NAN;
  }

  if(t < 2013.0){
    sup1rate = get_par("sup1rate_pre2013",pars,ranges);
    falloff1 = get_par("falloff1_pre2013",pars,ranges);
  } else if(t >= 2013.0 && t < 2015.0){
    sup1rate = get_par("sup1rate_pre2015",pars,ranges);
    falloff1 = get_par("falloff1_pre2015",pars,ranges);
  } else if(t >= 2015.0){
    sup1rate = get_par("sup1rate_post2015",pars,ranges);
    falloff1 = get_par("falloff1_post2015",pars,ranges);
  } else {
    sup1rate = NAN; /* to shut up the compiler... */
    falloff1 = NAN;
  }

  /* these are death rates (all cause and HIV specific mortality) */
  // double D = state[48];
  // double D_HIV = state[49];

  /* general population */
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

  /* forces of infection */

  /* general population */
  double con_circ_G = ((1.0-con_FG)+con_FG*(1.0-con_eff))*((1.0-circum) + circum*(1.0-circum_eff));
  double lambda_MG = (beta_MGI * (I_FG/N_FG) * con_circ_G) + /* FOI from I females  */
    (beta_MGV * (IV_FG/N_FG) * con_circ_G) + /* FOI from IV females */
    (beta_MGT * (T1_FG/N_FG) * con_circ_G) + /* FOI from T1 females */
    (beta_MGTS * (T1S_FG/N_FG) * con_circ_G) + /* FOI from T1 virally suppressed females */
    (beta_MGT * (T2_FG/N_FG) * con_circ_G) + /* FOI from T2 females */
    (beta_MGTS * (T2S_FG/N_FG) * con_circ_G) + /* FOI from T2 virally suppressed females */
    (beta_MGA * (A_FG/N_FG) * con_circ_G); /* FOI from AIDS females */

  double lambda_FG = (beta_FGI * (I_MG/N_MG) * con_circ_G) +
    (beta_FGV * (IV_MG/N_MG) * con_circ_G) +
    (beta_FGT * (T1_MG/N_MG) * con_circ_G) +
    (beta_FGTS * (T1S_MG/N_MG) * con_circ_G) +
    (beta_FGT * (T2_MG/N_MG) * con_circ_G) +
    (beta_FGTS * (T2S_MG/N_MG) * con_circ_G) +
    (beta_FGA * (A_MG/N_MG) * con_circ_G);

  /* FSW and C */
  double con_circ_FSW = ((1.0-con_FSW)+con_FSW*(1.0-con_eff))*((1.0-circum) + circum*(1.0-circum_eff));
  double lambda_MC = (beta_MCI*(I_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCV*(IV_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCT*(T1_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCTS*(T1S_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCT*(T2_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCTS*(T2S_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCA*(A_FSW/N_FSW) * con_circ_FSW);

  double lambda_FSW = (beta_FSWI * (I_MC/N_MC) * con_circ_FSW) +
    (beta_FSWV * (IV_MC/N_MC) * con_circ_FSW) +
    (beta_FSWT * (T1_MC/N_MC) * con_circ_FSW) +
    (beta_FSWTS * (T1S_MC/N_MC) * con_circ_FSW) +
    (beta_FSWT * (T2_MC/N_MC) * con_circ_FSW) +
    (beta_FSWTS * (T2S_MC/N_MC) * con_circ_FSW) +
    (beta_FSWA * (A_MC/N_MC) * con_circ_FSW);

  /* F2 and M2 */
  double con_circ_2 = ((1.0-con_F2)+con_F2*(1.0-con_eff))*((1.0-circum) + circum*(1.0-circum_eff));
  double lambda_M2 = (beta_M2I * (I_F2/N_F2) * con_circ_2) +
    (beta_M2V * (IV_F2/N_F2) * con_circ_2) +
    (beta_M2T * (T1_F2/N_F2) * con_circ_2) +
    (beta_M2TS * (T1S_F2/N_F2) * con_circ_2) +
    (beta_M2T * (T2_F2/N_F2) * con_circ_2) +
    (beta_M2TS * (T2S_F2/N_F2) * con_circ_2) +
    (beta_M2A * (A_F2/N_F2) * con_circ_2);

  double lambda_F2 = (beta_F2I * (I_M2/N_M2) * con_circ_2) +
    (beta_F2V * (IV_M2/N_M2) * con_circ_2) +
    (beta_F2T * (T1_M2/N_M2) * con_circ_2) +
    (beta_F2TS * (T1S_M2/N_M2) * con_circ_2) +
    (beta_F2T * (T2_M2/N_M2) * con_circ_2) +
    (beta_F2TS * (T2S_M2/N_M2) * con_circ_2) +
    (beta_F2A * (A_M2/N_M2) * con_circ_2);

  /* general females */
  dx[0] = (birth * N_FG) - ((lambda_FG + mu) * S_FG); /* S_FG */
  dx[1] = (lambda_FG * S_FG) - ((mu + Tau) * IV_FG); /* IV_FG */
  dx[2] = (Tau * IV_FG) + (rho * A_FG) + (falloff1 * T1_FG) + (falloff2 * T2_FG) - ((mu_I + mu + sigma + omega) * I_FG); /* I_FG */
  dx[3] = (sigma * I_FG) + (unsup1rate * T1S_FG) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_FG); /* T1_FG */
  dx[4] = (sup1rate * T1_FG) - ((mu + mu_T + unsup1rate) * T1S_FG); /* T1S_FG */
  dx[5] = (perc2ndline * T1_FG) + (unsup2rate * T2S_FG) - ((mu + mu_T + sup2rate + falloff2) * T2_FG); /* T2_FG */
  dx[6] = (sup2rate * T2_FG) - ((mu + mu_T + unsup2rate) * T2S_FG); /* T2S_FG */
  dx[7] = (omega * I_FG) - ((mu + mu_A + rho) * A_FG); /* A_FG */

  /* general males */
  dx[8] = (birth * N_MG) - ((lambda_MG + mu) * S_MG); /* S_MG */
  dx[9] = (lambda_MG * S_MG) - ((mu + Tau) * IV_MG); /* IV_FG */
  dx[10] = (Tau * IV_MG) + (rho * A_MG) + (falloff1 * T1_MG) + (falloff2 * T2_MG) - ((mu_I + mu + sigma + omega) * I_MG); /* I_MG */
  dx[11] = (sigma * I_MG) + (unsup1rate * T1S_MG) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_MG); /* T1_MG */
  dx[12] = (sup1rate * T1_MG) - ((mu + mu_T + unsup1rate) * T1S_MG); /* T1S_MG */
  dx[13] = (perc2ndline * T1_MG) + (unsup2rate * T2S_MG) - ((mu + mu_T + sup2rate + falloff2) * T2_MG); /* T2_MG */
  dx[14] = (sup2rate * T2_MG) - ((mu + mu_T + unsup2rate) * T2S_MG); /* T2S_MG */
  dx[15] = (omega * I_MG) - ((mu + mu_A + rho) * A_MG); /* A_MG */

  /* female sex workers */
  dx[16] = (birth * N_FSW) - ((lambda_FSW + mu) * S_FSW); /* S_FSW */
  dx[17] = (lambda_FSW * S_FSW) - ((mu + Tau) * IV_FSW); /* IV_FSW */
  dx[18] = (Tau * IV_FSW) + (rho * A_FSW) + (falloff1 * T1_FSW) + (falloff2 * T2_FSW) - ((mu_I + mu + sigma + omega) * I_FSW); /* I_FSW */
  dx[19] = (sigma * I_FSW) + (unsup1rate * T1S_FSW) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_FSW); /* T1_FSW */
  dx[20] = (sup1rate * T1_FSW) - ((mu + mu_T + unsup1rate) * T1_FSW); /* T1S_FSW */
  dx[21] = (perc2ndline * T1_FSW) + (unsup2rate * T2S_FSW) - ((mu + mu_T + sup2rate + falloff2) * T2_FSW); /* T2_FSW */
  dx[22] = (sup2rate * T2_FSW) - ((mu + mu_T + unsup2rate) * T2S_FSW); /* T2S_FSW */
  dx[23] = (omega * I_FSW) - ((mu + mu_A + rho) * A_FSW); /* A_FSW */

  /* male clients */
  dx[24] = (birth * N_MC) - ((lambda_MC + mu) * S_MC); /* S_MC */
  dx[25] = (lambda_MC * S_MC) - ((mu + Tau) * IV_MC); /* IV_MC */
  dx[26] = (Tau * IV_MC) + (rho * A_MC) + (falloff1 * T1_MC) + (falloff2 * T2_MC) - ((mu_I + mu + sigma + omega) * I_MC); /* I_MC */
  dx[27] = (sigma * I_MC) + (unsup1rate * T1S_MC) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_MC); /* T1_MC */
  dx[28] = (sup1rate * T1_MC) - ((mu + mu_T + unsup1rate) * T1_MC); /* T1S_MC */
  dx[29] = (perc2ndline * T1_MC) + (unsup2rate * T2S_MC) - ((mu + mu_T + sup2rate + falloff2) * T2_MC); /* T2_MC */
  dx[30] = (sup2rate * T2_MC) - ((mu + mu_T + unsup2rate) * T2S_MC); /* T2S_MC */
  dx[31] = (omega * I_MC) - ((mu + mu_A + rho) * A_MC); /* A_MC */

  /* female 2+ */
  dx[32] = (birth * N_F2) - ((lambda_F2 + mu) * S_F2); /* S_F2 */
  dx[33] = (lambda_F2 * S_F2) - ((mu + Tau) * IV_F2); /* IV_F2 */
  dx[34] = (Tau * IV_F2) + (rho * A_F2) + (falloff1 * T1_F2) + (falloff2 * T2_F2) - ((mu_I + mu + sigma + omega) * I_F2); /* I_F2 */
  dx[35] = (sigma * I_F2) + (unsup1rate * T1S_F2) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_F2); /* T1_F2 */
  dx[36] = (sup1rate * T1_F2) - ((mu + mu_T + unsup1rate) * T1_F2); /* T1S_F2 */
  dx[37] = (perc2ndline * T1_F2) + (unsup2rate * T2S_F2) - ((mu + mu_T + sup2rate + falloff2) * T2_F2); /* T2_F2 */
  dx[38] = (sup2rate * T2_F2) - ((mu + mu_T + unsup2rate) * T2S_F2); /* T2S_F2 */
  dx[39] = (omega * I_F2) - ((mu + mu_A + rho) * A_F2); /* A_F2 */

  /* male 2+ */
  dx[40] = (birth * N_M2) - ((lambda_M2 + mu) * S_M2); /* S_M2 */
  dx[41] = (lambda_M2 * S_M2) - ((mu + Tau) * IV_M2); /* IV_M2 */
  dx[42] = (Tau * IV_M2) + (rho * A_M2) + (falloff1 * T1_M2) + (falloff2 * T2_M2) - ((mu_I + mu + sigma + omega) * I_M2); /* I_M2 */
  dx[43] = (sigma * I_M2) + (unsup1rate * T1S_M2) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_M2); /* T1_M2 */
  dx[44] = (sup1rate * T1_M2) - ((mu + mu_T + unsup1rate) * T1_M2); /* T1S_M2 */
  dx[45] = (perc2ndline * T1_M2) + (unsup2rate * T2S_M2) - ((mu + mu_T + sup2rate + falloff2) * T2_M2); /* T2_M2 */
  dx[46] = (sup2rate * T2_M2) - ((mu + mu_T + unsup2rate) * T2S_M2); /* T2S_M2 */
  dx[47] = (omega * I_M2) - ((mu + mu_A + rho) * A_M2); /* A_M2 */

  // /* all-cause and HIV-specific mortality */
  // dx[48] = mu*(S_FG + IV_FG + I_FG + T1_FG + T1S_FG+ T2_FG + T2S_FG+ A_FG +
  //            S_MG + IV_MG + I_MG + T1_MG + T1S_MG + T2_MG + T2S_MG + A_MG +
  //            S_FSW + IV_FSW + I_FSW + A_FSW + T1_FSW + T1S_FSW + T2_FSW + T2S_FSW +
  //            S_MC + IV_MC + I_MC + T1_MC + T1S_MC + T2_MC + T2S_MC + A_MC +
  //            S_F2 + IV_F2 + I_F2 + T1_F2 + T1S_F2 + T2_F2 + T2S_F2 + A_F2 +
  //            S_M2 + IV_M2 + I_M2 + T1_M2 + T1S_M2 + T2_M2 + T2S_M2 + A_M2) +
  //   mu_I*(I_FG + I_MG  + I_FSW + I_MC + I_F2 + I_M2 +
  //           T1_FG + T1_MG + T1_FSW +T1_MC + T1_F2 + T1_M2 +
  //           T2_FG + T2_MG + T2_FSW +T2_MC + T2_F2 + T2_M2) +
  //   mu_T*(T1S_FG + T1S_MG + T1S_FSW + T1S_MC + T1S_F2 + T1S_M2 +
  //           T2S_FG + T2S_MG + T2S_FSW + T2S_MC + T2S_F2 + T2S_M2) +
  //   mu_A*(A_FG + A_MG + A_FSW + A_MC + A_F2 + A_M2);
  //
  // dx[49] = mu_I*(I_FG + I_MG  + I_FSW + I_MC + I_F2 + I_M2 +
  //                   T1_FG + T1_MG + T1_FSW +T1_MC + T1_F2 + T1_M2 +
  //                   T2_FG + T2_MG + T2_FSW +T2_MC + T2_F2 + T2_M2) +
  //   mu_T*(T1S_FG + T1S_MG + T1S_FSW + T1S_MC + T1S_F2 + T1S_M2 +
  //           T2S_FG + T2S_MG + T2S_FSW + T2S_MC + T2S_F2 + T2S_M2) +
  //   mu_A*(A_FG + A_MG + A_FSW + A_MC + A_F2 + A_M2);

  /* useful output (for model fitting) */
  double Prev_F = ((I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
             I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
             I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2)/N_F)*100;

  double Prev_M = ((I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
             I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
             I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2)/N_M)*100;

  double Prev = ((I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
             I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
             I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2) +
            (I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
               I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
               I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2))*100/(N_F + N_M);

  out[1] = N_F;
  out[2] = N_M;
  out[3] = Prev_F;
  out[4] = Prev_M;
  out[5] = Prev;

  out.names() = outnames;
  return out;
};


/* ################################################################################
#   model for testing
################################################################################ */

/* model */
// [[Rcpp::export]]
Rcpp::List hiv_fsw(const double t, const Rcpp::NumericVector& state, const Rcpp::List& pars){

  /* out */
  Rcpp::List out(6);

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

    sigma = RLIST("sigma.1");

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

    sigma = RLIST("sigma.2");

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
  } else {

    sigma = NAN;

    con_FG = NAN;
    con_MG = NAN;
    con_MC = NAN;
    con_FSW = NAN;
    con_F2 = NAN;
    con_M2 = NAN;
    circum = NAN;
    mu_A = NAN;
    mu_T = NAN;
    mu_I = NAN;
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
  } else {
    sup1rate = NAN; /* to shut up the compiler... */
    falloff1 = NAN;
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

  double con_eff = RLIST("con_eff");
  double circum_eff = RLIST("circum_eff");

  double birth = RLIST("birth");
  double rho = RLIST("rho");

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

  /* these are death rates (all cause and HIV specific mortality) */
  // double D = state[48];
  // double D_HIV = state[49];

  /* aggregated pop sizes */
  double N_FG = S_FG+IV_FG+I_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG;
  double N_MG = S_MG+IV_MG+I_MG+T1_MG+T2_FG+T1S_MG+T2S_MG+A_MG;

  double N_FSW = S_FSW+IV_FSW+I_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW;
  double N_MC = S_MC+IV_MC+I_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC;

  double N_F2 = S_F2+IV_F2+I_F2+T1_F2+T2_F2+T1S_F2+T2S_F2+A_F2;
  double N_M2 = S_M2+IV_M2+I_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2;

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

  /* general population */
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

  /* forces of infection */

  /* general population */
  double con_circ_G = ((1.0-con_FG)+con_FG*(1.0-con_eff))*((1.0-circum) + circum*(1.0-circum_eff));
  double lambda_MG = (beta_MGI * (I_FG/N_FG) * con_circ_G) + /* FOI from I females  */
    (beta_MGV * (IV_FG/N_FG) * con_circ_G) + /* FOI from IV females */
    (beta_MGT * (T1_FG/N_FG) * con_circ_G) + /* FOI from T1 females */
    (beta_MGTS * (T1S_FG/N_FG) * con_circ_G) + /* FOI from T1 virally suppressed females */
    (beta_MGT * (T2_FG/N_FG) * con_circ_G) + /* FOI from T2 females */
    (beta_MGTS * (T2S_FG/N_FG) * con_circ_G) + /* FOI from T2 virally suppressed females */
    (beta_MGA * (A_FG/N_FG) * con_circ_G); /* FOI from AIDS females */

  double lambda_FG = (beta_FGI * (I_MG/N_MG) * con_circ_G) +
    (beta_FGV * (IV_MG/N_MG) * con_circ_G) +
    (beta_FGT * (T1_MG/N_MG) * con_circ_G) +
    (beta_FGTS * (T1S_MG/N_MG) * con_circ_G) +
    (beta_FGT * (T2_MG/N_MG) * con_circ_G) +
    (beta_FGTS * (T2S_MG/N_MG) * con_circ_G) +
    (beta_FGA * (A_MG/N_MG) * con_circ_G);

  /* FSW and C */
  double con_circ_FSW = ((1.0-con_FSW)+con_FSW*(1.0-con_eff))*((1.0-circum) + circum*(1.0-circum_eff));
  double lambda_MC = (beta_MCI*(I_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCV*(IV_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCT*(T1_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCTS*(T1S_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCT*(T2_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCTS*(T2S_FSW/N_FSW) * con_circ_FSW) +
    (beta_MCA*(A_FSW/N_FSW) * con_circ_FSW);

  double lambda_FSW = (beta_FSWI * (I_MC/N_MC) * con_circ_FSW) +
    (beta_FSWV * (IV_MC/N_MC) * con_circ_FSW) +
    (beta_FSWT * (T1_MC/N_MC) * con_circ_FSW) +
    (beta_FSWTS * (T1S_MC/N_MC) * con_circ_FSW) +
    (beta_FSWT * (T2_MC/N_MC) * con_circ_FSW) +
    (beta_FSWTS * (T2S_MC/N_MC) * con_circ_FSW) +
    (beta_FSWA * (A_MC/N_MC) * con_circ_FSW);

  /* F2 and M2 */
  double con_circ_2 = ((1.0-con_F2)+con_F2*(1.0-con_eff))*((1.0-circum) + circum*(1.0-circum_eff));
  double lambda_M2 = (beta_M2I * (I_F2/N_F2) * con_circ_2) +
    (beta_M2V * (IV_F2/N_F2) * con_circ_2) +
    (beta_M2T * (T1_F2/N_F2) * con_circ_2) +
    (beta_M2TS * (T1S_F2/N_F2) * con_circ_2) +
    (beta_M2T * (T2_F2/N_F2) * con_circ_2) +
    (beta_M2TS * (T2S_F2/N_F2) * con_circ_2) +
    (beta_M2A * (A_F2/N_F2) * con_circ_2);

  double lambda_F2 = (beta_F2I * (I_M2/N_M2) * con_circ_2) +
    (beta_F2V * (IV_M2/N_M2) * con_circ_2) +
    (beta_F2T * (T1_M2/N_M2) * con_circ_2) +
    (beta_F2TS * (T1S_M2/N_M2) * con_circ_2) +
    (beta_F2T * (T2_M2/N_M2) * con_circ_2) +
    (beta_F2TS * (T2S_M2/N_M2) * con_circ_2) +
    (beta_F2A * (A_M2/N_M2) * con_circ_2);

  /* general females */
  dx[0] = (birth * N_FG) - ((lambda_FG + mu) * S_FG); /* S_FG */
  dx[1] = (lambda_FG * S_FG) - ((mu + Tau) * IV_FG); /* IV_FG */
  dx[2] = (Tau * IV_FG) + (rho * A_FG) + (falloff1 * T1_FG) + (falloff2 * T2_FG) - ((mu_I + mu + sigma + omega) * I_FG); /* I_FG */
  dx[3] = (sigma * I_FG) + (unsup1rate * T1S_FG) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_FG); /* T1_FG */
  dx[4] = (sup1rate * T1_FG) - ((mu + mu_T + unsup1rate) * T1S_FG); /* T1S_FG */
  dx[5] = (perc2ndline * T1_FG) + (unsup2rate * T2S_FG) - ((mu + mu_T + sup2rate + falloff2) * T2_FG); /* T2_FG */
  dx[6] = (sup2rate * T2_FG) - ((mu + mu_T + unsup2rate) * T2S_FG); /* T2S_FG */
  dx[7] = (omega * I_FG) - ((mu + mu_A + rho) * A_FG); /* A_FG */

  /* general males */
  dx[8] = (birth * N_MG) - ((lambda_MG + mu) * S_MG); /* S_MG */
  dx[9] = (lambda_MG * S_MG) - ((mu + Tau) * IV_MG); /* IV_FG */
  dx[10] = (Tau * IV_MG) + (rho * A_MG) + (falloff1 * T1_MG) + (falloff2 * T2_MG) - ((mu_I + mu + sigma + omega) * I_MG); /* I_MG */
  dx[11] = (sigma * I_MG) + (unsup1rate * T1S_MG) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_MG); /* T1_MG */
  dx[12] = (sup1rate * T1_MG) - ((mu + mu_T + unsup1rate) * T1S_MG); /* T1S_MG */
  dx[13] = (perc2ndline * T1_MG) + (unsup2rate * T2S_MG) - ((mu + mu_T + sup2rate + falloff2) * T2_MG); /* T2_MG */
  dx[14] = (sup2rate * T2_MG) - ((mu + mu_T + unsup2rate) * T2S_MG); /* T2S_MG */
  dx[15] = (omega * I_MG) - ((mu + mu_A + rho) * A_MG); /* A_MG */

  /* female sex workers */
  dx[16] = (birth * N_FSW) - ((lambda_FSW + mu) * S_FSW); /* S_FSW */
  dx[17] = (lambda_FSW * S_FSW) - ((mu + Tau) * IV_FSW); /* IV_FSW */
  dx[18] = (Tau * IV_FSW) + (rho * A_FSW) + (falloff1 * T1_FSW) + (falloff2 * T2_FSW) - ((mu_I + mu + sigma + omega) * I_FSW); /* I_FSW */
  dx[19] = (sigma * I_FSW) + (unsup1rate * T1S_FSW) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_FSW); /* T1_FSW */
  dx[20] = (sup1rate * T1_FSW) - ((mu + mu_T + unsup1rate) * T1_FSW); /* T1S_FSW */
  dx[21] = (perc2ndline * T1_FSW) + (unsup2rate * T2S_FSW) - ((mu + mu_T + sup2rate + falloff2) * T2_FSW); /* T2_FSW */
  dx[22] = (sup2rate * T2_FSW) - ((mu + mu_T + unsup2rate) * T2S_FSW); /* T2S_FSW */
  dx[23] = (omega * I_FSW) - ((mu + mu_A + rho) * A_FSW); /* A_FSW */

  /* male clients */
  dx[24] = (birth * N_MC) - ((lambda_MC + mu) * S_MC); /* S_MC */
  dx[25] = (lambda_MC * S_MC) - ((mu + Tau) * IV_MC); /* IV_MC */
  dx[26] = (Tau * IV_MC) + (rho * A_MC) + (falloff1 * T1_MC) + (falloff2 * T2_MC) - ((mu_I + mu + sigma + omega) * I_MC); /* I_MC */
  dx[27] = (sigma * I_MC) + (unsup1rate * T1S_MC) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_MC); /* T1_MC */
  dx[28] = (sup1rate * T1_MC) - ((mu + mu_T + unsup1rate) * T1_MC); /* T1S_MC */
  dx[29] = (perc2ndline * T1_MC) + (unsup2rate * T2S_MC) - ((mu + mu_T + sup2rate + falloff2) * T2_MC); /* T2_MC */
  dx[30] = (sup2rate * T2_MC) - ((mu + mu_T + unsup2rate) * T2S_MC); /* T2S_MC */
  dx[31] = (omega * I_MC) - ((mu + mu_A + rho) * A_MC); /* A_MC */

  /* female 2+ */
  dx[32] = (birth * N_F2) - ((lambda_F2 + mu) * S_F2); /* S_F2 */
  dx[33] = (lambda_F2 * S_F2) - ((mu + Tau) * IV_F2); /* IV_F2 */
  dx[34] = (Tau * IV_F2) + (rho * A_F2) + (falloff1 * T1_F2) + (falloff2 * T2_F2) - ((mu_I + mu + sigma + omega) * I_F2); /* I_F2 */
  dx[35] = (sigma * I_F2) + (unsup1rate * T1S_F2) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_F2); /* T1_F2 */
  dx[36] = (sup1rate * T1_F2) - ((mu + mu_T + unsup1rate) * T1_F2); /* T1S_F2 */
  dx[37] = (perc2ndline * T1_F2) + (unsup2rate * T2S_F2) - ((mu + mu_T + sup2rate + falloff2) * T2_F2); /* T2_F2 */
  dx[38] = (sup2rate * T2_F2) - ((mu + mu_T + unsup2rate) * T2S_F2); /* T2S_F2 */
  dx[39] = (omega * I_F2) - ((mu + mu_A + rho) * A_F2); /* A_F2 */

  /* male 2+ */
  dx[40] = (birth * N_M2) - ((lambda_M2 + mu) * S_M2); /* S_M2 */
  dx[41] = (lambda_M2 * S_M2) - ((mu + Tau) * IV_M2); /* IV_M2 */
  dx[42] = (Tau * IV_M2) + (rho * A_M2) + (falloff1 * T1_M2) + (falloff2 * T2_M2) - ((mu_I + mu + sigma + omega) * I_M2); /* I_M2 */
  dx[43] = (sigma * I_M2) + (unsup1rate * T1S_M2) - ((mu + mu_I + perc2ndline + sup1rate + falloff1) * T1_M2); /* T1_M2 */
  dx[44] = (sup1rate * T1_M2) - ((mu + mu_T + unsup1rate) * T1_M2); /* T1S_M2 */
  dx[45] = (perc2ndline * T1_M2) + (unsup2rate * T2S_M2) - ((mu + mu_T + sup2rate + falloff2) * T2_M2); /* T2_M2 */
  dx[46] = (sup2rate * T2_M2) - ((mu + mu_T + unsup2rate) * T2S_M2); /* T2S_M2 */
  dx[47] = (omega * I_M2) - ((mu + mu_A + rho) * A_M2); /* A_M2 */

  // /* all-cause and HIV-specific mortality */
  // dx[48] = mu*(S_FG + IV_FG + I_FG + T1_FG + T1S_FG+ T2_FG + T2S_FG+ A_FG +
  //            S_MG + IV_MG + I_MG + T1_MG + T1S_MG + T2_MG + T2S_MG + A_MG +
  //            S_FSW + IV_FSW + I_FSW + A_FSW + T1_FSW + T1S_FSW + T2_FSW + T2S_FSW +
  //            S_MC + IV_MC + I_MC + T1_MC + T1S_MC + T2_MC + T2S_MC + A_MC +
  //            S_F2 + IV_F2 + I_F2 + T1_F2 + T1S_F2 + T2_F2 + T2S_F2 + A_F2 +
  //            S_M2 + IV_M2 + I_M2 + T1_M2 + T1S_M2 + T2_M2 + T2S_M2 + A_M2) +
  //   mu_I*(I_FG + I_MG  + I_FSW + I_MC + I_F2 + I_M2 +
  //           T1_FG + T1_MG + T1_FSW +T1_MC + T1_F2 + T1_M2 +
  //           T2_FG + T2_MG + T2_FSW +T2_MC + T2_F2 + T2_M2) +
  //   mu_T*(T1S_FG + T1S_MG + T1S_FSW + T1S_MC + T1S_F2 + T1S_M2 +
  //           T2S_FG + T2S_MG + T2S_FSW + T2S_MC + T2S_F2 + T2S_M2) +
  //   mu_A*(A_FG + A_MG + A_FSW + A_MC + A_F2 + A_M2);
  //
  // dx[49] = mu_I*(I_FG + I_MG  + I_FSW + I_MC + I_F2 + I_M2 +
  //                   T1_FG + T1_MG + T1_FSW +T1_MC + T1_F2 + T1_M2 +
  //                   T2_FG + T2_MG + T2_FSW +T2_MC + T2_F2 + T2_M2) +
  //   mu_T*(T1S_FG + T1S_MG + T1S_FSW + T1S_MC + T1S_F2 + T1S_M2 +
  //           T2S_FG + T2S_MG + T2S_FSW + T2S_MC + T2S_F2 + T2S_M2) +
  //   mu_A*(A_FG + A_MG + A_FSW + A_MC + A_F2 + A_M2);

  /* useful output (for model fitting) */
  double Prev_F = ((I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
             I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
             I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2)/N_F)*100;

  double Prev_M = ((I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
             I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
             I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2)/N_M)*100;

  double Prev = ((I_MG+IV_MG+T1_MG+T2_MG+T1S_MG+T2S_MG+A_MG+
             I_MC+IV_MC+T1_MC+T2_MC+T1S_MC+T2S_MC+A_MC+
             I_M2+IV_M2+T1_M2+T2_M2+T1S_M2+T2S_M2+A_M2) +
            (I_FG+IV_FG+T1_FG+T2_FG+T1S_FG+T2S_FG+A_FG+
               I_FSW+IV_FSW+T1_FSW+T2_FSW+T1S_FSW+T2S_FSW+A_FSW+
               I_F2+IV_F2+T1_F2+T1S_F2+T2_F2+T2S_F2+A_F2))*100/(N_F + N_M);

  out[1] = N_F;
  out[2] = N_M;
  out[3] = Prev_F;
  out[4] = Prev_M;
  out[5] = Prev;

  out.names() = outnames;
  return out;
};
