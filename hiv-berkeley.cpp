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
    sup1rate = RLIST(pars["sup.1.rate.pre2013"]);
    falloff1 = RLIST(pars["fall.off.1.pre2013"]);
  } else if(t >= 2013.0 && t < 2015.0){
    sup1rate = RLIST(pars["sup.1.rate.pre2015"]);
    falloff1 = RLIST(pars["fall.off.1.pre2015"]);
  } else if(t >= 2015.0){
    sup1rate = RLIST(pars["sup.1.rate.post2015"]);
    falloff1 = RLIST(pars["fall.off.1.post2015"]);
  }

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

  
  //            S_FSW=NF_0*pars$Prop_FSW, IV_FSW=0, I_FSW=NF_0*pars$Prop_FSW*pars$Prop_IFSW, T1_FSW=0,T1S_FSW=0,T2_FSW=0,T2S_FSW=0, A_FSW=0,
  //            S_MC=NM_0*pars$Prop_MC, IV_MC=0, I_MC=NM_0*pars$Prop_MC*pars$Prop_IMC, T1_MC=0,T1S_MC=0,T2_MC=0,T2S_MC=0, A_MC=0,
  //            S_F2=NF_0*pars$Prop_F2, IV_F2=0 ,I_F2=NF_0*pars$Prop_F2*pars$Prop_IF2,T1_F2=0,T1S_F2=0,T2_F2=0,T2S_F2=0, A_F2 = 0,
  //            S_M2=NM_0*pars$Prop_M2, IV_M2=0, I_M2=NM_0*pars$Prop_IM2*pars$Prop_M2, T1_M2=0,T1S_M2=0,T2_M2=0,T2S_M2=0, A_M2 = 0,D=0, D_HIV


};
