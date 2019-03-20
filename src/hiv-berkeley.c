// /* ################################################################################
// #
// #   Translation of the basic model (by Partow Imani) into C (for deSolve)
// #   Sean Wu
// #   March 2019
// #
// ################################################################################ */
//
// /* R includes */
// #include <R.h>
// #include <Rinternals.h>
// #include <Rdefines.h>
// #include <R_ext/Rdynload.h>
//
//
// /* ################################################################################
// # pull something out of a VECSEXP (list with named elements)
// # list: the VECSEXP list
// # str: name of the element to get
// ################################################################################ */
//
// SEXP getListElement(SEXP list, const char *str) {
//     SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
//     int i;
//     for (i = 0; i < length(list); i++){
//       if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0){
//           elmt = VECTOR_ELT(list, i);
//           break;
//       }
//     }
//     if (elmt == R_NilValue){
//       return NULL;
//     } else {
//       return elmt;
//     }
// }
//
//
// /* ################################################################################
// # parameters
// ################################################################################ */
//
// /* 1: t < 2015, 2: t >= 2015 */
// static double sigma_1; /* rate infecteds get 1st line treatment */
// static double sigma_2;
//
// static double con_FG_1; /* P(condom use | MG encounters with FG) */
// static double con_FG_2;
// static double con_MG_1; /* P(condom use | FG encounters with MG) */
// static double con_MG_2;
// static double con_MC_1; /* P(condom use | FSW encounters with MC) */
// static double con_MC_2;
// static double con_FSW_1; /* P(condom use | MC encounters with FSW) */
// static double con_FSW_2;
// static double con_F2_1; /* P(condom use | M2 encounters with F2) */
// static double con_F2_2;
// static double con_M2_1; /* P(condom use | F2 encounters with M2) */
// static double con_M2_2;
//
// static double circum_1; /* prop. circumcision */
// static double circum_2;
//
// static double mu_A_1; /* additional mortality due to AIDS */
// static double mu_A_2;
// static double mu_T_1; /* additional mortality when on treatment */
// static double mu_T_2;
// static double mu_I_1; /* additional mortality in chronic infection period */
// static double mu_I_2;
//
// /* 1: t < 2013, 2: 2013 <= t < 2015, 3: t >= 2015 */
// static double sup1rate_1; /* rate at which individuals become suppressed on T1 within 6 months of starting */
// static double sup1rate_2;
// static double sup1rate_3;
// static double falloff1_1; /* rate at which individuals fall off treatment on 1st line */
// static double falloff1_2;
// static double falloff1_3;
//
// /* constants */
// static double Tau;
// static double mu;
// static double omega;
// static double falloff2;
// static double sup2rate;
// static double unsup2rate;
// static double unsup1rate;
// static double perc2ndline;
//
// static double con_eff;
// static double circum_eff;
//
// static double birth;
// static double rho;
//
// /* contact constraints */
// static double C_FGMG;
// static double C_FSWC;
// static double C_F2M2;
//
// /* transmission rates */
// static double P_transmission;
// static double HighV_factor;
// static double T_factor;
//
//
// /* ################################################################################
// #   initializer
// ################################################################################ */
//
// void initmod(void (* odeparms)(int *, double *)){
//
//   DL_FUNC get_deSolve_gparms;
//   get_deSolve_gparms = R_GetCCallable("deSolve","get_deSolve_gparms");
//   SEXP gparms = get_deSolve_gparms();
//
//   /* rate infecteds get 1st line treatment */
//   SEXP sigma_1_R = getListElement(gparms, "sigma.1");
//   sigma_1 = asReal(sigma_1_R);
//
//   SEXP sigma_2_R = getListElement(gparms, "sigma.2");
//   sigma_2 = asReal(sigma_2_R);
//
//   /* P(condom use | MG encounters with FG) */
//   SEXP con_FG_1_R = getListElement(gparms, "con_FG.1");
//   con_FG_1 = asReal(con_FG_1);
//
//   SEXP con_FG_2_R = getListElement(gparms, "con_FG.2");
//   con_FG_2 = asReal(con_FG_2);
//
//
//
//
//   static double con_MG_1; /* P(condom use | FG encounters with MG) */
//   static double con_MG_2;
//   static double con_MC_1; /* P(condom use | FSW encounters with MC) */
//   static double con_MC_2;
//   static double con_FSW_1; /* P(condom use | MC encounters with FSW) */
//   static double con_FSW_2;
//   static double con_F2_1; /* P(condom use | M2 encounters with F2) */
//   static double con_F2_2;
//   static double con_M2_1; /* P(condom use | F2 encounters with M2) */
//   static double con_M2_2;
//
//   static double circum_1; /* prop. circumcision */
//   static double circum_2;
//
//   static double mu_A_1; /* additional mortality due to AIDS */
//   static double mu_A_2;
//   static double mu_T_1; /* additional mortality when on treatment */
//   static double mu_T_2;
//   static double mu_I_1; /* additional mortality in chronic infection period */
//   static double mu_I_2;
//
//   /* 1: t < 2013, 2: 2013 <= t < 2015, 3: t >= 2015 */
//   static double sup1rate_1; /* rate at which individuals become suppressed on T1 within 6 months of starting */
//   static double sup1rate_2;
//   static double sup1rate_3;
//   static double falloff1_1; /* rate at which individuals fall off treatment on 1st line */
//   static double falloff1_2;
//   static double falloff1_3;
//
//   /* constants */
//   static double Tau;
//   static double mu;
//   static double omega;
//   static double falloff2;
//   static double sup2rate;
//   static double unsup2rate;
//   static double unsup1rate;
//   static double perc2ndline;
//
//   static double con_eff;
//   static double circum_eff;
//
//   static double birth;
//   static double rho;
//
//   /* contact constraints */
//   static double C_FGMG;
//   static double C_FSWC;
//   static double C_F2M2;
//
//   /* transmission rates */
//   static double P_transmission;
//   static double HighV_factor;
//   static double T_factor;
//
//
//
// }
