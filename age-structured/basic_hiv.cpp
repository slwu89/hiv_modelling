/* ################################################################################
# 
#   Basic HIV/AIDS model with demography
#   Using model described in:
#   May, R. M., Anderson, R. M., & McLean, A. R. (1988). Possible demographic consequences of HIV/AIDS epidemics. I. assuming HIV infection always leads to AIDS. Mathematical Biosciences, 90(1–2), 475–505. https://doi.org/10.1016/0025-5564(88)90079-X
# 
################################################################################ */

#include <Rcpp.h>
#include <algorithm>
 
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

/* fertility function */
double beta(const Rcpp::NumericMatrix& U, const double t);

/* mortality function */
double mu(const double a, const double t);

/* force of infection */
double lambda(const Rcpp::NumericMatrix& X, const Rcpp::NumericMatrix& Y, const double a, const double t);

/* partner function: P(person of age a chooses sexual partner of age a_prime) */
double p(const double a, const double a_prime);

/* simulation */
// [[Rcpp::export]]
void basicHIV_upwind(Rcpp::NumericMatrix& X, Rcpp::NumericMatrix& Y, const double dt, const double da, const bool prog = true){
  
  size_t t_cells_n = X.nrow();
  size_t a_cells_n = X.ncol();
  
  if(t_cells_n != Y.nrow() || a_cells_n != Y.ncol()){
    Rcpp::stop("X and Y must be the same dimension");
  }
  
  double v = dt/da;
  
  Progress p(t_cells_n, prog);
  
  for(size_t t = 0; t < t_cells_n-1; t++){
    
    /* boundary condition at u(0,t) */
    size_t a = 0;
    // U.at(t+1,a) = beta(U,t)*dt;
    // 
    // for(a = 1; a < a_cells_n; a++){
    //   
    //   U.at(t+1,a) = U.at(t,a) - (v*(U.at(t,a) - U.at(t,(a-1)))) - (dt*mu(a,t)*U.at(t,a));
    //   
    // }
    // 
    p.increment();
  }
  
};
