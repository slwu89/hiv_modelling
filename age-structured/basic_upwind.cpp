/* ################################################################################
# 
#   Testing a basic first-order upwind scheme for McKendrick-Von Foerster equation
#   see https://en.wikipedia.org/wiki/Von_Foerster_equation for details.
#   Using implementation described in:
#   Abia, L. M., Angulo, O., & Lòpez-Marcos, J. C. (2002). Numerical Solution of Stuctured Population Models, (May 2014).
# 
################################################################################ */

#include <Rcpp.h>
#include <algorithm>
 
 // [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

/* fertility function */
double beta(const Rcpp::NumericMatrix& U, const double t);

double beta(const Rcpp::NumericMatrix& U, const double t){
  double N = std::accumulate(U.row(t).begin(),U.row(t).end(),0.0);
  return (1.0/52.0) * N;
};

/* mortality function */
double mu(const double a, const double t);

double mu(const double a, const double t){
  return (1.0/52.0);
};

/* simulation */
// [[Rcpp::export]]
void vonFoerster_upwind(Rcpp::NumericMatrix& U, const double dt, const double da, const bool prog = true){
  
  size_t t_cells_n = U.nrow();
  size_t a_cells_n = U.ncol();
  
  double v = dt/da;
  
  Progress p(t_cells_n, prog);
  
  for(size_t t = 0; t < t_cells_n-1; t++){
    
    /* boundary condition at u(0,t) */
    size_t a = 0;
    U.at(t+1,a) = beta(U,t)*dt;
    
    for(a = 1; a < a_cells_n; a++){
      
      U.at(t+1,a) = U.at(t,a) - (v*(U.at(t,a) - U.at(t,(a-1)))) - (dt*mu(a,t)*U.at(t,a));
      
    }
    
    p.increment();
  }
};
