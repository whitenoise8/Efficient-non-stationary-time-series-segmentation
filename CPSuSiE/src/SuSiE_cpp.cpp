// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>

using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List SER(arma::vec y,arma::cube X, arma::cube XX,double sigma2,double sigma2_0,arma::vec pi_0) {
  double n = y.n_elem;
  double p = X.n_slices;
  
  arma::cube S_tilde = zeros(2,2,p);
  arma::mat mu_tilde = zeros(2,p);
  
  arma::vec lBF = zeros(p);
  for (int j = 0; j < p; j++) {
    arma::mat X_j = X.slice(j);
    arma::mat XX_j = XX.slice(j);
    
    arma::mat Omega_tilde = 1/sigma2*XX_j + 1/sigma2_0*eye(2,2);
    S_tilde.slice(j) = inv_sympd(Omega_tilde);
    mu_tilde.col(j) = 1/sigma2*S_tilde.slice(j)*X_j.t()*y;
    
    lBF(j) = -log(sigma2_0)-0.5*log(det(Omega_tilde))+0.5*as_scalar(mu_tilde.col(j).t()*Omega_tilde*mu_tilde.col(j));
  }
  
  double maxlbf = max(lBF);
  arma::vec w = exp(lBF - maxlbf);
  arma::vec w_weighted = w%pi_0;
  double weighted_sum_w = sum(w_weighted);
  arma::vec alpha = w_weighted/weighted_sum_w;
  
  arma::mat S_hat = S_tilde.slice(index_max(alpha));
  arma::vec mu_hat = mu_tilde.col(index_max(alpha));
  
  double lbf_model = maxlbf + log(weighted_sum_w);
  double loglik = lbf_model + sum(-0.5*log(2*M_PI)-0.5*log(sigma2)-0.5*y%y/sigma2);
  
  return Rcpp::List::create(
    Rcpp::Named("mu_b") = mu_hat,
    Rcpp::Named("S_b") = S_hat,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("loglik") = loglik
  );
}

// [[Rcpp::export]]
Rcpp::List SuSiE_group(arma::vec y,arma::mat X1,arma::mat X2,double L,double delta=0.001,int maxIt=100) {
  
  double n = X1.n_rows;
  double p = X1.n_cols;
  
  double sigma2 = 1;
  double sigma2_0 = 1;
  arma::vec pi_0 = 1/p*ones(p);
  
  arma::mat B1 = zeros(p,L);
  arma::mat B2 = zeros(p,L);
  arma::mat B1_sq = zeros(p,L);
  arma::mat B2_sq = zeros(p,L);
  arma::mat A = zeros(p,L);
  
  arma::vec elbo = zeros(maxIt);
  arma::vec KL = zeros(L);
  
  arma::mat X1_2 = X1%X1;
  arma::mat X2_2 = X2%X2;
  arma::cube X = zeros(n,2,p);
  arma::cube XX = zeros(2,2,p);
  for (int j = 0; j < p; j++) {
    X.slice(j) = join_rows(X1.col(j),X2.col(j));
    XX.slice(j) = X.slice(j).t()*X.slice(j);
  }
  
  arma::vec r = y - X1*sum(B1,1) - X2*sum(B2,1);
  
  int it = 0;
  int conv = 0;
  while (conv == 0) {
    it = it+1;

    if (L > 0) {
      for (int l = 0; l < L; l++) {
        
        arma::vec r_l = r + X1*B1.col(l) + X2*B2.col(l);
        
        Rcpp::List ser = SER(r_l,X,XX,sigma2,sigma2_0,pi_0);
        arma::vec alpha_ser = ser["alpha"];
        arma::vec mu_b_ser = ser["mu_b"];
        arma::mat S_b_ser = ser["S_b"];
        double loglik_ser = ser["loglik"];
        
        A.col(l) = alpha_ser;
        
        B1.col(l) = mu_b_ser(0)*alpha_ser;
        B2.col(l) = mu_b_ser(1)*alpha_ser;
        
        B1_sq.col(l) = (S_b_ser(0,0)+mu_b_ser(0)*mu_b_ser(0))*alpha_ser;
        B2_sq.col(l) = (S_b_ser(1,1)+mu_b_ser(1)*mu_b_ser(1))*alpha_ser;
        
        KL[l] = - loglik_ser +
          - 0.5*n*log(2*M_PI*sigma2) +
          - 0.5/sigma2*(sum(y%y) - 2*sum(y%(X1*B1.col(l)+X2*B2.col(l))) + sum(X1_2*B1_sq.col(l)+X2_2*B2_sq.col(l)));
          
          r = r_l - X1*B1.col(l) - X2*B2.col(l);
      }
    }
    
    arma::mat XB = (X1*B1+X2*B2);
    double ERSS = sum(r%r); //sum(r%r) + sum(sum(X1_2*B1_sq+X2_2*B2_sq,0)) - sum(sum(XB%XB,0));
    sigma2 = ERSS/n;
    
    elbo(it-1) = - 0.5*n*log(2*M_PI*sigma2) - 0.5/sigma2*ERSS - sum(KL);
    
    if (it > 1) {
      if (abs(elbo(it-1)-elbo(it-2)) < delta) conv = 1;
      if (it == maxIt) conv = 1;
    }
    
  }
  
  arma::vec elbo_iter = elbo(span(0,it-1));
  double elbo_final = elbo(it-1);
  
  arma::vec B1_out = sum(B1,1);
  arma::vec B2_out = sum(B2,1);
  arma::mat alpha_out = A;
  
  arma::vec yhat = X1*B1_out + X2*B2_out;
  
  return Rcpp::List::create(
    Rcpp::Named("B1") = B1_out,
    Rcpp::Named("B2") = B2_out,
    Rcpp::Named("sigma2") = sigma2,
    Rcpp::Named("yhat") = yhat,
    Rcpp::Named("alpha") = alpha_out,
    Rcpp::Named("elbo_iter") = elbo_iter,
    Rcpp::Named("elbo") = elbo_final
  );
}
