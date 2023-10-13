// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>

using namespace arma;
using namespace Rcpp;

Rcpp::List SER(arma::vec y,arma::cube X, arma::cube XX,double sigma2,double sigma2_0,arma::vec pi_0) {
  double n = y.n_elem;
  double p = X.n_slices;
  
  arma::cube S_tilde = zeros(2,2,p);
  arma::mat s2_tilde = zeros(2,p);
  arma::mat mu_tilde = zeros(2,p);
  
  arma::vec lBF = zeros(p);
  for (int j = 0; j < p; j++) {
    arma::mat X_j = X.slice(j);
    arma::mat XX_j = XX.slice(j);
    
    arma::mat Omega_tilde = 1/sigma2*XX_j + 1/sigma2_0*eye(2,2);
    S_tilde.slice(j) = inv_sympd(Omega_tilde);
    s2_tilde.col(j) = S_tilde.slice(j).diag();
    mu_tilde.col(j) = 1/sigma2*S_tilde.slice(j)*X_j.t()*y;
    
    lBF(j) = -log(sigma2_0)-0.5*log(det(Omega_tilde))+0.5*as_scalar(mu_tilde.col(j).t()*Omega_tilde*mu_tilde.col(j));
  }
  
  double maxlbf = max(lBF);
  arma::vec w = exp(lBF - maxlbf);
  arma::vec w_weighted = w%pi_0;
  double weighted_sum_w = sum(w_weighted);
  arma::vec alpha = w_weighted/weighted_sum_w;
  
  double lbf_model = maxlbf + log(weighted_sum_w);
  double loglik = lbf_model + sum(-0.5*log(2*M_PI)-0.5*log(sigma2)-0.5*y%y/sigma2);
  
  return Rcpp::List::create(
    Rcpp::Named("mu_tilde") = mu_tilde.t(),
    Rcpp::Named("s2_tilde") = s2_tilde.t(),
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("loglik") = loglik
  );
}

// [[Rcpp::export]]
Rcpp::List SuSiE_group(arma::vec y,arma::mat X1,arma::mat X2,double L,double f=1,double delta=0.001,int maxIt=100) {
  
  if (L > 0) {
    double n = X1.n_rows;
    double p = X1.n_cols;
    
    double A_sigma = 0.001;
    double B_sigma = 0.001;
    double sigma2_0 = 1;
    arma::vec pi_0 = 1/p*ones(p);
    
    arma::mat Mu_q_beta1 = zeros(p,L);
    arma::mat Mu_q_beta2 = zeros(p,L);
    arma::mat Sigma_q_beta1 = zeros(p,L);
    arma::mat Sigma_q_beta2 = zeros(p,L);
    arma::mat Alpha = zeros(p,L);
    double mu_q_sigma2 = 1;
    double mu_q_sigma2inv = 1;
    
    arma::vec elbo = zeros(maxIt);
    arma::vec KL = zeros(L);
    
    arma::vec X1_2 = sum(X1%X1,0).t();
    arma::vec X2_2 = sum(X2%X2,0).t();
    arma::cube X = zeros(n,2,p);
    arma::cube XX = zeros(2,2,p);
    for (int j = 0; j < p; j++) {
      X.slice(j) = join_rows(X1.col(j),X2.col(j));
      XX.slice(j) = X.slice(j).t()*X.slice(j);
    }
    
    arma::vec r = y - X1*sum(Mu_q_beta1,1) - X2*sum(Mu_q_beta2,1);
    
    int it = 0;
    int conv = 0;
    while (conv == 0) {
      it = it+1;
      
      for (int l = 0; l < L; l++) {
        
        arma::vec r_l = r + X1*Mu_q_beta1.col(l) + X2*Mu_q_beta2.col(l);
        
        Rcpp::List ser = SER(r_l,X,XX,1/(f*mu_q_sigma2inv),sigma2_0,pi_0);
        arma::vec A_ser = ser["alpha"];
        arma::mat M_ser = ser["mu_tilde"];
        arma::mat S_ser = ser["s2_tilde"];
        double loglik_ser = ser["loglik"];
        
        Alpha.col(l) = A_ser;
        
        Mu_q_beta1.col(l) = M_ser.col(0)%A_ser;
        Mu_q_beta2.col(l) = M_ser.col(1)%A_ser;
        
        Sigma_q_beta1.col(l) = A_ser%S_ser.col(0) + A_ser%(1-A_ser)%M_ser.col(0)%M_ser.col(0);
        Sigma_q_beta2.col(l) = A_ser%S_ser.col(1) + A_ser%(1-A_ser)%M_ser.col(1)%M_ser.col(1);
        
        r = r_l - X1*Mu_q_beta1.col(l) - X2*Mu_q_beta2.col(l);
        
        KL(l) = - loglik_ser +
          - 0.5*n*log(2*M_PI*mu_q_sigma2) +
          - 0.5/mu_q_sigma2*(sum(r%r) + sum(X1_2%sum(Sigma_q_beta1,1)) + sum(X2_2%sum(Sigma_q_beta2,1)));
          
      }
      
      double ERSS = sum(r%r) + sum(X1_2%sum(Sigma_q_beta1,1)) + sum(X2_2%sum(Sigma_q_beta2,1));
      double B_q_sigma2 = B_sigma + 0.5*f*ERSS;
      double A_q_sigma2 = 0.5*n*f + A_sigma;
      mu_q_sigma2 = B_q_sigma2/(A_q_sigma2-1);
      mu_q_sigma2inv = A_q_sigma2/B_q_sigma2;
      
      elbo(it-1) = - 0.5*n*log(2*M_PI*mu_q_sigma2) - 0.5/mu_q_sigma2*ERSS - sum(KL);
      
      if (it > 1) {
        if (abs(elbo(it-1)-elbo(it-2)) < delta) conv = 1;
      }
      if (it == maxIt) conv = 1;
      
    }
    
    arma::vec elbo_iter = elbo(span(0,it-1));
    double elbo_final = elbo(it-1);
    
    arma::vec B1_out = sum(Mu_q_beta1,1);
    arma::vec B2_out = sum(Mu_q_beta2,1);
    arma::mat Alpha_out = Alpha;
    
    arma::vec yhat = X1*B1_out + X2*B2_out;
    
    return Rcpp::List::create(
      Rcpp::Named("B1") = B1_out,
      Rcpp::Named("B2") = B2_out,
      Rcpp::Named("alpha") = Alpha_out,
      Rcpp::Named("sigma2") = mu_q_sigma2,
      Rcpp::Named("yhat") = yhat,
      Rcpp::Named("elbo_iter") = elbo_iter,
      Rcpp::Named("elbo") = elbo_final
    );
  }
  
  if (L == 0) {
    double n = y.n_elem;
    
    double A_sigma = 0.001;
    double B_sigma = 0.001;
    
    double mu_q_sigma2 = 1;
    double mu_q_sigma2inv = 1;
    
    double B_q_sigma2 = B_sigma + 0.5*f*sum(y%y);
    double A_q_sigma2 = 0.5*n*f + A_sigma;
    mu_q_sigma2 = B_q_sigma2/(A_q_sigma2-1);
    
    double elbo = - 0.5*n*log(2*M_PI*mu_q_sigma2) - 0.5/mu_q_sigma2*sum(y%y);
    
    arma::vec yhat = zeros(n);
    
    return Rcpp::List::create(
      Rcpp::Named("yhat") = yhat,
      Rcpp::Named("sigma2") = mu_q_sigma2,
      Rcpp::Named("elbo") = elbo
    );
  }
  
}
