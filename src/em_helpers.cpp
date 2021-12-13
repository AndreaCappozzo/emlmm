#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


const double log2pi = std::log(2.0 * M_PI);

double dmvnrm_arma(arma::mat x,
                   arma::rowvec mean,
                   arma::mat sigma,
                   bool logd = false) {
  //int n = x.n_rows;
  int xdim = x.n_cols;
  double out;
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  arma::vec z = rooti * arma::trans( x - mean) ;
  out      = constants - 0.5 * arma::sum(z%z) + rootisum;

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
Rcpp::List estep_lmm_cpp(arma::vec res_fixed, arma::mat Z,
                      arma::vec group_indicator, arma::mat inv_Omega,
                      double sigma2, int J)
{
  int N = Z.n_rows;
  int q = Z.n_cols;

  // Containers
  arma::mat est_second_moment(q,q);
  double est_second_moment_error = 0.0;
  arma::vec raneff_i(N);
  arma::mat mu_raneff(J,q);

  //  Fill containers
  raneff_i.zeros();
  est_second_moment.zeros();

  arma::mat Omega=inv_Omega.i();
  for ( int j = 1; j < (J+1); j++ ) {
    arma::uvec rows_j = find(group_indicator==j);
    arma::mat Z_j=Z.rows(rows_j);
    arma::vec res_fixed_j=res_fixed(rows_j);
    arma::mat first_piece = (Z_j.t()*Z_j)/ sigma2 + inv_Omega;
    arma::mat Gamma_j = first_piece.i();
    arma::vec mu_j = (Gamma_j*Z_j.t()*res_fixed_j)/sigma2;
    mu_raneff.row((j-1))=mu_j.t();
    raneff_i(rows_j)=Z_j * mu_j;
    est_second_moment += Gamma_j + mu_j * mu_j.t();
    // Rcout << est_second_moment;
    est_second_moment_error += as_scalar(arma::trace(Z_j*Gamma_j*Z_j.t()));
  }

  return Rcpp::List::create( Named("est_second_moment") = est_second_moment,
                             Named("est_second_moment_error") = est_second_moment_error,
                             Named("mu_raneff") = mu_raneff,
                             Named("raneff_i") = raneff_i);
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
double log_lik_lmm_cpp(arma::vec y, arma::mat Z, arma::mat X,
                      arma::vec group_indicator, arma::vec beta, arma::mat Omega,
                      double sigma2, int J)
{

  // Containers
  double log_lik_lmm = 0.0;

  for ( int j = 1; j < (J+1); j++ ) {
    arma::uvec rows_j = find(group_indicator==j);
    arma::mat Z_j=Z.rows(rows_j);
    arma::mat X_j=X.rows(rows_j);
    arma::vec y_j=y(rows_j);
    int n_j = X_j.n_rows;
    arma::mat diag_sigma(n_j,n_j);
    diag_sigma.eye();
    arma::mat G_j=Z_j * Omega * Z_j.t()+diag_sigma/pow(sigma2,-1);
    arma::vec mu_j=X_j*beta;
    log_lik_lmm+=dmvnrm_arma(y_j.t(),mu_j.t(),G_j,true);
  }

  return log_lik_lmm;
}


// MULTIVARIATE RESPONSE

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
Rcpp::List estep_mlmm_cpp(arma::vec vec_res_fixed, arma::mat Z,
                         arma::vec group_indicator,
                         arma::vec vec_group_indicator,
                         arma::mat inv_Psi,
                         arma::mat I_r,
                         arma::mat inv_Sigma, int r, int J)
{
  int N = Z.n_rows;
  int q = Z.n_cols;

  // Containers
  arma::mat est_second_moment(q*r,q*r);
  arma::mat est_second_moment_error(r,r);
  arma::vec vec_raneff_i(N*r);
  arma::mat mat_mu_raneff(q*r,J);

  //  Fill containers
  vec_raneff_i.zeros();
  est_second_moment.zeros();

  for ( int j = 1; j < (J+1); j++ ) {
    arma::uvec rows_j = find(group_indicator==j);
    arma::uvec vec_rows_j = find(vec_group_indicator==j);
    arma::mat Z_j=Z.rows(rows_j);
    int n_j = Z_j.n_rows;
    arma::mat I_nj(n_j,n_j);
    I_nj.eye();

    arma::vec vec_res_fixed_j=vec_res_fixed(vec_rows_j);
    arma::mat common_component_j = kron(I_r,Z_j.t())*(kron(inv_Sigma,I_nj));
    arma::mat first_piece= common_component_j * kron(I_r,Z_j) + inv_Psi;
    arma::mat Gamma_j = first_piece.i();
    arma::vec vec_Delta_j= Gamma_j * common_component_j * vec_res_fixed_j;
    mat_mu_raneff.col((j-1))=vec_Delta_j;
    vec_raneff_i(vec_rows_j)=kron(I_r,Z_j) * vec_Delta_j;
    est_second_moment += Gamma_j + vec_Delta_j * vec_Delta_j.t();
    arma::mat var_ei = kron(I_r,Z_j) * Gamma_j * kron(I_r,Z_j.t());
    // arma::uvec slice_rows=arma::regspace(0, 1, n_j);
    arma::uvec slice_rows = regspace<arma::uvec>(0, 1, (n_j-1));

    for(int row=0; row<r; row++ ){
    arma::uvec slice_cols = regspace<arma::uvec>(0, 1, (n_j-1));
      for(int col=0; col<r; col++ ){
        est_second_moment_error(row,col)+=trace(var_ei.submat(slice_rows,slice_cols));
        slice_cols+=(n_j);
      }
      slice_rows+= (n_j);
    }
  }

  return Rcpp::List::create( Named("est_second_moment") = est_second_moment,
                             Named("est_second_moment_error") = est_second_moment_error,
                             Named("mat_mu_raneff") = mat_mu_raneff,
                             Named("vec_raneff_i") = vec_raneff_i);
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
double log_lik_mlmm_cpp(arma::vec vec_Y, arma::mat Z, arma::vec vec_XB,
                       arma::vec group_indicator, arma::vec vec_group_indicator, arma::mat PSI,
                       arma::mat SIGMA, arma::mat I_r, int J)
{

  // Containers
  double log_lik_lmm = 0.0;

  for ( int j = 1; j < (J+1); j++ ) {
    arma::uvec rows_j = find(group_indicator==j);
    arma::uvec vec_rows_j = find(vec_group_indicator==j);
    arma::mat Z_j=Z.rows(rows_j);
    arma::vec vec_XB_j=vec_XB(vec_rows_j);
    arma::vec vec_y_j=vec_Y(vec_rows_j);
    int n_j = Z_j.n_rows;
    arma::mat I_nj(n_j,n_j);
    I_nj.eye();
    arma::mat G_j=kron(I_r,Z_j)*PSI*kron(I_r, Z_j.t()) + kron(SIGMA,I_nj);
    log_lik_lmm+=dmvnrm_arma(vec_y_j.t(),vec_XB_j.t(),G_j,true);
  }

  return log_lik_lmm;
}
