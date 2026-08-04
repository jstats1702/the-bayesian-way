#include <RcppArmadillo.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
#include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace R;
using namespace Rcpp;
using namespace std;

uword get_k (const int& i, const int& ii, const double& I)
{
        return(0.5*(I*(I-1.0) - (I-i)*(I-i-1.0)) + ii - i - 1.0);
} 

double expit (const double& x) 
{
        return(1.0/(1.0 + exp(-x)));
}

void tunning (double& del2, int& n_tun, const double& mix_rate, const int& b)
{ 
        // tunning paramter calibration
        double mr = 0.35, eps = 0.05, diff = mix_rate - mr;
        if ( b%n_tun == 0 ) {
                if ( abs(diff) > eps ) {
                        double tmp = del2, cont = 1.0;
                        do {
                                tmp = del2 + (0.1/cont) * diff;
                                cont++;
                        } while ( (tmp <= 0.0) && (cont <= 1000000.0) );
                        if ( tmp > 0.0 ) del2 = tmp;
                        n_tun = 100;
                } else {
                        n_tun += 100;
                }
        }
}

// [[Rcpp::export]]
double loglik (const double& I, const double& zeta, const mat& U, const vec& Y)
{
        double out = 0.0;
        for (uword i = 0; i < I-1; i++) {
                for (uword ii = i+1; ii < I; ii++) {
                        out += R::dbinom(Y[get_k(i, ii, I)], 1, expit(zeta - norm(U.row(i) - U.row(ii))), 1);
                }
        }
        return(out);
}

double lfcd_U (const rowvec& x, const uword& i, const double& I, const double& K, const double& sigsq, const double& zeta, const mat& U, const vec& Y)
{
        double out = -pow(norm(x), 2.0)/(2.0*sigsq);
        if (i < I-1) for (uword ii = i+1; ii < I; ii++) out += R::dbinom(Y[get_k(i, ii, I)], 1, expit(zeta - norm(x - U.row(ii))), 1);
        if (i > 0  ) for (uword ii = 0;   ii < i; ii++) out += R::dbinom(Y[get_k(ii, i, I)], 1, expit(zeta - norm(x - U.row(ii))), 1);
        return(out);
}

// [[Rcpp::export]]
List sample_U (const double& b, int n_tun_U, double del2_U, int n_U, const int& n_burn, const double& I, const double& K, const double& sigsq, const double& zeta, mat U, const vec& Y)
{
        // Metropolis step
        rowvec u_c(K), u_p(K);
        for (uword i = 0; i < I; i++) {
                u_c = U.row(i);
                u_p = u_c + sqrt(del2_U)*randn<rowvec>(K);
                if (R::runif(0, 1) < exp(lfcd_U(u_p, i, I, K, sigsq, zeta, U, Y) - lfcd_U(u_c, i, I, K, sigsq, zeta, U, Y))) {
                        U.row(i) = u_p;
                        n_U++;
                }
        }
        if (b < n_burn ) {
                double mix_rate = n_U/(b*I);
                tunning(del2_U, n_tun_U, mix_rate, b);
        }
        return List::create(Named("U")       = U,
                            Named("del2_U")  = del2_U,
                            Named("n_U")     = n_U,
                            Named("n_tun_U") = n_tun_U);
}

double lfcd_zeta (const double& x, const double& I, const double& omesq, const mat& U, const vec& Y)
{
        double out = -pow(x, 2.0)/(2.0*omesq);
        for (uword i = 0; i < I-1; i++) {
                for (uword ii = i+1; ii < I; ii++) {
                        out += R::dbinom(Y[get_k(i, ii, I)], 1, expit(x - norm(U.row(i) - U.row(ii))), 1);
                }
        }
        return(out);
}

// [[Rcpp::export]]
List sample_zeta (const double& b, int n_tun_zeta, double del2_zeta, int n_zeta, const int& n_burn, const double& I, const double& omesq, double zeta, const mat& U, const vec& Y)
{
        // Metropolis step
        double zeta_p = R::rnorm(zeta, sqrt(del2_zeta));
        if (R::runif(0, 1) < exp(lfcd_zeta(zeta_p, I, omesq, U, Y) - lfcd_zeta(zeta, I, omesq, U, Y))) {
                zeta = zeta_p; 
                n_zeta++;
        }
        if (b < n_burn) {
                double mix_rate = n_zeta/b;
                tunning(del2_zeta, n_tun_zeta, mix_rate, b);
        }
        return List::create(Named("zeta")       = zeta,
                            Named("del2_zeta")  = del2_zeta,
                            Named("n_zeta")     = n_zeta,
                            Named("n_tun_zeta") = n_tun_zeta);
}

// [[Rcpp::export]]
double sample_sigsq (const double& I, const double &K, const double& a_sig, const double& b_sig, const mat& U)
{
        return(1.0/R::rgamma(a_sig + 0.5*K*I, 1.0/(b_sig + 0.5*accu(pow(U, 2)))));
}

// [[Rcpp::export]]
double sample_omesq (const double& a_ome, const double& b_ome, const double& zeta)
{
        return(1.0/R::rgamma(a_ome + 0.5, 1.0/(b_ome + 0.5*pow(zeta, 2))));
}

// [[Rcpp::export]]
vec sample_Y (const double& I, const double& zeta, const mat& U, const vec& na_indices, vec Yna)
{
        // sample NA values in Y
        uword k;
        for (uword i = 0; i < I-1; i++) {
                for (uword ii = i+1; ii < I; ii++) {
                        k = get_k(i, ii, I);
                        if (na_indices[k] == true) Yna[k] = R::rbinom(1, expit(zeta - norm(U.row(i) - U.row(ii))));
                }
        }
        return(Yna);
}

// [[Rcpp::export]]
rowvec interaction_probs0 (const double& I, const double& K, const double& B, const vec& zeta_chain, const mat& U_chain) 
{
        double zeta;
        rowvec u_i(K), u_ii(K), out(0.5*I*(I-1.0), fill::zeros);
        for (uword b = 0; b < B; b++) {
                zeta = zeta_chain[b];
                for (uword i = 0; i < I-1; i++) {
                        for (uword ii = i+1; ii < I; ii++) {
                                for (uword k = 0; k < K; k++) {
                                        u_i [k] = U_chain.at(b, k*I + i );
                                        u_ii[k] = U_chain.at(b, k*I + ii);
                                }
                                out[get_k(i, ii, I)] += expit(zeta - norm(u_i - u_ii))/B;        
                        }
                }
        }
        return(out);
}

// [[Rcpp::export]]
mat simulate_data (const double& I, const double& zeta, const mat& U)
{
        mat Y(I, I, fill::zeros);
        for (uword i = 0; i < I-1; i++) {
                for (uword ii = i+1; ii < I; ii++) {
                         Y.at(i, ii) = R::rbinom(1, expit(zeta - norm(U.row(i) - U.row(ii))));
                }
        }
        Y = Y + Y.t();
        return(Y);
}