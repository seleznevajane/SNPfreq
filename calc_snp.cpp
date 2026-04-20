#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/gamma.hpp>

// [[Rcpp::depends(BH)]]

using namespace Rcpp;

//' Calculate scaled exponential integral using asymptotic expansion for large x:
//' exp(x) * Ei(-x) = 1/x * (1 - 1/x + 2/x^2 - 6/x^3 + 24/x^4)
//' @param x In exponential growth models, (j * (j - 1) / 2) / (r * Ne0)
double scaled_expint(double x) {
  if (x > 50.0) {
    double invX = 1.0 / x;
    return invX * (1.0 - invX + (2.0 * invX * invX) - (6.0 * invX * invX * invX) +
                   (24.0 * invX * invX * invX * invX));
  } else {
    return -std::exp(x) * boost::math::expint(-x);
  }
}

//' Calculate expectations of coalescence times
 //' @param n Sample size
 //' @param population_type Model of populations size:
 //'  0 - constant - need to specify Ne0, 
 //'  1 - exponential - need to specify kappa = r * Ne0,
 //'  2 - biphasic - need to specify Ne0, NeT and T,
 //'  3 - power - need to specify kappa and beta,
 //'  4 - stepwise constant - need to specify tk and eta.
 //' @param kappa The product of current population size and the growth rate (Ne0 * r)
 //' @param Ne0 Current effective population size
 //' @param T In a biphasic model, time in history when the growth model changed
 //'  from constant to exponential
 //' @param NeT In a biphasic model, effective population size at moment T
 //' @param beta In a power model, the power in the model growth equation
 //' @param tk In a stepwise constant model, times when the population changed
 //' @param eta In a stepwise constant model, effective population sizes 
 //' at time periods tk[i-1] to tk[i]
 // [[Rcpp::export]]
 NumericVector Calculate_e_j(int n, 
                             int population_type = 0, 
                             double kappa = 1.0,
                             int Ne0 = 100,
                             int T = 10, 
                             int NeT = 10,
                             double beta = 0.5,
                             NumericVector tk = NumericVector(),
                             NumericVector eta = NumericVector()) {
   NumericVector e_j(n - 1, 0.0);
   
   if (population_type == 0) {
     for (int j = 2; j <= n; j++) {
       e_j[j - 2] = Ne0 * 2.0 / (j * (j - 1.0));
     }
   }
   
   if (population_type == 1) {
     for (int j = 2; j <= n; j++) {
       double x = (j * (j - 1.0)) / (2.0 * kappa);
       e_j[j - 2] = scaled_expint(x);
     }
   }
   
   if (population_type == 2) {
     double r = std::log(Ne0 * 1.0 / NeT) / (T * 1.0);
     for (int j = 2; j <= n; j++) {
       double x0 = (j * (j - 1.0)) / (2.0 * r * Ne0);
       double xT = x0 * std::exp(r * T);
       e_j[j - 2] = std::exp(x0) * (boost::math::expint(-xT) - boost::math::expint(-x0)) +
                    std::exp(x0 - xT) / (r * xT);
     }
   }
   
   if (population_type == 3) {
     for (int j = 2; j <= n; j++) {
       double C = (j * (j - 1.0)) / (2.0 * kappa * (1.0 - beta));
       e_j[j - 2] = std::exp(C + std::log(boost::math::tgamma(1.0 / (1.0 - beta), C)) -
                    std::log((1.0 - beta) * std::pow(C, 1.0 / (1.0 - beta))));
     }
   }
   
   if (population_type == 4) {
     int K = tk.size();
     std::vector <double> tau(K, 0.0);
     for (int k = 1; k < K - 2; k++) {
       tau[k] = tau[k-1] + 1.0 * (tk[k] - tk[k-1]) / eta[k];
     }
     for (int j = 2; j <= n; j++) {
       double C = (j * (j - 1.0)) / 2.0;
       for (int k = 1; k <= K - 1; k++) {
         e_j[j - 2] += eta[k] * std::exp(-C * tau[k-1]) * 
           (1.0 - std::exp(-C * (tk[k] - tk[k-1]) / eta[k])) / C;
       }
     }
   }
   
   return e_j;
 }
 
 double single_snp_prob(int n, int b, NumericVector e_j){
   std::vector<double> W(n + 1, 0.0);
   std::vector<double> V(n + 1, 0.0);
   
   // --- Compute W_j coefficients---
   W[2] = 6.0 / (n + 1.0);
   W[3] = 30.0 * (n - 2.0 * b) / ((n + 1.0) * (n + 2.0));
   
   for (int j = 2; j <= (n - 2); ++j) {
     W[j + 2] = - W[j] * (1.0 + j) * (3.0 + 2.0 * j) * (n - j + 0.0) /
       (j * (2.0 * j - 1) * (n + j + 1.0)) +
         W[j + 1] * (3.0 + 2.0 * j) * (n - 2.0 * b) /
           (j * (n + j + 1.0));
   }
   
   // --- Compute V_j coefficients ---
   V[2] = 6.0 * (n - 1.0) / (n + 1.0);
   
   for (int j = 4; j <= n; j += 2) {
     V[j] = (2.0 * j - 1.0) * (n - j + 1.0) * (n - j + 2.0) * V[j-2] /
       ((2.0 * j - 5.0) * (n + j - 1.0) * (n + j - 2.0));
   }
   
   // --- Final Calculation ---
   double numerator = 0.0;
   double denominator = 0.0;
   
   for (int j = 2; j <= n; ++j) {
     numerator += e_j[j - 2] * W[j];
     denominator += e_j[j - 2] * V[j];
   }
   
   return numerator / denominator;
 }
 
 
 //' Calculate SNP Probabilities using recursive coefficients
 //' @param n Sample size
 //' @param b Mutation count (frequency)
 //' @param e_j Vector of expected coalescence times (length n-1, for j=2 to n)
 // [[Rcpp::export]]
 NumericVector Calculate_SNP_Probabilities(int n, NumericVector b, NumericVector e_j) {
   int m = b.size();
   NumericVector out(m);
   for(int i = 0; i < m; ++i) {
     out[i] = single_snp_prob(n, b[i], e_j);
   }
   return out;
 }
 
 
 
 