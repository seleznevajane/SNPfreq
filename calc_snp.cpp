#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/gamma.hpp>

// [[Rcpp::depends(BH)]]

using namespace Rcpp;

//' Calculate expectations of coalescence times
 //' @param n Sample size
 //' @param population_type Model of populations size (0 - constant, 1 - exponential)
 // [[Rcpp::export]]
 NumericVector Calculate_e_j(int n, 
                             int population_type = 0, 
                             int Ne0 = 1, 
                             int K = 10, 
                             double kappa = 1.0,
                             double epsilon = 0.01) {
   NumericVector e_j(n - 1, 0.0);
   
   if (population_type == 0) {
     for (int j = 2; j <= n; j++) {
       e_j[j - 2] = Ne0 * 2.0 / (j * (j - 1.0));
     }
   }
   
   if (population_type == 1) {
     for (int j = 2; j <= n; j++) {
       double x = (j * (j - 1.0)) / (2.0 * kappa);
       e_j[j - 2] = - std::exp(x) * boost::math::expint(-x);
     }
   }
   
   return e_j;
 }
 
 double single_snp_prob(int n, int b, NumericVector e_j){
   // Note: R vectors are 0-indexed. e_j[0] corresponds to j=2.
   // int num_coeffs = n - 1; 
   std::vector<double> W(n + 1, 0.0);
   std::vector<double> V(n + 1, 0.0);
   
   // --- 1. Compute W_j coefficients using Recursion (Eq 15 style) ---
   // Initializing base cases for the recursion
   W[2] = 6.0 / (n + 1.0);
   W[3] = 30.0 * (n - 2.0 * b) / ((n + 1.0) * (n + 2.0));
   
   for (int j = 2; j <= (n - 2); ++j) {
     // Placeholder for Equation 15: 
     // Typically: W_{j+2} = -( (coeff1 * W_j) + (coeff2 * W_{j+1}) ) / coeff3
     // Replace this logic with your specific recursive formula:
     // double j_val = static_cast<double>(j);
     // W[j + 2] = ( (j_val + 2.0) * (n - j_val) * W[j] ) / ((j_val + 1.0) * (n + j_val + 1.0)); 
     W[j + 2] = - W[j] * (1.0 + j) * (3.0 + 2.0 * j) * (n - j + 0.0) /
       (j * (2.0 * j - 1) * (n + j + 1.0)) +
         W[j + 1] * (3.0 + 2.0 * j) * (n - 2.0 * b) /
           (j * (n + j + 1.0));
   }
   
   // --- 2. Compute Denominator coefficients (V_j) (Eq 12 style) ---
   // Initializing base case for the recursion
   V[2] = 6.0 * (n - 1.0) / (n + 1.0);
   
   for (int j = 4; j <= n; j += 2) {
     // Placeholder for Equation 12:
     // V[j] = j * (j - 1.0); // Example: Kingman's coalescent weights
     V[j] = (2.0 * j - 1.0) * (n - j + 1.0) * (n - j + 2.0) * V[j-2] /
       ((2.0 * j - 5.0) * (n + j - 1.0) * (n + j - 2.0));
   }
   
   // --- 3. Final Calculation ---
   double numerator = 0.0;
   double denominator = 0.0;
   
   for (int j = 2; j <= n; ++j) {
     // Map j index to e_j index (j=2 -> index 0)
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
 
 
 
 