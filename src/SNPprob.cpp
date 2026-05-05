#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/gamma.hpp>

// [[Rcpp::depends(BH)]]

// Calculate scaled exponential integral using asymptotic expansion for large x:
// exp(x) * Ei(-x) = 1/x * (1 - 1/x + 2/x^2 - 6/x^3 + 24/x^4)
double scaled_expint(double x) {
  if (x > 50.0) {
    double invX = 1.0 / x;
    return invX * (1.0 - invX + (2.0 * invX * invX) - (6.0 * invX * invX * invX) +
                   (24.0 * invX * invX * invX * invX));
  } else {
    return -std::exp(x) * boost::math::expint(-x);
  }
}

// Caclulate expected coalescence times
//' @rdname coaltimes
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector coaltimes(int n, 
                              int population_type = 0, 
                              double kappa = 1.0,
                              int Ne0 = 100,
                              int T = 10, 
                              int NeT = 10,
                              double beta = 0.5,
                              Rcpp::NumericVector tk = Rcpp::NumericVector::create(),
                              Rcpp::NumericVector eta = Rcpp::NumericVector::create()) {
  Rcpp::NumericVector e_j(n - 1, 0.0);
  
  // --- Constant population model ---
  if (population_type == 0 || 
      (population_type == 1 && kappa == 0) ||
      (population_type == 3 && beta == 0)) {
    for (int j = 2; j <= n; j++) {
      e_j[j - 2] = 2.0 / (j * (j - 1.0));
    }
  }
   
  // --- Exponential population model ---
  if (population_type == 1) {
    if (kappa < 0) {
      Rcpp::stop("Exponential model (type 1) requires 'kappa' to be a positive real number.");
    }
    double x;
    for (int j = 2; j <= n; j++) {
      x = (j * (j - 1.0)) / (2.0 * kappa);
      e_j[j - 2] = scaled_expint(x);
    }
  }
  
  // --- Biphasic population model ---
  if (population_type == 2) {
    if (T <= 0) {
      Rcpp::stop("Biphasic model (type 2) requires 'T' to be a positive integer number.");
    }
    if (Ne0 <= 0) {
      Rcpp::stop("Biphasic model (type 2) requires 'Ne0' to be a positive integer number.");
    }
    if (NeT <= 0) {
      Rcpp::stop("Biphasic model (type 2) requires 'NeT' to be a positive integer number.");
    }
    if (Ne0 <= NeT) {
      Rcpp::stop("Biphasic model (type 2) requires 'NeT' to be less than 'Ne0' as the second phase is exponential growth.");
    }
    double r = std::log(Ne0 * 1.0 / NeT) / (T * 1.0), x0, xT;
    for (int j = 2; j <= n; j++) {
      x0 = (j * (j - 1.0)) / (2.0 * r * Ne0);
      xT = x0 * std::exp(r * T);
      e_j[j - 2] = std::exp(x0) * (boost::math::expint(-xT) - boost::math::expint(-x0)) +
                   std::exp(x0 - xT) / (r * xT);
    }
  }
   
  // --- Power population model ---
  if (population_type == 3) {
    if (kappa <= 0) {
      Rcpp::stop("Power model (type 3) requires 'kappa' to be a positive real number.");
    }
    if (beta >= 1) {
      Rcpp::stop("Power model (type 3) requires 'beta' to be less than 1.");
    }
    double C;
    for (int j = 2; j <= n; j++) {
      C = (j * (j - 1.0)) / (2.0 * kappa * (1.0 - beta));
      e_j[j - 2] = std::exp(C + std::log(boost::math::tgamma(1.0 / (1.0 - beta), C)) -
                   std::log((1.0 - beta) * std::pow(C, 1.0 / (1.0 - beta))));
    }
  }
  
  // --- Stepwise constant population model ---
  if (population_type == 4) {
    if (tk.size() == 0 || eta.size() == 0) {
      Rcpp::stop("Stepwise constant model (type 4) requires 'tk' and 'eta' vectors.");
    }
    if (tk[0] != 0) {
      Rcpp::stop("Stepwise constant model (type 4) requires 'tk[0]' to be 0.");
    }
    int K = tk.size();
    std::vector <double> tau(K, 0.0);
    for (int k = 1; k < K - 2; k++) {
      tau[k] = tau[k - 1] + 1.0 * (tk[k] - tk[k - 1]) / eta[k];
    }
    double C;
    for (int j = 2; j <= n; j++) {
      C = (j * (j - 1.0)) / 2.0;
      for (int k = 1; k <= K - 1; k++) {
        e_j[j - 2] += eta[k] * std::exp(-C * tau[k - 1]) * 
          (1.0 - std::exp(- C * (tk[k] - tk[k - 1]) / eta[k])) / C;
      }
    }
  }
   
  return e_j;
}

// Calculate SNP probability for a single number of mutant bases
double single_snp_prob(int n, 
                       int b, 
                       Rcpp::NumericVector e,
                       double* W, 
                       double* V) {
   
  std::fill(W, W + n + 1, 0.0);
  std::fill(V, V + n + 1, 0.0);
   
  // --- Compute W_j coefficients---
  W[2] = 6.0 / (n + 1.0);
  W[3] = 30.0 * (n - 2.0 * b) / ((n + 1.0) * (n + 2.0));
   
  for (int j = 2; j <= (n - 2); ++j) {
    W[j + 2] = - W[j] * (1.0 + j) * (3.0 + 2.0 * j) * (n - j + 0.0) /
               (j * (2.0 * j - 1) * (n + j + 1.0)) +
               W[j + 1] * (3.0 + 2.0 * j) * (n - 2.0 * b) / (j * (n + j + 1.0));
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
    numerator += e[j - 2] * W[j];
    denominator += e[j - 2] * V[j];
  }
   
  return numerator / denominator;
}
 
// Calculate SNP probabilities for a range of mutant bases
//' @rdname SNPprob
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector SNPprob(int n, 
                            Rcpp::IntegerVector b, 
                            Rcpp::NumericVector e) {
  // --- Catch incorrect input ---
  if (e.size() != n - 1) {
    Rcpp::stop("Expected coalescence times e must be a vector of size n - 1 (use SNPprob::coaltimes)");
  }
  if (!std::all_of(b.begin(), b.end(), [n] (double val) 
                  {return val >= 1 && val <= (n - 1);})) {
    Rcpp::stop("Mutation counts b must be integer numbers from the interval [1, n - 1]");
  }
     
  // --- Allocate memory for results and recursive coefficients ---
  int m = b.size();
  Rcpp::NumericVector out(m);
  
  std::vector<double> W(n + 1);
  std::vector<double> V(n + 1);
   
  double* pW = W.data();
  double* pV = V.data();
   
  // --- Calculate probabilities for each mutation count ---
  for(int i = 0; i < m; ++i) {
    out[i] = single_snp_prob(n, b[i], e, pW, pV);
  }
   
  return out;
}
 
 
 
 