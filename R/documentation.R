#' SNP Probabilities
#'
#' @description
#' Computes mutation frequency probabilities for a given set of mutation counts 
#' using an optimized recursive coefficient method.
#'
#' @param n Integer sample size.
#' @param b Scalar integer or integer vector of mutation counts (values between 1 and n-1).
#' @param e Numeric vector of expected coalescence times calculated using [coaltimes()].
#'
#' @details
#' We consider the process of coalescence with time-changing effective population size.
#' For a sample size of \eqn{n} DNA sequences, let \eqn{T_k}, \eqn{k = 2, 3, \dots, n}, 
#' be cumulative times to coalescence from a sample of size \eqn{n} to a sample of 
#' size \eqn{k - 1}.
#'
#' Under the infinite-sites mutation model, the probability \eqn{q_{nb}} that a 
#' SNP site in a sample of \eqn{n} chromosomes has \eqn{b} mutant bases 
#' (\eqn{0 < b < n}) is given by Polanski & Kimmel (2003):
#' 
#' \deqn{q_{nb} = \frac{\sum_{j=2}^n e_j W_{bj}^n}{\sum_{j=2}^n e_j V_j^n}}
#' 
#' where the coefficients \eqn{V_j^n} and \eqn{W_{bj}^n} are computed recursively 
#' and are independent of the population growth model.
#'
#' The values \eqn{e_j} represent expectations of coalescence times under the 
#' effective population size history \eqn{N_e(t), t \in [0, \infty)}. To compute 
#' \eqn{e_j}, use the function [coaltimes()].
#'
#' @return A numeric vector of probabilities.
#'
#' @references
#' Polanski, A., & Kimmel, M. (2003). New explicit expressions for relative 
#' frequencies of single-nucleotide polymorphisms with application to 
#' statistical inference on population growth. Genetics, 165(1), 427-436.
#' <https://pmc.ncbi.nlm.nih.gov/articles/PMC1462751/>
#'
#' @seealso [coaltimes()]
#'
#' @examples
#' library(SNPprob)
#' 
#' n = 30
#' 
#' # Constant population size model
#' e <- coaltimes(n, population_type = 0)
#' 
#' # Single mutation count
#' prob <- SNPprob(n, 10, e)
#' 
#' # Range of mutation counts
#' probs <- SNPprob(n, c(1:15), e)
#' 
#' @name SNPprob
NULL

#' Expectations of Coalescence Times
#'
#' @description
#' Computes a vector of expected coalescence times up to a multiplier for 
#' various population growth models.
#'
#' @usage 
#' coaltimes(n, population_type = 0)
#' coaltimes(n, population_type = 1, kappa)
#' coaltimes(n, population_type = 2, Ne0, NeT, T)
#' coaltimes(n, population_type = 3, kappa, beta)
#' coaltimes(n, population_type = 4, tk, eta)
#'
#' @param n Sample size.
#' @param population_type Integer code for the population model:
#' * `0`: Constant.
#' * `1`: Exponential (requires `kappa`).
#' * `2`: Biphasic (requires `Ne0`, `NeT`, and `T`).
#' * `3`: Power (requires `kappa` and `beta`).
#' * `4`: Stepwise Constant (requires `tk` and `eta`).
#' @param kappa Product of current effective population size and growth rate, 
#' `kappa = Ne0 * r` (Exponential, Power models).
#' @param Ne0 Current effective population size (Biphasic model).
#' @param T Time in generations when the model changed (Biphasic model).
#' @param NeT Effective population size at moment T (Biphasic model).
#' @param beta Growth equation power (Power model). If `beta < 0`, increasing 
#' population; if `0 < beta < 1`, decreasing population.
#' @param tk Vector of times of population changes (Stepwise Constant model). See details.
#' @param eta Vector of effective population sizes for intervals in `tk` 
#' (Stepwise Constant model). See details.
#'
#' @details
#' We consider the process of coalescence with time-changing effective population 
#' size \eqn{N_e(t)}. For a sample size \eqn{n}, the function computes expectations 
#' of coalescence times:
#' \deqn{e_j = \int_0^\infty t \, q_j(t) \, dt}
#' where the distribution \eqn{q_j(t)} is defined as:
#' \deqn{q_j(t) = \frac{\binom{j}{2}}{N_e(t)} \exp\left( -\int_0^t \frac{\binom{j}{2} d\sigma}{N_e(\sigma)} \right)}
#'
#' The following population models are supported:
#' * **Constant:** \eqn{N_e(t) = N_{e0}}
#' * **Exponential:** \eqn{N_e(t) = N_{e0} \exp(-rt)}
#' * **Biphasic (constant, then exponential):**
#'   \deqn{N_e(t) = \begin{cases} N_{e0} \exp(-rt), & 0 \le t \le T \\ N_{e0} \exp(-rT), & t > T \end{cases}}
#' * **Power:** \eqn{N_e(t) = N_{e0}(1+rt)^\beta}
#' * **Stepwise Constant:** #'   \deqn{N_e(t) = \begin{cases} \eta_1, \quad & 0 \le t < t_1 \\
#'                                \eta_2, \quad & t_1 \le t < t_2 \\ \vdots \\
#'                                \eta_K = 0, \quad & t \ge t_{K-1} \end{cases}}
#' 
#' For stepwise constant model, you must provide vectors `tk` = \eqn{(0, t_1,\dots, t_{K-1})} 
#' and `eta` = \eqn{(\eta_1,\dots, \eta_K)}.
#'
#' This function is designed to be used in conjunction with [SNPprob()] 
#' and computes expected times only up to a multiplier.
#'
#' @return A numeric vector of expected coalescence times.
#'
#' @references
#' Polanski, A., & Kimmel, M. (2003). New explicit expressions for relative 
#' frequencies of single-nucleotide polymorphisms with application to 
#' statistical inference on population growth. Genetics, 165(1), 427-436.
#' <https://pmc.ncbi.nlm.nih.gov/articles/PMC1462751/>
#'
#' @seealso [SNPprob()]
#'
#' @examples
#' library(SNPprob)
#' n <- 30
#'
#' # Constant population model
#' e_const <- coaltimes(n, population_type = 0)
#'
#' # Exponential population model
#' e_exp <- coaltimes(n, population_type = 1, kappa = 10)
#'
#' # Biphasic population model
#' e_biph <- coaltimes(n, population_type = 2, Ne0 = 10, NeT = 1, T = 1)
#' 
#' # Power population model
#' e_power <- coaltimes(n, population_type = 3, kappa = 1, beta = 0.5)
#'
#' # Stepwise constant population model
#' e_step <- coaltimes(n, population_type = 4, 
#'                     tk = c(0, 9, 100, 600, 10000, 290000), 
#'                     eta = c(1e8, 1e6, 1e5, 4000, 5000, 0))
#' 
#' @name coaltimes
NULL

#' @useDynLib SNPprob, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL