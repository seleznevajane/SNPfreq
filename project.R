library(Rcpp)

# 1. Compile the C++ code
sourceCpp("calc_snp.cpp")

# 2. Define your population history solver in R (or C++)
# This example assumes a simple constant population where e_j = 2 / (j*(j-1))
solve_ej <- function(n) {
  j <- 2:n
  return(2 / (j * (j - 1)))
}

# 3. Execute
n_val <- 100
b_val <- 5
expected_times <- Calculate_e_j(n_val, 0) # solve_ej(n)

result <- Calculate_SNP_Probabilities(n_val, b_val, expected_times)
print(paste("Probability:", result))