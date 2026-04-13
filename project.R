library(Rcpp)
library(ggplot2)

# 1. Compile the C++ code
sourceCpp("calc_snp.cpp")

# 2. Define your population history solver in R (or C++)
# This example assumes a simple constant population where e_j = 2 / (j*(j-1))
solve_ej <- function(n) {
  j <- 2:n
  return(2 / (j * (j - 1)))
}

# 3. Execute
# n_val <- 30
# b_val <- 5
# expected_times <- Calculate_e_j(n_val, 1, kappa = 1)

# result <- Calculate_SNP_Probabilities(n_val, b_val, expected_times)
# print(paste("Probability:", result))

n_val <- 30
b_list <- c(1:15)

# expected_times <- Calculate_e_j(n_val, 3)
# print(expected_times)
# res = Calculate_SNP_Probabilities(n_val, b_list, expected_times)
# print(res)

# Calculate e_j for three different models (example types)
# 1: Constant, 2: Exponential Growth, 3: Population Bottleneck (Hypothetical)
ej_0  <- Calculate_e_j(n_val, population_type = 0, Ne0 = 10)
ej_1  <- Calculate_e_j(n_val, population_type = 1, kappa = 1)
ej_10 <- Calculate_e_j(n_val, population_type = 1, kappa = 10)
ej_2  <- Calculate_e_j(n_val, population_type = 2, Ne0 = 10, NeT = 1, T = 1)
ej_3  <- Calculate_e_j(n_val, population_type = 3, kappa = 1, beta = 0.5)

# Compute probabilities using your vectorized C++ function
prob_0  <- Calculate_SNP_Probabilities(n_val, b_list, ej_0)
prob_1  <- Calculate_SNP_Probabilities(n_val, b_list, ej_1)
prob_10 <- Calculate_SNP_Probabilities(n_val, b_list, ej_10)
prob_2  <- Calculate_SNP_Probabilities(n_val, b_list, ej_2)
prob_3  <- Calculate_SNP_Probabilities(n_val, b_list, ej_3)

# Combine into a data frame for plotting
plot_data <- data.frame(
  b = rep(b_list, 5),
  probability = c(prob_0, prob_1, prob_10, prob_2, prob_3),
  model = rep(c("constant", "exp(1)", "exp(10)", "biphasic", "power"), each = length(b_list))
)

ggplot(plot_data, aes(x = b, y = probability, color = model, shape = model)) +
  geom_point(size = 3) +
  geom_line(alpha = 0.3) +                # Optional: adds faint lines to see trends
  scale_y_log10() +                       # Sets y-axis to logarithmic
  labs(
    title = "SNP Probabilities by Population Model",
    x = "Mutation Count (b)",
    y = "Probability (Log Scale)",
    color = "Population Model",
    shape = "Population Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
