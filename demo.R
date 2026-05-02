library(Rcpp)
library(ggplot2)
library(scales)

# Compile the C++ code
sourceCpp("SNPfreq.cpp")


# Calculate a single SNP probability
n_val <- 30
b_val <- 5
expected_times <- coaltimes(n_val, 4, tk = c(0, 10, 100), eta = c(100, 20, 50))

result <- SNPprob(n_val, b_val, expected_times)
print(paste("Probability:", result))

# Calculate SNP probabilities for all population models
n_val <- 30
b_list <- c(1:15)

# Calculate e_j for different models
ej_0  <- coaltimes(n_val, population_type = 0, Ne0 = 10)
ej_1  <- coaltimes(n_val, population_type = 1, kappa = 10)
ej_2  <- coaltimes(n_val, population_type = 2, Ne0 = 10, NeT = 1, T = 1)
ej_3  <- coaltimes(n_val, population_type = 3, kappa = 1, beta = 0.5)
ej_4  <- coaltimes(n_val, population_type = 4, 
                       tk = c(0, 9, 100, 600, 10000, 290000), 
                       eta = c(100000000, 1000000, 100000, 4000, 5000, 0))

# Compute probabilities
prob_0  <- SNPprob(n_val, b_list, ej_0)
prob_1  <- SNPprob(n_val, b_list, ej_1)
prob_2  <- SNPprob(n_val, b_list, ej_2)
prob_3  <- SNPprob(n_val, b_list, ej_3)
prob_4  <- SNPprob(n_val, b_list, ej_4)

# Visualize the SNP probabilities for different models
plot_data <- data.frame(
  b = rep(b_list, 5),
  probability = c(prob_0, prob_1, prob_2, prob_3, prob_4),
  model = factor(rep(c("Constant", "Exponential", "Biphasic", "Power", "Stepwise"), 
                     each = length(b_list)),
                 levels = c("Constant", "Exponential", "Biphasic", "Power", "Stepwise"))
)

ggplot(plot_data, aes(x = b, y = probability, color = model, shape = model)) +
  geom_line(linewidth = 0.5, alpha = 0.5) +
  geom_point(size = 3, stroke = 0.8) + 
  
  scale_y_log10(
    labels = trans_format("log10", math_format(10^.x)),
    breaks = trans_breaks("log10", function(x) 10^x)
  ) +
  annotation_logticks(sides = "l", color = "grey80") +
  
  scale_color_viridis_d(option = "plasma", end = 0.8) +
  
  labs(
    title = "SNP Probabilities for All Population Growth Model",
    subtitle = "Results for sample size n = 30",
    x = "Mutation Count (b)",
    y = "Probability q(n, b)",
    color = "Model Type:",
    shape = "Model Type:"
  ) +
  
  theme_minimal(base_size = 10) + 
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(), 
    legend.background = element_rect(fill = "white", color = NA)
  )
