# SNPfreq: SNP Probabilities for Various Population Models

`SNPfreq` is an R package for the efficient and precise analytical calculation 
of mutation frequency spectra (SNP probabilities) for various demographic models. 

This package implements the recursive coefficient method described by 
[**Polanski and Kimmel (2003)**](https://doi.org/10.1093/genetics/165.1.427). 
It provides a significant speedup over Monte Carlo simulations for calculating 
expected coalescence times and mutation frequency probabilities under constant, 
exponential, biphasic, power-law, and stepwise constant population growth models.

## Theoretical Overview

We consider the process of coalescence with time changing effective population size.
Sample size is $n$ DNA sequences, and $T_k$, $k = 2, 3,..., n$, are cumulative times 
to coalescence from sample of size $n$ to sample of size $k - 1$. Under the infinite-sites 
mutation model, the probability $q_{nb}$ that a SNP site in a sample of $n$ 
chromosomes has $b$ mutant bases ($0 < b < n$) is given by 
[Polanski and Kimmel (2003)](https://doi.org/10.1093/genetics/165.1.427):

$$q_{nb} = \dfrac{\sum_{j=2}^n e_jW_{bj}^n}{\sum_{j=2}^ne_jV_j^n}$$

where coefficients $V_j^n$ and $W_{bj}^n$ are computed recursively 
and are independent of the population growth model. And

$$e_j = \int_0^\infty t\,q_j(t)\,dt,$$

are expectations of coalescence times distributed as

$$q_j(t) = \frac{\dbinom{j}{2}}{N_e(t)} \exp\left( -\int_0^t \frac{\dbinom{j}{2}d\sigma}{N_e(\sigma)} \right),$$

with the effective population size history described as a function of reverse time, $N_e(t)$, $t\in[0,\infty)$.

The following population models are supported:
 * **Constant:**

$$N_e(t) = N_{e0}$$
 
 * **Exponential:**

$$N_e(t) = N_{e0} \exp(-rt)$$
   
 * **Biphasic (constant, then exponential):**

$$N_e(t) = \begin{cases} N_{e0} \exp(-rt), & 0 \le t \le T \\ 
N_{e0} \exp(-rT), & t > T \end{cases}$$

 * **Power:**

$$N_e(t) = N_{e0}(1+rt)^\beta$$
   
 * **Stepwise Constant:**   

$$N_e(t) = \begin{cases} \eta_1, \quad & 0 \le t < t_1 \\
\eta_2, \quad & t_1 \le t < t_2 \\ \vdots \\
\eta_K = 0, \quad & t \ge t_{K-1} \end{cases}$$



## Installation

You can install the development version directly from GitHub using the remotes
package:

```R
# Install remotes if you haven't already
install.packages("remotes")
# Install SNPfreq
remotes::install_github("seleznevajane/SNPfreq")
```

## Quick Start / Usage

The core workflow involves calculating expected coalescence times and then deriving
the SNP probabilities.

```R
library(SNPfreq)
# 1. Define sample size
n <- 30
# 2. Calculate expected coalescence times under exponential growth
# (population_type = 1 for exponential)
times <- coaltimes(n, population_type = 1, kappa = 1)
# 3. Calculate probability of finding 1 to 5 mutant bases
probs <- SNPprob(n, b = 1:5, e = times)
# 4. View results
print(probs)
```

## Main Functions

`coaltimes(n, population_type, ...)` : Computes the expected coalescence times $e_j$.

`SNPprob(n, b, e)` : Computes the probability of $b$ mutants in a sample of
size $n$ given expected times $e$.

## References
Polanski, A., & Kimmel, M. (2003). New explicit expressions for relative 
frequencies of single-nucleotide polymorphisms with application to statistical 
inference on population growth. Genetics, 165(1), 427-434. 
[doi:10.1093/genetics/165.1.427](https://doi.org/10.1093/genetics/165.1.427)

## Authors
Evgeniia Selezneva - Maintainer

Thy Luong

Andrzej Polanski - Original Algorithm Design

Marek Kimmel - Original Algorithm Design

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
