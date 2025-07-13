# GRAND: Graph Release with Assured Node Differential Privacy

## Overview

GRAND (Graph Release with Assured Node Differential Privacy) is an R package that implements a novel method for privatizing network data using differential privacy. The package provides functions for generating synthetic networks, applying differential privacy to network latent positions, and evaluating the utility of privatized networks through various network statistics.

## Installation

### From GitHub (Development Version)

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install GRAND from GitHub
devtools::install_github("lsq0000/GRAND")
```

### Dependencies

The package requires the following R packages:
- EnvStats
- rmutil  
- RSpectra
- diffpriv
- truncnorm
- randnet
- igraph
- HCD
- transport

## Quick Start

```r
library(GRAND)

# Generate a sample network using Latent Space Model
# Note: Use larger networks (n >= 1000) for better stability
network <- LSM.Gen(n = 1000, k = 3, K = 5)

# Privatize the first 500 nodes with different privacy budgets
result <- GRAND.privatize(
  A = network$A, 
  K = 3, 
  idx = 1:500, 
  eps = c(1, 2, 5, 10), 
  model = "LSM"
)

# Evaluate the privatization results
evaluation <- GRAND.evaluate(result)
print(evaluation)

# Evaluate specific statistics only
evaluation_subset <- GRAND.evaluate(result, statistics = c("degree", "triangle"))
print(evaluation_subset)
```

## Main Functions

### `LSM.Gen(n, k, K, avg.d = NULL)`
Generates a random network following the Latent Space Model (LSM) with specified parameters.

- `n`: Integer. Number of nodes in the network.
- `k`: Integer. Dimension of the latent space.
- `K`: Integer. Number of communities/groups.
- `avg.d`: Numeric. Target average degree. If NULL, no degree adjustment is performed.

**Returns**: A list containing:
- `A`: Adjacency matrix of the generated network
- `g`: Graph object of the generated network
- `P`: Probability matrix of the generated network
- `alpha`: Node-specific intercept parameters
- `Z`: Latent positions in k-dimensional space
- `idx`: Community assignments for each node

### `GRAND.privatize(A, K, idx, eps = 1, oracle.dt = NULL, model = "LSM", niter = 500, rho = 0.05)`
Applies the GRAND (Graph Release with Assured Node Differential Privacy) method to privatize network data using differential privacy.

- `A`: Matrix. Adjacency matrix of the input network.
- `K`: Integer. Dimension of the latent space for network embedding.
- `idx`: Integer vector. Indices of nodes to be privatized.
- `eps`: Numeric or vector. Privacy budget parameter(s) for differential privacy. Default is 1.
- `oracle.dt`: List. Optional oracle data containing true parameters for comparison. Default is NULL.
- `model`: Character. Model type, either "LSM" (Latent Space Model) or "RDPG" (Random Dot Product Graph). Default is "LSM".
- `niter`: Integer. Number of iterations for the optimization algorithm. Default is 500.
- `rho`: Numeric. Parameter controlling the neighborhood size for conditional distributions. Default is 0.05.

**Returns**: A list containing:
- `non.private.result`: Results without privacy (original and estimated data)
- `GRAND.result`: Results from GRAND privatization method
- `Laplace.result`: Results from baseline Laplace mechanism
- `eps`: Privacy parameters used
- `oracle.result`: Oracle comparison results (if oracle.dt provided)

### `GRAND.evaluate(result, statistics = c("degree", "triangle", "vshape", "eigen", "harmonic"))`
Evaluates the quality of GRAND privatization results by comparing various network statistics between the original and privatized networks using Wasserstein distance.

- `result`: List. Output from GRAND.privatize function containing privatization results.
- `statistics`: Character vector. Network statistics to evaluate. Options include: "degree", "triangle", "vshape", "eigen", "harmonic". Default is all statistics.

**Returns**: A data frame containing evaluation results with columns:
- `metric`: Type of network statistic evaluated
- `eps`: Privacy parameter used
- `Hat`: Wasserstein distance for non-private estimation
- `Hat2`: Wasserstein distance for holdout set estimation
- `GRAND`: Wasserstein distance for GRAND privatization
- `Laplace`: Wasserstein distance for Laplace mechanism
- `Oracle`: Wasserstein distance for oracle method (if available)

## Features

- **Network Generation**: Generate synthetic networks using Latent Space Models
- **Differential Privacy**: Apply node-level differential privacy to network data
- **Multiple Models**: Support for both LSM and RDPG models
- **Comprehensive Evaluation**: Evaluate utility through multiple network statistics
- **Flexible Privacy Budgets**: Support for multiple privacy levels in a single run

## Methodology

GRAND uses a two-step approach:
1. **Latent Position Estimation**: Estimates latent positions from the network structure
2. **Multivariate Differential Privacy**: Applies DIP (Differential Privacy) mechanism to protect latent positions while preserving network utility

## License

GPL (>= 3)

## Citation

If you use GRAND in your research, please cite:

```
@misc{grand2024,
  title={GRAND: Graph Release with Assured Node Differential Privacy},
  author={Your Name},
  year={2024},
  howpublished={R package}
}
```

## Issues and Contributions

Please report issues or contribute to the package through GitHub.