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
Generates a random network following the Latent Space Model (LSM).

- `n`: Number of nodes
- `k`: Dimension of latent space  
- `K`: Number of communities
- `avg.d`: Target average degree (optional)

### `GRAND.privatize(A, K, idx, eps = 1, oracle.dt = NULL, model = "LSM", niter = 500, rho = 0.05)`
Applies GRAND privatization to network data.

- `A`: Adjacency matrix
- `K`: Latent space dimension
- `idx`: Indices of nodes to privatize
- `eps`: Privacy budget parameter(s)
- `model`: "LSM" or "RDPG"

### `GRAND.evaluate(result, statistics = c("degree", "triangle", "vshape", "eigen", "harmonic"))`
Evaluates privatization quality using various network statistics.

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