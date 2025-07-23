\name{GRAND-package}
\alias{GRAND-package}
\alias{_PACKAGE}
\docType{package}
\title{
\packageTitle{GRAND}
}
\description{
\packageDescription{GRAND}
}
\details{
GRAND (Graph Release with Assured Node Differential Privacy) is an R package that implements a novel method for privatizing network data using differential privacy. The package provides functions for generating synthetic networks, applying differential privacy to network latent positions, and evaluating the utility of privatized networks through various network statistics.

The package implements a two-step approach:
\enumerate{
  \item \strong{Latent Position Estimation}: Estimates latent positions from the network structure using either Latent Space Models (LSM) or Random Dot Product Graphs (RDPG)
  \item \strong{Multivariate Differential Privacy}: Applies the DIP (Differential Privacy) mechanism to protect latent positions while preserving network utility
}

\strong{Main Functions:}
\itemize{
  \item \code{\link{LSM.Gen}}: Generate synthetic networks using Latent Space Models
  \item \code{\link{GRAND.privatize}}: Apply GRAND privatization to network data
  \item \code{\link{GRAND.evaluate}}: Evaluate the utility of privatized networks
}

\strong{Key Features:}
\itemize{
  \item Network Generation: Generate synthetic networks using Latent Space Models
  \item Differential Privacy: Apply node-level differential privacy to network data
  \item Multiple Models: Support for both LSM and RDPG models
  \item Comprehensive Evaluation: Evaluate utility through multiple network statistics
  \item Flexible Privacy Budgets: Support for multiple privacy levels in a single run
}

The package is designed for researchers and practitioners working with sensitive network data who need to balance privacy protection with data utility.
}
\author{
\packageAuthor{GRAND}

Maintainer: \packageMaintainer{GRAND}
}
\references{
P. D. Hoff, A. E. Raftery, and M. S. Handcock. Latent space approaches to social network analysis. 
Journal of the American Statistical Association, 97(460):1090–1098, 2002.

S. J. Young and E. R. Scheinerman. Random dot product graph models for social networks. 
In International Workshop on Algorithms and Models for the Web-Graph, pages 138–149. Springer, 2007.

S. Liu, X. Bi, and T. Li. GRAND: Graph Release with Assured Node Differential Privacy. 
arXiv preprint arXiv:2507.00402, 2025.
}
\examples{
# Generate a sample network
dt <- LSM.Gen(n = 500, k = 2, K = 3)

# Privatize the first 250 nodes with multiple epsilon values
result <- GRAND.privatize(A = dt$A, K = 2, idx = 1:250, eps = c(1, 2, 5, 10))

# Evaluate the privatization results
evaluation <- GRAND.evaluate(dt, result)
print(evaluation)
}
\keyword{package}