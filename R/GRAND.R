#' @importFrom RSpectra eigs_sym
#' @importFrom diffpriv DPParamsEps DPMechLaplace releaseResponse
#' @importFrom truncnorm rtruncnorm
#' @importFrom stats glm binomial lm
#' @importFrom randnet LSM.PGD
#' @importFrom igraph graph.adjacency degree count_triangles eigen_centrality harmonic_centrality
#' @importFrom HCD gen.A.from.P
#' @importFrom transport wasserstein1d

f <- function(X) X

ase <- function(A, K) {
  A2 <- t(A) %*% A
  eig <- eigs_sym(A2, k = K)
  Xhat <- t(t(eig$vectors[, 1:K, drop = FALSE]) * eig$values[1:K]^(1/4))

  return(Xhat)
}

Add.Laplace <- function(X, eps = 1) {
  p <- ncol(X)
  n <- nrow(X)
  pparams <- DPParamsEps(epsilon = eps / p)
  X_lap <- X

  for (j in 1:p) {
    mechanism <- DPMechLaplace(target = f, sensitivity = max(abs(X[, j])), dims = n)
    r_lap <- releaseResponse(mechanism, privacyParams = pparams, X = X[, j])
    X_lap[, j] <- r_lap$response
  }

  return(X_lap)
}

#' Generate Latent Space Model Network
#'
#' @title Generate Latent Space Model Network
#' @description Generates a random network following the Latent Space Model (LSM) with specified parameters.
#' @param n Integer. Number of nodes in the network.
#' @param k Integer. Dimension of the latent space.
#' @param K Integer. Number of communities/groups.
#' @param avg.d Numeric. Target average degree. If NULL, no degree adjustment is performed.
#' @return A list containing:
#' \itemize{
#'   \item A: Adjacency matrix of the generated network
#'   \item g: Community assignments for each node
#'   \item P: Probability matrix
#'   \item alpha: Node-specific intercept parameters
#'   \item Z: Latent positions in k-dimensional space
#' }
#' @export
#' @examples
#' # Generate a network with 100 nodes, 2D latent space, 3 communities
#' result <- LSM.Gen(n = 100, k = 2, K = 3)
#' # Generate with target average degree of 10
#' result <- LSM.Gen(n = 100, k = 2, K = 3, avg.d = 10)
LSM.Gen <- function(n, k, K, avg.d = NULL) {
  alpha <- runif(n, 1, 3)
  alpha <- -alpha / 2

  mu <- matrix(runif(K * k, -1, 1), K, k)
  idx <- sample(1:K, n, replace = TRUE)
  mu <- mu[idx, ]

  Z <- mu + matrix(rtruncnorm(n * k, -2, 2), n, k)
  J <- diag(n) - rep(1, n) %*% t(rep(1, n)) / n
  Z <- J %*% Z
  G <- Z %*% t(Z)
  Z <- Z / sqrt(sqrt(sum(G^2)) / n)

  theta <- alpha %*% t(rep(1, n)) + rep(1, n) %*% t(alpha) + Z %*% t(Z)
  P <- 1 / (1 + exp(-theta))

  if (!is.null(avg.d)) {
    for (i in 1:10) {
      ratio <- mean(rowSums(P)) / avg.d
      beta <- -log(ratio)
      alpha <- alpha + beta / 2
      theta <- alpha %*% t(rep(1, n)) + rep(1, n) %*% t(alpha) + Z %*% t(Z)
      P <- 1 / (1 + exp(-theta))
    }
  }

  upper.index <- which(upper.tri(P))
  upper.p <- P[upper.index]
  upper.u <- runif(length(upper.p))
  upper.A <- rep(0, length(upper.p))
  upper.A[upper.u < upper.p] <- 1
  A <- matrix(0, n, n)
  A[upper.index] <- upper.A
  A <- A + t(A)
  diag(A) <- 0

  return(list(A = A, g = idx, P = P, alpha = alpha, Z = Z))
}

GRAND.one.node <- function(A, given.index, new.index, given.Z, given.alpha = NULL, model = c("LSM", "RDPG")) {
  model <- match.arg(model)

  Y <- A[new.index, given.index]

  if (model == "LSM") {
    if (is.null(given.alpha)) {
      stop("given.alpha must be provided for LSM model.")
    }

    fit <- glm(Y ~ given.Z, family = binomial(), offset = given.alpha, maxit = 1000)
  } else {
    fit <- lm(Y ~ given.Z - 1)
  }

  return(fit)
}

GRAND.estimate <- function(A, K, holdout.index, release.index, model = c("LSM", "RDPG"), niter = 500) {
  model <- match.arg(model)

  m <- length(holdout.index)
  n <- length(release.index)
  A22 <- A[holdout.index, holdout.index]

  if (model == "LSM") {
    fitting.count <- 0

    while (fitting.count < 5) {
      fit.holdout <- LSM.PGD(A22, k = K, niter = niter)
      fitting.count <- fitting.count + 1

      if (sum(is.na(fit.holdout$Z)) == 0) {
        break
      }
    }

    if (sum(colSums(abs(fit.holdout$Z)) < 1e-4) > 0) {
      print(paste0("Potential deficiency of ", sum(colSums(abs(fit.holdout$Z)) < 1e-4), "."))
    }

    print(paste0("PGD used ", length(fit.holdout$obj), " iterations."))

    Z.holdout <- fit.holdout$Z
    alpha.holdout <- fit.holdout$alpha
  } else {
    Z.holdout <- ase(A22, K)
    alpha.holdout <- NULL
  }

  Z.hat <- matrix(0, nrow = m + n, ncol = K)
  alpha.hat <- rep(0, m + n)
  Z.hat[holdout.index, ] <- Z.holdout
  if (!is.null(alpha.holdout)) {
    alpha.hat[holdout.index] <- alpha.holdout
  }

  for (j in 1:n) {
    fit <- GRAND.one.node(A = A,
                          given.index = holdout.index,
                          new.index = release.index[j],
                          given.Z = Z.holdout,
                          given.alpha = alpha.holdout,
                          model = model)

    if (model == "LSM") {
      Z.hat[release.index[j], ] <- fit$coefficients[2:(K + 1)]
      alpha.hat[release.index[j]] <- fit$coefficients[1]
    } else {
      Z.hat[release.index[j], ] <- fit$coefficients[1:K]
    }
  }

  if (model == "LSM") {
    return(list(Z = Z.hat, alpha = alpha.hat))
  } else {
    return(list(Z = Z.hat))
  }
}

#' GRAND: Graph Release with Assured Node Differential Privacy
#'
#' @title GRAND Privatization of Network Data
#' @description Applies the GRAND (Graph Release with Assured Node Differential Privacy) method
#' to privatize network data using differential privacy. The method estimates latent positions
#' from the network and applies multivariate differential privacy to protect sensitive information.
#' @param A Matrix. Adjacency matrix of the input network.
#' @param K Integer. Dimension of the latent space for network embedding.
#' @param idx Integer vector. Indices of nodes to be privatized.
#' @param eps Numeric or vector. Privacy budget parameter(s) for differential privacy. Default is 1.
#' @param oracle.dt List. Optional oracle data containing true parameters for comparison. Default is NULL.
#' @param model Character. Model type, either "LSM" (Latent Space Model) or "RDPG" (Random Dot Product Graph). Default is "LSM".
#' @param niter Integer. Number of iterations for the optimization algorithm. Default is 500.
#' @param rho Numeric. Parameter controlling the neighborhood size for conditional distributions. Default is 0.05.
#' @return A list containing:
#' \itemize{
#'   \item non.private.result: Results without privacy (original and estimated data)
#'   \item GRAND.result: Results from GRAND privatization method
#'   \item Laplace.result: Results from baseline Laplace mechanism
#'   \item eps: Privacy parameters used
#'   \item oracle.result: Oracle comparison results (if oracle.dt provided)
#' }
#' @export
#' @examples
#' # Generate a sample network
#' dt <- LSM.Gen(n = 1000, k = 3, K = 5)
#' # Privatize the first 500 nodes with epsilon = 1, 2, 5, 10
#' result <- GRAND.privatize(A = dt$A, K = 3, idx = 1:500, eps = c(1, 2, 5, 10))
GRAND.privatize <- function(A, K, idx, eps = 1, oracle.dt = NULL, model = c("LSM", "RDPG"), niter = 500, rho = 0.05) {
  model <- match.arg(model)

  n <- nrow(A)
  n1 <- length(idx)
  n2 <- n - n1
  holdout.idx <- setdiff(1:n, idx)
  A11 <- A[idx, idx]
  g1 <- graph.adjacency(A11, "undirected")
  A22 <- A[holdout.idx, holdout.idx]
  g2 <- graph.adjacency(A22, "undirected")

  fit <- GRAND.estimate(A = A,
                        K = K,
                        holdout.index = holdout.idx,
                        release.index = idx,
                        model = model,
                        niter = niter)
  X.trans <- if (model == "LSM") cbind(fit$alpha, fit$Z) else fit$Z

  X1.hat <- X.trans[idx, , drop = FALSE]
  X2.hat <- X.trans[holdout.idx, , drop = FALSE]

  if (model == "LSM") {
    alpha1.hat <- X1.hat[, 1, drop = FALSE]
    Z1.hat <- X1.hat[, -1, drop = FALSE]
    theta1.hat <- alpha1.hat %*% t(rep(1, n1)) + rep(1, n1) %*% t(alpha1.hat) + Z1.hat %*% t(Z1.hat)
    P1.hat <- 1 / (1 + exp(-theta1.hat))
  } else {
    Z1.hat <- X1.hat
    theta1.hat <- Z1.hat %*% t(Z1.hat)
    P1.hat <- pmin(pmax(theta1.hat, 1e-5), 1 - 1e-5)
  }
  A1.hat <- gen.A.from.P(P1.hat)
  g1.hat <- graph.adjacency(A1.hat, "undirected")

  if (model == "LSM") {
    alpha2.hat <- X2.hat[, 1, drop = FALSE]
    Z2.hat <- X2.hat[, -1, drop = FALSE]
    theta2.hat <- alpha2.hat %*% t(rep(1, n2)) + rep(1, n2) %*% t(alpha2.hat) + Z2.hat %*% t(Z2.hat)
    P2.hat <- 1 / (1 + exp(-theta2.hat))
  } else {
    Z2.hat <- X2.hat
    theta2.hat <- Z2.hat %*% t(Z2.hat)
    P2.hat <- pmin(pmax(theta2.hat, 1e-5), 1 - 1e-5)
  }
  A2.hat <- gen.A.from.P(P2.hat)
  g2.hat <- graph.adjacency(A2.hat, "undirected")

  non.private.result <- list(A1 = A11, A2 = A22,
                             A1.hat = A1.hat, A2.hat = A2.hat,
                             P1.hat = P1.hat, P2.hat = P2.hat,
                             g1 = g1, g2 = g2,
                             g1.hat = g1.hat, g2.hat = g2.hat,
                             X1.hat = X1.hat, X2.hat = X2.hat)
  L <- length(eps)
  GRAND.result <- list()
  oracle.result <- list()

  if (!is.null(oracle.dt)) {
    if (model == "LSM") {
      oracle.X1 <- cbind(oracle.dt$alpha[idx], oracle.dt$Z[idx, , drop = FALSE])
      oracle.X2 <- cbind(oracle.dt$alpha[holdout.idx], oracle.dt$Z[holdout.idx, , drop = FALSE])
    } else {
      ASE.oracle <- ase(oracle.dt$P, K)
      oracle.X1 <- ASE.oracle[idx, , drop = FALSE]
      oracle.X2 <- ASE.oracle[holdout.idx, , drop = FALSE]
    }
  }

  for (jj in 1:L) {
    cat(paste("Calling GRAND with \u03B5=", eps[jj], ".\n", sep = ""))
    X1.dip <- DIP.multivariate(X1.hat, eps[jj], X2.hat, rho)
    cat("Finish GRAND.\n")

    if (model == "LSM") {
      alpha1.dip <- X1.dip[, 1, drop = FALSE]
      Z1.dip <- X1.dip[, -1, drop = FALSE]
      theta1.dip <- alpha1.dip %*% t(rep(1, n1)) + rep(1, n1) %*% t(alpha1.dip) + Z1.dip %*% t(Z1.dip)
      P1.dip <- 1 / (1 + exp(-theta1.dip))
    } else {
      Z1.dip <- X1.dip
      theta1.dip <- Z1.dip %*% t(Z1.dip)
      P1.dip <- pmin(pmax(theta1.dip, 1e-5), 1 - 1e-5)
    }
    A1.dip <- gen.A.from.P(P1.dip)
    g1.dip <- graph.adjacency(A1.dip, "undirected")
    GRAND.result[[jj]] <- list(A1.grand = A1.dip,
                               P1.grand = P1.dip,
                               g1.grand = g1.dip,
                               X1.grand = X1.dip)

    if (!is.null(oracle.dt)) {
      cat("Call the oracle version.\n")
      cat(paste("Calling GRAND with \u03B5=", eps[jj], ".\n", sep = ""))
      oracle.X1.dip <- DIP.multivariate(oracle.X1, eps[jj], oracle.X2, rho)
      cat("Finish GRAND.\n")

      if (model == "LSM") {
        oracle.alpha1.dip <- oracle.X1.dip[, 1, drop = FALSE]
        oracle.Z1.dip <- oracle.X1.dip[, -1, drop = FALSE]
        oracle.theta1.dip <- oracle.alpha1.dip %*% t(rep(1, n1)) + rep(1, n1) %*% t(oracle.alpha1.dip) + oracle.Z1.dip %*% t(oracle.Z1.dip)
        oracle.P1.dip <- 1 / (1 + exp(-oracle.theta1.dip))
      } else {
        oracle.Z1.dip <- oracle.X1.dip
        oracle.theta1.dip <- oracle.Z1.dip %*% t(oracle.Z1.dip)
        oracle.P1.dip <- pmin(pmax(oracle.theta1.dip, 1e-5), 1 - 1e-5)
      }
      oracle.A1.dip <- gen.A.from.P(oracle.P1.dip)
      oracle.g1.dip <- graph.adjacency(oracle.A1.dip, "undirected")
      oracle.result[[jj]] <- list(A1.oracle = oracle.A1.dip,
                                  P1.oracle = oracle.P1.dip,
                                  g1.oracle = oracle.g1.dip,
                                  X1.oracle = oracle.X1.dip)
    }
  }

  Laplace.result <- list()
  for (jj in 1:L) {
    cat(paste("Calling Laplace with \u03B5=", eps[jj], ".\n", sep = ""))
    X1.Lap <- Add.Laplace(X = X1.hat, eps = eps[jj])
    cat("Finish Laplace.\n")

    if (model == "LSM") {
      alpha1.Lap <- X1.Lap[, 1, drop = FALSE]
      Z1.Lap <- X1.Lap[, -1, drop = FALSE]
      theta1.Lap <- alpha1.Lap %*% t(rep(1, n1)) + rep(1, n1) %*% t(alpha1.Lap) + Z1.Lap %*% t(Z1.Lap)
      P1.Lap <- 1 / (1 + exp(-theta1.Lap))
    } else {
      Z1.Lap <- X1.Lap
      theta1.Lap <- Z1.Lap %*% t(Z1.Lap)
      P1.Lap <- pmin(pmax(theta1.Lap, 1e-5), 1 - 1e-5)
    }
    A1.Lap <- gen.A.from.P(P1.Lap)
    g1.Lap <- graph.adjacency(A1.Lap, "undirected")
    Laplace.result[[jj]] <- list(A1.Lap = A1.Lap,
                                 P1.Lap = P1.Lap,
                                 g1.Lap = g1.Lap,
                                 X1.Lap = X1.Lap)
  }

  return(list(non.private.result = non.private.result, GRAND.result = GRAND.result, Laplace.result = Laplace.result, eps = eps, oracle.result = oracle.result))
}

GRAND.evaluate.degree <- function(result) {
  degree.true <- log(1 + degree(result$non.private.result$g1))
  degree.hat <- log(1 + degree(result$non.private.result$g1.hat))
  degree.true2 <- log(1 + degree(result$non.private.result$g2))
  degree.hat2 <- log(1 + degree(result$non.private.result$g2.hat))
  degree.grand <- lapply(result$GRAND.result, function(x) log(1 + degree(x$g1.grand)))
  degree.lap <- lapply(result$Laplace.result, function(x) log(1 + degree(x$g1.Lap)))
  degree.mat <- data.frame(metric = rep("Node Degree", length(result$eps)),
                           eps = result$eps,
                           Hat = rep(wasserstein1d(degree.true, degree.hat), length(result$eps)),
                           Hat2 = rep(wasserstein1d(degree.true2, degree.hat2), length(result$eps)),
                           GRAND = unlist(lapply(degree.grand, function(x) wasserstein1d(degree.true, x))),
                           Laplace = unlist(lapply(degree.lap, function(x) wasserstein1d(degree.true, x))))

  if (length(result$oracle.result) > 0) {
    degree.oracle <- lapply(result$oracle.result, function(x) log(1 + degree(x$g1.oracle)))
    degree.mat$Oracle <- unlist(lapply(degree.oracle, function(x) wasserstein1d(degree.true, x)))
  } else {
    degree.mat$Oracle <- NA
  }

  return(degree.mat)
}

GRAND.evaluate.triangle <- function(result) {
  tri.true <- log(1 + count_triangles(result$non.private.result$g1))
  tri.hat <- log(1 + count_triangles(result$non.private.result$g1.hat))
  tri.true2 <- log(1 + count_triangles(result$non.private.result$g2))
  tri.hat2 <- log(1 + count_triangles(result$non.private.result$g2.hat))
  tri.grand <- lapply(result$GRAND.result, function(x) log(1 + count_triangles(x$g1.grand)))
  tri.lap <- lapply(result$Laplace.result, function(x) log(1 + count_triangles(x$g1.Lap)))
  tri.mat <- data.frame(metric = rep("Triangle Count", length(result$eps)),
                        eps = result$eps,
                        Hat = rep(wasserstein1d(tri.true, tri.hat), length(result$eps)),
                        Hat2 = rep(wasserstein1d(tri.true2, tri.hat2), length(result$eps)),
                        GRAND = unlist(lapply(tri.grand, function(x) wasserstein1d(tri.true, x))),
                        Laplace = unlist(lapply(tri.lap, function(x) wasserstein1d(tri.true, x))))

  if (length(result$oracle.result) > 0) {
    tri.oracle <- lapply(result$oracle.result, function(x) log(1 + count_triangles(x$g1.oracle)))
    tri.mat$Oracle <- unlist(lapply(tri.oracle, function(x) wasserstein1d(tri.true, x)))
  } else {
    tri.mat$Oracle <- NA
  }

  return(tri.mat)
}

GRAND.evaluate.vshape <- function(result) {
  vs.true <- log(1 + get.v(result$non.private.result$g1))
  vs.hat <- log(1 + get.v(result$non.private.result$g1.hat))
  vs.true2 <- log(1 + get.v(result$non.private.result$g2))
  vs.hat2 <- log(1 + get.v(result$non.private.result$g2.hat))
  vs.grand <- lapply(result$GRAND.result, function(x) log(1 + get.v(x$g1.grand)))
  vs.lap <- lapply(result$Laplace.result, function(x) log(1 + get.v(x$g1.Lap)))
  vs.mat <- data.frame(metric = rep("V-Shape Count", length(result$eps)),
                       eps = result$eps,
                       Hat = rep(wasserstein1d(vs.true, vs.hat), length(result$eps)),
                       Hat2 = rep(wasserstein1d(vs.true2, vs.hat2), length(result$eps)),
                       GRAND = unlist(lapply(vs.grand, function(x) wasserstein1d(vs.true, x))),
                       Laplace = unlist(lapply(vs.lap, function(x) wasserstein1d(vs.true, x))))

  if (length(result$oracle.result) > 0) {
    vs.oracle <- lapply(result$oracle.result, function(x) log(1 + get.v(x$g1.oracle)))
    vs.mat$Oracle <- unlist(lapply(vs.oracle, function(x) wasserstein1d(vs.true, x)))
  } else {
    vs.mat$Oracle <- NA
  }

  return(vs.mat)
}

GRAND.evaluate.eigen <- function(result) {
  eigen.true <- eigen_centrality(result$non.private.result$g1)$vector
  eigen.hat <- eigen_centrality(result$non.private.result$g1.hat)$vector
  eigen.true2 <- eigen_centrality(result$non.private.result$g2)$vector
  eigen.hat2 <- eigen_centrality(result$non.private.result$g2.hat)$vector
  eigen.grand <- lapply(result$GRAND.result, function(x) eigen_centrality(x$g1.grand)$vector)
  eigen.lap <- lapply(result$Laplace.result, function(x) eigen_centrality(x$g1.Lap)$vector)
  eigen.mat <- data.frame(metric = rep("Eigen Centrality", length(result$eps)),
                          eps = result$eps,
                          Hat = rep(wasserstein1d(eigen.true, eigen.hat), length(result$eps)),
                          Hat2 = rep(wasserstein1d(eigen.true2, eigen.hat2), length(result$eps)),
                          GRAND = unlist(lapply(eigen.grand, function(x) wasserstein1d(eigen.true, x))),
                          Laplace = unlist(lapply(eigen.lap, function(x) wasserstein1d(eigen.true, x))))

  if (length(result$oracle.result) > 0) {
    eigen.oracle <- lapply(result$oracle.result, function(x) eigen_centrality(x$g1.oracle)$vector)
    eigen.mat$Oracle <- unlist(lapply(eigen.oracle, function(x) wasserstein1d(eigen.true, x)))
  } else {
    eigen.mat$Oracle <- NA
  }

  return(eigen.mat)
}

GRAND.evaluate.harmonic <- function(result) {
  harmonic.true <- harmonic_centrality(result$non.private.result$g1)
  harmonic.hat <- harmonic_centrality(result$non.private.result$g1.hat)
  harmonic.true2 <- harmonic_centrality(result$non.private.result$g2)
  harmonic.hat2 <- harmonic_centrality(result$non.private.result$g2.hat)
  harmonic.grand <- lapply(result$GRAND.result, function(x) harmonic_centrality(x$g1.grand))
  harmonic.lap <- lapply(result$Laplace.result, function(x) harmonic_centrality(x$g1.Lap))
  harmonic.mat <- data.frame(metric = rep("Harmonic Centrality", length(result$eps)),
                             eps = result$eps,
                             Hat = rep(wasserstein1d(harmonic.true, harmonic.hat), length(result$eps)),
                             Hat2 = rep(wasserstein1d(harmonic.true2, harmonic.hat2), length(result$eps)),
                             GRAND = unlist(lapply(harmonic.grand, function(x) wasserstein1d(harmonic.true, x))),
                             Laplace = unlist(lapply(harmonic.lap, function(x) wasserstein1d(harmonic.true, x))))

  if (length(result$oracle.result) > 0) {
    harmonic.oracle <- lapply(result$oracle.result, function(x) harmonic_centrality(x$g1.oracle))
    harmonic.mat$Oracle <- unlist(lapply(harmonic.oracle, function(x) wasserstein1d(harmonic.true, x)))
  } else {
    harmonic.mat$Oracle <- NA
  }

  return(harmonic.mat)
}

get.v <- function(g) {
  deg_vec <- degree(g)
  twostar_counts <- choose(deg_vec, 2)

  return(twostar_counts)
}

#' Evaluate GRAND Results
#'
#' @title Evaluate GRAND Privatization Results
#' @description Evaluates the quality of GRAND privatization results by comparing
#' various network statistics between the original and privatized networks using
#' Wasserstein distance.
#' @param result List. Output from GRAND.privatize function containing privatization results.
#' @param statistics Character vector. Network statistics to evaluate. Options include:
#' "degree", "triangle", "vshape", "eigen", "harmonic". Default is all statistics.
#' @return A data frame containing evaluation results with columns:
#' \itemize{
#'   \item metric: Type of network statistic evaluated
#'   \item eps: Privacy parameter used
#'   \item Hat: Wasserstein distance for non-private estimation
#'   \item Hat2: Wasserstein distance for holdout set estimation
#'   \item GRAND: Wasserstein distance for GRAND privatization
#'   \item Laplace: Wasserstein distance for Laplace mechanism
#'   \item Oracle: Wasserstein distance for oracle method (if available)
#' }
#' @export
#' @examples
#' # Generate and privatize a network
#' dt <- LSM.Gen(n = 1000, k = 3, K = 5)
#' result <- GRAND.privatize(A = dt$A, K = 3, idx = 1:500, eps = c(1, 2, 5, 10))
#' # Evaluate results for all statistics
#' eval_results <- GRAND.evaluate(result)
#' # Evaluate only degree and triangle statistics
#' eval_results <- GRAND.evaluate(result, statistics = c("degree", "triangle"))
GRAND.evaluate <- function(result, statistics = c("degree", "triangle", "vshape", "eigen", "harmonic")) {
  statistic_funcs <- list(degree = GRAND.evaluate.degree,
                          triangle = GRAND.evaluate.triangle,
                          vshape = GRAND.evaluate.vshape,
                          eigen = GRAND.evaluate.eigen,
                          harmonic = GRAND.evaluate.harmonic)

  selected_funcs <- statistic_funcs[statistics]
  results <- lapply(selected_funcs, function(f) f(result))
  output <- do.call(rbind, results)
  rownames(output) <- NULL

  return(output)
}
