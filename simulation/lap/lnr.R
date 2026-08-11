# =============================================================================
# File: lnr.R
# Purpose: Local network regression (LNR) for graph Laplacian responses.
#          Source: Zhou, Y. and Müller, H.-G. (2022). Network regression with
#          graph Laplacians. Journal of Machine Learning Research, 23(1), 1-41.
#          The code here is adapted from the supplementary material of that paper.
# Dependencies: osqp (QP solver for Laplacian projection), pracma (for internal
#               CV helpers)
# =============================================================================

# -----------------------------------------------------------------------------
# Main function: local linear network regression
# -----------------------------------------------------------------------------

#' Local linear regression for graph Laplacian responses
#'
#' @description
#' Fits a local linear regression model where the response is a graph Laplacian
#' matrix and the covariates are scalar or bivariate (p = 1 or 2). The fitted
#' value at each observation (or prediction) point is obtained by:
#' (1) computing local-linear kernel weights centered at the query point,
#' (2) forming a weighted linear combination of the vectorized Laplacians,
#' (3) projecting the result onto the set of valid graph Laplacians via OSQP.
#' Bandwidth can be user-supplied or selected by leave-one-out cross-validation.
#'
#' @param gl    List of n graph Laplacian matrices (each m x m), OR a 3D array
#'              (m x m x n) which is converted to a list internally.
#' @param x     n x p numeric matrix (or data.frame/vector) of covariate values;
#'              p must be 1 or 2.
#' @param xOut  Optional n_out x p matrix of prediction points; if NULL, only
#'              in-sample fits are returned.
#' @param optns Named list of options:
#'   - metric: "frobenius" (default) or "power" (Euclidean power metric).
#'   - alpha:  Numeric >= 0; power exponent used when metric = "power" (default 1).
#'   - kernel: Kernel type: "gauss" (default), "rect", "epan", "gausvar", "quar".
#'   - bw:     Numeric vector of length p; bandwidth(s). If NA (default), selected
#'             by leave-one-out cross-validation over a 20-point log-spaced grid
#'             per dimension.
#'   - digits: Integer; if not NA, rounds the predicted Laplacian entries.
#'
#' @return An object of class "nr" (a list) with:
#'   - fit:       List of length n; fitted Laplacian at each observation point.
#'   - predict:   (only if xOut is not NULL) List of length n_out; predictions.
#'   - residuals: Numeric vector length n; Frobenius norm of fit residuals.
#'   - gl, x, xOut, optns: inputs/options (returned for reference).
lnr <- function(gl = NULL, x = NULL, xOut = NULL, optns = list()) {
  if (is.null(gl) | is.null(x)) {
    stop("requires the input of both gl and x")
  }
  if (is.null(optns$metric)) {
    optns$metric <- "frobenius"
  }
  if (!(optns$metric %in% c("frobenius", "power"))) {
    stop("metric choice not supported")
  }
  if (is.null(optns$alpha)) {
    optns$alpha <- 1
  }
  if (optns$alpha < 0) {
    stop("alpha must be non-negative")
  }
  if (is.null(optns$kernel)) {
    optns$kernel <- "gauss"
  }
  if (is.null(optns$bw)) {
    optns$bw <- NA
  }
  if (is.null(optns$digits)) {
    optns$digits <- NA
  }
  if (!is.matrix(x)) {
    if (is.data.frame(x) | is.vector(x)) {
      x <- as.matrix(x)
    } else {
      stop("x must be a matrix or a data frame or a vector")
    }
  }
  n <- nrow(x) # number of observations
  p <- ncol(x) # number of covariates
  if (p > 2) {
    stop("local method is designed to work in low dimensional case (p is either 1 or 2)")
  }
  if (!is.na(sum(optns$bw))) {
    if (sum(optns$bw <= 0) > 0) {
      stop("bandwidth must be positive")
    }
    if (length(optns$bw) != p) {
      stop("dimension of bandwidth does not agree with x")
    }
  }
  if (!is.list(gl)) {
    if (is.array(gl)) {
      gl <- lapply(seq(dim(gl)[3]), function(i) gl[, , i])
    } else {
      stop("gl must be a list or an array")
    }
  }
  if (length(gl) != n) {
    stop("the number of rows in x must be the same as the number of graph Laplacians in gl")
  }
  if (!is.null(xOut)) {
    if (!is.matrix(xOut)) {
      if (is.data.frame(xOut)) {
        xOut <- as.matrix(xOut)
      } else if (is.vector(xOut)) {
        if (p == 1) {
          xOut <- as.matrix(xOut)
        } else {
          xOut <- t(xOut)
        }
      } else {
        stop("xOut must be a matrix or a data frame or a vector")
      }
    }
    if (ncol(xOut) != p) {
      stop("x and xOut must have the same number of columns")
    }
    nOut <- nrow(xOut) # number of predictions
  } else {
    nOut <- 0
  }
  nodes <- colnames(gl[[1]])
  m <- ncol(gl[[1]]) # number of nodes in the graph
  if (is.null(nodes)) nodes <- 1:m
  glVec <- matrix(unlist(gl), ncol = m^2, byrow = TRUE) # n by m^2: vectorized Laplacians
  if (substr(optns$metric, 1, 1) == "p") {
    # For power metric: pre-compute L^alpha for each observation using eigendecomposition
    glAlpha <- lapply(gl, function(gli) {
      eigenDecom <- eigen(gli)
      Lambda <- pmax(Re(eigenDecom$values), 0) # exclude 0i
      U <- eigenDecom$vectors
      U %*% diag(Lambda^optns$alpha) %*% t(U)
    })
    glAlphaVec <- matrix(unlist(glAlpha), ncol = m^2, byrow = TRUE) # n by m^2
  }
  
  # --- OSQP solver initialization for Laplacian projection ---
  # The projection of a symmetric matrix onto the set of valid graph Laplacians
  # is a quadratic program: minimize ||L_vec - q||^2 s.t. L is symmetric,
  # off-diagonal entries <= 0, and row sums = 0.
  W <- 2^32 # large upper bound on edge weights (effectively unconstrained above)
  nConsts <- m^2 # total number of linear constraints
  # Constraint lower bounds: symmetry (=0), row-sum (=0), nonpos off-diag (<= 0 encoded as >= -W)
  l <- c(rep.int(0, m * (m + 1) / 2), rep.int(-W, m * (m - 1) / 2))
  u <- rep.int(0, nConsts)  # all upper bounds = 0
  q <- rep.int(0, m^2)      # linear part of QP (updated at each query point)
  P <- diag(m^2)            # identity quadratic part: minimize ||L_vec - q||^2
  consts <- matrix(0, nrow = nConsts, ncol = m^2)
  k <- 0
  # Symmetry constraints: L[j,i] - L[i,j] = 0 for all i < j
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      k <- k + 1
      consts[k, (j - 1) * m + i] <- 1
      consts[k, (i - 1) * m + j] <- -1
    }
  }
  # Row-sum constraints: sum_j L[i,j] = 0 for each row i (Laplacian property)
  for (i in 1:m) {
    consts[k + i, ((i - 1) * m + 1):(i * m)] <- rep(1, m)
  }
  k <- k + m
  # Non-positivity constraints on off-diagonal entries: L[j,i] <= 0
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      k <- k + 1
      consts[k, (j - 1) * m + i] <- 1
    }
  }
  model <- osqp::osqp(P, q, consts, l, u, osqp::osqpSettings(verbose = FALSE))
  
  # --- Kernel function setup ---
  # Kern is a univariate kernel; K is the product kernel for p-dimensional x.
  Kern <- kerFctn(optns$kernel)
  K <- function(x, h) {
    k <- 1
    for (i in 1:p) {
      k <- k * Kern(x[i] / h[i])  # product kernel: K(x, h) = prod_l Kern(x_l / h_l)
    }
    return(as.numeric(k))
  }

  # --- Bandwidth selection by leave-one-out cross-validation ---
  # SetBwRange computes a safe [min, max] candidate range for each dimension
  SetBwRange <- function(xin, xout, kernel_type) {
    xinSt <- sort(unique(as.numeric(xin)))
    if (length(xinSt) < 2L) {
      # Degenerate case: only one unique sample point
      rng <- max(1e-3, diff(range(c(xin, xout))))
      return(list(min = rng * 0.1, max = rng * 0.5))
    }
    # Use nearest gap and boundary distances to set lower bound
    dx_min  <- min(diff(xinSt))
    left    <- xinSt[2] - min(xout)
    right   <- max(xout) - xinSt[length(xinSt)-1]
    s       <- max(dx_min, left, right) * 1.1
    
    # Kernel-dependent scaling
    denom <- (ifelse(kernel_type == "gauss", 3, 1) *
                ifelse(kernel_type == "gausvar", 2.5, 1))
    bw.min <- s / denom
    
    bw.max <- diff(range(xin)) / 3
    if (bw.max < bw.min) {
      bw.max <- if (bw.min > bw.max * 1.5) bw.min * 1.01 else bw.max * 1.5
    }
    list(min = bw.min, max = bw.max)
  }
  
  # ---------- Bandwidth selection by cross-validation (uses SetBwRange & K) ----------
  if (is.na(sum(optns$bw))) {
    p <- ncol(x); n <- nrow(x)
    
    # (1) Per-dimension bandwidth ranges from your SetBwRange()
    bw_ranges <- lapply(seq_len(p), function(l) SetBwRange(x[, l], x[, l], optns$kernel))
    
    # (2) Log-spaced candidate grid: 20 per dimension (safe lower/upper bounds)
    hs <- matrix(0, p, 20)
    for (l in 1:p) {
      lo <- max(bw_ranges[[l]]$min, 1e-6)
      hi <- max(bw_ranges[[l]]$max, lo * 1.05)
      hs[l, ] <- exp(seq(from = log(lo), to = log(hi), length.out = 20))
    }
    
    # (3) CV loop over all combinations
    ngrid   <- 20^p
    cv      <- rep(0, ngrid)
    ridge   <- 1e-8    # diagonal ridge for mu2 inversion
    w_floor <- 1e-12   # minimum neighborhood weight sum
    
    for (k in 0:(ngrid - 1)) {
      # Extract candidate bandwidth vector h
      h <- numeric(p)
      for (l in 1:p) {
        kl <- floor((k %% (20^l)) / (20^(l - 1))) + 1
        h[l] <- hs[l, kl]
      }
      
      bad_candidate <- FALSE
      
      # Leave-one-out CV
      for (j in 1:n) {
        a  <- x[j, , drop = FALSE]      # 1 x p
        xx <- x[-j, , drop = FALSE]     # (n-1) x p
        
        # Centered covariates Z_i = x_i - a
        Z <- sweep(xx, 2, as.numeric(a), FUN = "-")
        
        # Pure kernel weights w0_i = K(Z_i, h) using your K()
        w0 <- apply(Z, 1, function(z) K(z, h))
        sw <- sum(w0)
        if (!is.finite(sw) || sw < w_floor) {
          cv[k + 1] <- Inf
          bad_candidate <- TRUE
          break
        }
        
        # Local-linear moments (vectorized)
        # mu1 = E[K(Z) * Z], mu2 = E[K(Z) * (Z Z^T)]
        mu1 <- colSums(sweep(Z, 1, w0, `*`)) / length(w0)
        Zw  <- sweep(Z, 1, sqrt(w0), `*`)
        mu2 <- crossprod(Zw) / length(w0)               # p x p
        
        # Stabilize mu2 and invert
        mu2 <- as.matrix(mu2); diag(mu2) <- diag(mu2) + ridge
        inv_mu2 <- tryCatch(solve(mu2), error = function(e) NULL)
        if (is.null(inv_mu2) || any(!is.finite(inv_mu2))) {
          cv[k + 1] <- Inf
          bad_candidate <- TRUE
          break
        }
        
        # Local-linear correction weights: w_i = K(Z_i)*(1 - wc %*% Z_i)
        wc   <- drop(mu1 %*% inv_mu2)                   # length p
        corr <- 1 - as.numeric(Z %*% wc)                # length n-1
        w    <- w0 * corr
        if (!all(is.finite(w)) || sum(abs(w)) < w_floor) {
          cv[k + 1] <- Inf
          bad_candidate <- TRUE
          break
        }
        
        # Accumulate CV loss
        if (substr(optns$metric, 1, 1) == "f") {
          # Frobenius metric branch
          qNew <- apply(glVec[-j, ], 2, weighted.mean, w = w)  # vectorized target
          model$Update(q = -qNew)
          fitj <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
          fitj <- (fitj + t(fitj)) / 2                         # symmetrize
          if (!is.na(optns$digits)) fitj <- round(fitj, digits = optns$digits)
          # Enforce Laplacian structure
          fitj[fitj > 0] <- 0
          diag(fitj) <- 0
          diag(fitj) <- -colSums(fitj)
          cv[k + 1] <- cv[k + 1] + sum((gl[[j]] - fitj)^2) / n
          
        } else if (substr(optns$metric, 1, 1) == "p") {
          # Euclidean power metric branch
          bAlpha <- matrix(apply(glAlphaVec[-j, ], 2, weighted.mean, w = w), ncol = m)
          ed     <- eigen(bAlpha, symmetric = TRUE)
          Lambda <- pmax(Re(ed$values), 0)
          U      <- ed$vectors
          qNew   <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U))
          model$Update(q = -qNew)
          fitj <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
          fitj <- (fitj + t(fitj)) / 2
          if (!is.na(optns$digits)) fitj <- round(fitj, digits = optns$digits)
          fitj[fitj > 0] <- 0
          diag(fitj) <- 0
          diag(fitj) <- -colSums(fitj)
          cv[k + 1] <- cv[k + 1] + sum((gl[[j]] - fitj)^2) / n
        }
      } # end for j
      
      if (bad_candidate) next
    } # end for k
    
    # (4) Pick best bandwidth; if all invalid, fallback to mid-range per dimension
    if (all(!is.finite(cv))) {
      optns$bw <- sapply(seq_len(p), function(l) {
        exp((log(hs[l, 1]) + log(hs[l, 20])) / 2)
      })
      warning("All CV candidates invalid; using mid-range bandwidths per dimension.")
    } else {
      bwi <- which.min(cv)
      optns$bw <- numeric(p)
      for (l in 1:p) {
        kl <- floor((bwi %% (20^l)) / (20^(l - 1))) + 1
        optns$bw[l] <- hs[l, kl]
      }
    }
  }
  
  
  # --- In-sample fitting and out-of-sample prediction ---
  # For each query point a, compute local-linear weights w_i and project the
  # weighted mean of vectorized Laplacians onto the valid Laplacian cone via OSQP.

  fit <- vector(mode = "list", length = n)
  residuals <- rep.int(0, n)
  if (substr(optns$metric, 1, 1) == "f") {
    # --- Frobenius metric branch ---
    for (i in 1:n) {
      a <- x[i, ]
      # Compute local-linear moments mu1 (p-vector) and mu2 (p x p matrix)
      if (p > 1) {
        mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
      } else {
        mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
      }
      wc <- t(mu1) %*% solve(mu2) # local-linear bias correction vector (1 x p)
      w <- apply(x, 1, function(xi) {
        K(xi - a, optns$bw) * (1 - wc %*% (xi - a))  # local-linear effective weight
      })
      qNew <- apply(glVec, 2, weighted.mean, w) # weighted mean of vectorized Laplacians (m^2)
      model$Update(q = -qNew)  # QP: minimize ||L - qNew||^2 over Laplacian cone
      temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
      temp <- (temp + t(temp)) / 2 # symmetrize to correct numerical asymmetry
      if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # optional rounding
      temp[temp > 0] <- 0 # enforce: off-diagonal Laplacian entries must be <= 0
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)  # enforce row-sum = 0 (Laplacian diagonal)
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((gl[[i]] - temp)^2))  # Frobenius residual
    }
    if (nOut > 0) {
      # --- Out-of-sample predictions at xOut query points ---
      predict <- vector(mode = "list", length = nOut)
      for (i in 1:nOut) {
        a <- xOut[i, ]
        if (p > 1) {
          mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(x, 1, function(xi) {
          K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
        }) # local-linear effective weight
        qNew <- apply(glVec, 2, weighted.mean, w) # m^2
        model$Update(q = -qNew)
        temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
        temp[temp > 0] <- 0 # off diagonal should be negative
        diag(temp) <- 0
        diag(temp) <- -colSums(temp)
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, residuals = residuals, gl = gl, x = x, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, residuals = residuals, gl = gl, x = x, optns = optns)
    }
  } else if (substr(optns$metric, 1, 1) == "p") {
    # --- Euclidean power metric branch ---
    # Uses L^alpha as the Euclidean embedding; inverts the power mapping after
    # taking the weighted mean, then projects onto the Laplacian cone via OSQP.
    for (i in 1:n) {
      a <- x[i, ]
      if (p > 1) {
        mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
      } else {
        mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
        mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
      }
      wc <- t(mu1) %*% solve(mu2) # 1 by p
      w <- apply(x, 1, function(xi) {
        K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
      }) # weight
      bAlpha <- matrix(apply(glAlphaVec, 2, weighted.mean, w), ncol = m) # weighted mean of L^alpha (m x m)
      eigenDecom <- eigen(bAlpha)
      Lambda <- pmax(Re(eigenDecom$values), 0) # project eigenvalues to M_m (nonneg PSD cone)
      U <- eigenDecom$vectors
      qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power: (bAlpha)^{1/alpha}
      model$Update(q = -qNew)
      temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
      temp <- (temp + t(temp)) / 2 # symmetrize
      if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
      temp[temp > 0] <- 0 # off diagonal should be negative
      diag(temp) <- 0
      diag(temp) <- -colSums(temp)
      fit[[i]] <- temp
      residuals[i] <- sqrt(sum((gl[[i]] - temp)^2))
      # Alternative: compute residuals in alpha-power metric space (commented out)
      # eigenDecom <- eigen(fit[[i]])
      # Lambda <- pmax(Re(eigenDecom$values), 0)
      # U <- eigenDecom$vectors
      # fitiAlpha <- U%*%diag(Lambda^optns$alpha)%*%t(U)
      # residuals[i] <- sqrt(sum((glAlpha[[i]]-fitiAlpha)^2))# using Euclidean power metric
    }
    if (nOut > 0) {
      predict <- vector(mode = "list", length = nOut)
      for (i in 1:nOut) {
        a <- xOut[i, ]
        if (p > 1) {
          mu1 <- rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- matrix(rowMeans(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a)))), ncol = p)
        } else {
          mu1 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * (xi - a)))
          mu2 <- mean(apply(x, 1, function(xi) K(xi - a, optns$bw) * ((xi - a) %*% t(xi - a))))
        }
        wc <- t(mu1) %*% solve(mu2) # 1 by p
        w <- apply(x, 1, function(xi) {
          K(xi - a, optns$bw) * (1 - wc %*% (xi - a))
        }) # weight
        bAlpha <- matrix(apply(glAlphaVec, 2, weighted.mean, w), ncol = m) # m x m
        eigenDecom <- eigen(bAlpha)
        Lambda <- pmax(Re(eigenDecom$values), 0) # projection to M_m
        U <- eigenDecom$vectors
        qNew <- as.vector(U %*% diag(Lambda^(1 / optns$alpha)) %*% t(U)) # inverse power
        model$Update(q = -qNew)
        temp <- matrix(model$Solve()$x, ncol = m, dimnames = list(nodes, nodes))
        temp <- (temp + t(temp)) / 2 # symmetrize
        if (!is.na(optns$digits)) temp <- round(temp, digits = optns$digits) # round
        temp[temp > 0] <- 0 # off diagonal should be negative
        diag(temp) <- 0
        diag(temp) <- -colSums(temp)
        predict[[i]] <- temp
      }
      res <- list(fit = fit, predict = predict, residuals = residuals, gl = gl, x = x, xOut = xOut, optns = optns)
    } else {
      res <- list(fit = fit, residuals = residuals, gl = gl, x = x, optns = optns)
    }
  }
  class(res) <- "nr"
  res
}

# -----------------------------------------------------------------------------
# Kernel factory
# -----------------------------------------------------------------------------

#' Return a univariate kernel function by name
#'
#' @description
#' Returns a closure implementing the specified kernel function for use in
#' local regression weight computation. The returned function takes a scalar
#' (or vector) x and returns the kernel evaluation K(x).
#'
#' @param kernel_type Character; one of "gauss", "rect", "epan", "gausvar", "quar".
#'
#' @return A function of one argument x returning the kernel evaluated at x.
kerFctn <- function(kernel_type) {
  if (kernel_type == "gauss") {
    ker <- function(x) {
      dnorm(x) # standard Gaussian kernel: exp(-x^2/2) / sqrt(2*pi)
    }
  } else if (kernel_type == "rect") {
    ker <- function(x) {
      as.numeric((x <= 1) & (x >= -1))
    }
  } else if (kernel_type == "epan") {
    ker <- function(x) {
      n <- 1
      (2 * n + 1) / (4 * n) * (1 - x^(2 * n)) * (abs(x) <= 1)
    }
  } else if (kernel_type == "gausvar") {
    ker <- function(x) {
      dnorm(x) * (1.25 - 0.25 * x^2)
    }
  } else if (kernel_type == "quar") {
    ker <- function(x) {
      (15 / 16) * (1 - x^2)^2 * (abs(x) <= 1)
    }
  } else {
    stop("Unavailable kernel")
  }
  return(ker)
}

