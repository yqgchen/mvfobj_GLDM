# =============================================================================
# File: simulate_data.R
# Purpose: Generate synthetic SPD matrix-valued functional data under the
#          GLDM model with amplitude, phase (subject + component + nuisance),
#          and multiplicative perturbation noise.
# Dependencies: base R only (no external packages required)
# =============================================================================

#' Simulate SPD matrix-valued functional data under the GLDM model
#'
#' @description
#' Generates n subjects x p=4 components of 3x3 SPD-valued trajectories on
#' obsGrid, following the GLDM decomposition:
#'   X_{ij}(t) = omega0 + A_{ij} * (tau(R_{ij}(Psi_j(H_i(t)))) - omega0)
#' Subject-level warps H_i, component warps Psi_j, nuisance warps R_{ij}, and
#' amplitude scalars A_{ij} are all randomly generated according to the sigma
#' parameters. Observed data Y_{ij}(t) adds multiplicative log-normal noise.
#'
#' @param n            Integer; number of subjects.
#' @param sigma_warp   Numeric; SD of z_i ~ N(0, sigma_warp) controlling
#'                     subject-level exponential warp H_i.
#' @param sigma_dist   Numeric; SD of r ~ N(0, sigma_dist) controlling the
#'                     nuisance (within-component) warp R_{ij}.
#' @param sigma_pert   Numeric; scale for log-normal multiplicative SPD noise.
#' @param axis_variation Numeric; controls the spread of eigenvalues along the
#'                       SPD template trajectory (larger = more variation).
#' @param obsGrid      Numeric vector; common time grid in [0,1].
#' @param seed         Integer or NULL; if non-NULL, passed to set.seed().
#'
#' @return A named list containing grids, warping maps, amplitude matrix,
#'   noiseless trajectories (X), noisy observations (Y), and oracle-aligned
#'   versions (X_Ginv, X_aligned). See "Return" section below for details.
simulate_data <- function(
    n = 15,               # number of subjects
    sigma_warp = 0.5,     # paramater for generating subject-level deformation functions
    sigma_dist = 0.3,     # paramater for generating nuisance time distortion
    sigma_pert = 0.5,     # paramater for generating SPD perturbation
    axis_variation = 5,   # axis variation parameter for the generation of SPD eigenvalues
    obsGrid = seq(0, 1, 0.05), # time grid on [0,1] for data generation
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  # -----------------------------------------------------------------------------
  # Utility functions (local, used internally)
  # -----------------------------------------------------------------------------
  normalize_vec <- function(v) v / sqrt(sum(v^2))  # L2-normalize a vector
  #FDist2 <- function(A, B) sum((A - B)^2)
  
  ## Template tau(t): deterministic SPD trajectory on [0,1]
  #
  #' @description
  #' Constructs a 3x3 SPD matrix at time tt by building an orthogonal frame
  #' (Q = [v1, v2, v3]) with time-varying columns and diagonal eigenvalue matrix D.
  #' Eigenvalues: lambda1 grows linearly, lambda2 and lambda3 decrease, so the
  #' matrix spreads in one direction over time.
  #'
  #' @param tt Scalar in [0,1]: the time point.
  #' @param axis_variation Numeric: controls eigenvalue spread (slope of lambda1).
  #' @return A 3x3 symmetric positive definite matrix.
  generate_spd_matrix_3d <- function(tt, axis_variation) {
    # Define two time-varying unit vectors in R^3 (not necessarily orthogonal before normalization)
    v2_fun <- function(s) c(cos(pi/4 + pi*s) / sqrt(2),
                            cos(pi/4 + pi*s) / sqrt(2),
                            sin(pi/4 + pi*s))
    v3_fun <- function(s) c(cos(3*pi/4 + pi*s) / sqrt(2),
                            cos(3*pi/4 + pi*s) / sqrt(2),
                            sin(3*pi/4 + pi*s))
    v2 <- normalize_vec(v2_fun(tt))
    v3 <- normalize_vec(v3_fun(tt))
    v1 <- normalize_vec(c(1, -1, 0))  # fixed eigenvector direction
    Q <- cbind(v1, v2, v3)            # orthogonal frame (columns = eigenvectors)
    # Q <- qr.Q(qr(cbind(v1, v2, v3)))  # ensure orthonormal frame
    # if (det(Q) < 0) Q[, 1] <- -Q[, 1]
    # Eigenvalues: first grows, remaining two shrink reciprocally
    D <- diag(c(
      1+axis_variation*tt,
      1/(1+axis_variation*tt),
      1/(1+axis_variation*tt)
    ))
    Q %*% D %*% t(Q)  # spectral decomposition: Q D Q^T
  }

  # Linear interpolation / composition of warping maps on a grid
  compose_map <- function(f_xy, x_grid, xout) approx(x_grid, f_xy, xout = xout, rule = 2)$y

  ## Calculate a SPD matrix to a power
  #
  # Symmetrizes A, then applies spectral decomposition to compute A^pow.
  # pmax(..., 0) clamps eigenvalues to be non-negative for numerical safety.
  mat_sqrt <- function(A, pow = 0.5) {
    eig <- eigen((A + t(A))/2, symmetric = TRUE)
    V <- eig$vectors; d <- pmax(eig$values, 0)
    V %*% diag(d^pow) %*% t(V)
  }

  ## Random orthogonal matrix (Haar-uniform)
  #
  # Generates a random d×d orthogonal matrix via QR decomposition of a
  # Gaussian random matrix. The sign flip on the first column ensures det(Q) = +1.
  rand_orth <- function(d = 3) {
    M <- matrix(rnorm(d*d), d, d)
    Q <- qr.Q(qr(M))
    if (det(Q) < 0) Q[, 1] <- -Q[, 1]
    Q
  }

  ## Multiplicative SPD-preserving noise (mean-preserving)
  #
  # Adds multiplicative log-normal noise to an SPD matrix X:
  #   Y = X^{1/2} S X^{1/2},  where S = U diag(xi) U^T
  # xi ~ LogNormal with mean-log = -sigma^2/2 so E[xi] = 1 (mean-preserving).
  # sigma = sigma_pert/10 ensures small noise when sigma_pert is moderate.
  makeY <- function(X, U, sigma_pert = 0.2) {
    d <- 3
    sigma <- sigma_pert/10          # convert perturbation parameter to log-scale SD
    mu <- - sigma^2 / 2             # mean-log ensuring E[xi_k] = 1
    xi <- rlnorm(d, meanlog = mu, sdlog = sigma)   # d log-normal eigenvalue multipliers
    S  <- U %*% diag(xi) %*% t(U)  # random SPD noise matrix in the eigenbasis of U
    Xh <- mat_sqrt(X, 0.5)         # X^{1/2}
    Y  <- Xh %*% S %*% Xh          # congruence transform: X^{1/2} S X^{1/2}
    (Y + t(Y)) / 2                  # symmetrize for numerical stability
  }
  
  # # Calculate true c
  # true_c <- function(X, omega0) {
  #   # sup_t distance for one 3D array (d x d x T)
  #   sup_dist <- function(arr, omega0) {
  #     Tlen <- dim(arr)[3]
  #     max(vapply(seq_len(Tlen), function(tk) sqrt(FDist2(arr[,,tk], omega0)), numeric(1)))
  #   }
  #   
  #   c_star <- Inf
  #   for (j in seq_along(X)) {
  #     for (i in seq_along(X[[j]])) {
  #       val <- sup_dist(X[[j]][[i]], omega0)
  #       if (val < c_star) c_star <- val
  #     }
  #   }
  #   c_star
  # }
  
  # -----------------------------------------------------------------------------
  # Common grids and template
  # -----------------------------------------------------------------------------
  # Evaluate the deterministic SPD template tau(t) on all time points in obsGrid
  tau <- lapply(obsGrid, function(tt) generate_spd_matrix_3d(tt, axis_variation))
  tau <- simplify2array(tau)  # converts list of 3x3 matrices to a 3x3xT array
  omega0 <- diag(0, 3)        # reference matrix (zero matrix = "origin" of the SPD space)

  # -----------------------------------------------------------------------------
  # Warping maps H, Psi, R
  # -----------------------------------------------------------------------------

  # Subject-level H_i, H_i^{-1}
  # H_i is an exponential warp: H_i^{-1}(t) = (e^{z_i t} - 1)/(e^{z_i} - 1)
  # When z_i ~ 0, the warp degenerates to identity (obsGrid).
  m <- length(obsGrid)
  hinvmat <- hmat <- matrix(0, n, m)
  z_subj <- rnorm(n, 0, sigma_warp)  # one random scalar per subject
  for (i in 1:n) {
    # Inverse warp H_i^{-1}: identity when z_i ~ 0, exponential otherwise
    hinvmat[i, ] <- if (abs(z_subj[i]) < .Machine$double.eps)
      obsGrid else (exp(obsGrid * z_subj[i]) - 1)/(exp(z_subj[i]) - 1)
    # Forward warp H_i: numerical inverse of H_i^{-1} via interpolation
    hmat[i, ] <- approx(hinvmat[i, ], obsGrid, obsGrid, rule = 2)$y
  }

  # Component-level Psi_j, Psi_j^{-1}
  # p=4 components; Psi_j is a beta-CDF mixture warp for j=1,2;
  # Psi_{j+2} is its "reflected" counterpart: Psi_{j+2}^{-1}(t) = 2t - Psi_j^{-1}(t)
  p <- 4
  alpha <- c(2, 1)    # beta distribution shape parameters (one per pair)
  beta <- c(2, 1/2)   # beta distribution rate parameters
  lambda <- 0.5       # mixing weight between beta-CDF warp and identity
  psi <- psi_inv <- vector("list", p)
  for (j in 1:(p/2)) {
    # Psi_j(t) = lambda * F_beta(t; alpha_j, beta_j) + (1-lambda) * t
    psi[[j]] <- lambda * pbeta(obsGrid, alpha[j], beta[j]) + (1 - lambda)*obsGrid
    psi_inv[[j]] <- approx(psi[[j]], obsGrid, obsGrid, rule = 2)$y
    # Reflected inverse: Psi_{j+2}^{-1}(t) = 2t - Psi_j^{-1}(t)
    psi_inv[[j+2]] <- 2*obsGrid - psi_inv[[j]]
    psi[[j+2]] <- approx(psi_inv[[j+2]], obsGrid, obsGrid, rule = 2)$y
  }

  # R_{i,j} and R_{i,j}^{-1} (nuisance within-component warp)
  # Same exponential family as H_i, but independently drawn per (i,j) pair.
  rmat <- rinv <- vector("list", n)
  for (i in 1:n) {
    rmat[[i]] <- vector("list", p)
    rinv[[i]] <- vector("list", p)
    for (j in 1:p) {
      r <- rnorm(1, 0, sigma_dist)  # random nuisance warp parameter for (i,j)
      # inverse map R_{ij}^{-1}
      rinv_ij <- if (abs(r) < .Machine$double.eps)
        obsGrid else (exp(obsGrid * r) - 1)/(exp(r) - 1)
      rinv[[i]][[j]] <- rinv_ij
      # forward map R_{ij} obtained as numerical inverse of R_{ij}^{-1}
      rmat[[i]][[j]] <- approx(rinv_ij, obsGrid, obsGrid, rule = 2)$y
    }
  }

  # -----------------------------------------------------------------------------
  # Component tempo trajectories (true aligned component means)
  # -----------------------------------------------------------------------------
  # eta_j(t) = tau(Psi_j(t)): template warped by component warp Psi_j
  Leta <- vector("list", p)
  for (j in seq_len(p)) {
    stopifnot(is.numeric(psi[[j]]))
    eta_j <- lapply(psi[[j]], function(tt) generate_spd_matrix_3d(tt, axis_variation))
    Leta[[j]] <- simplify2array(eta_j)  # d x d x T array for component j
  }


  # -----------------------------------------------------------------------------
  # Amplitudes A_{i,j}
  # -----------------------------------------------------------------------------
  amat <- matrix(runif(n*p, 1, 2), n, p)  # random amplitude scalars, uniform on [1,2]
  
  # -----------------------------------------------------------------------------
  # Generate noiseless SPD trajectories X_{i,j}(t) for t in obsGrid
  # -----------------------------------------------------------------------------

  #' Build a single noiseless SPD trajectory for subject i, component j.
  #'
  #' @description
  #' Evaluates tau at the composed warp s = R_{ij}(Psi_j(H_i(t))), then
  #' rescales: X_{ij}(t) = omega0 + A_ij * (tau(s(t)) - omega0).
  #'
  #' @param A_ij    Scalar amplitude for (i,j).
  #' @param omega0  d x d reference matrix.
  #' @param rmat_ij Numeric vector of length T: R_{ij}(t) on obsGrid.
  #' @param h_row   Numeric vector of length T: H_i(t) on obsGrid.
  #' @param psi_j   Numeric vector of length T: Psi_j(t) on obsGrid.
  #' @param obsGrid Numeric vector (time grid).
  #' @return d x d x T array.
  make_X <- function(A_ij, omega0, rmat_ij, h_row, psi_j, obsGrid) {
    d <- nrow(omega0)
    arr <- array(0, c(d,d,length(obsGrid)))
    # s = (R_{i,j} ∘ Psi_j ∘ H_i)(t) for t in obsGrid (composed warps)
    s <- compose_map(
      rmat_ij, obsGrid,
      compose_map(psi_j, obsGrid, h_row)
    )
    for (k in seq_along(obsGrid)) {
      # Amplitude-scaled and origin-shifted SPD matrix at time t_k
      arr[,,k] <- A_ij * ( generate_spd_matrix_3d(s[k], axis_variation) - omega0 ) + omega0
    }
    arr
  }

  # Build noiseless X: component-first list [[j]][[i]], each d x d x T
  simdata_X <- vector("list", p)
  for (j in 1:p) {
    simdata_X[[j]] <- vector("list", n)
    for (i in 1:n)
      simdata_X[[j]][[i]] <- make_X(amat[i,j], omega0, rmat[[i]][[j]], hmat[i,], psi[[j]], obsGrid)
  }

  # -----------------------------------------------------------------------------
  # Two versions of "true aligned" trajectories (for evaluation)
  # -----------------------------------------------------------------------------

  d    <- nrow(omega0)   # matrix dimension (3)

  # Fully aligned: X_{ij}(G_{ij}^{-1}(t)) = omega0 + A_ij*(tau(t) - omega0)
  # (all warps removed; oracle alignment target)
  make_X_aligned_array <- function(A_ij, omega0, obsGrid) {
    d <- nrow(omega0)
    arr <- array(0, c(d,d,length(obsGrid)))
    for (k in seq_along(obsGrid)) {
      arr[,,k] <- A_ij * ( generate_spd_matrix_3d(obsGrid[k], axis_variation) - omega0 ) + omega0
    }
    arr
  }

  # Partially aligned: X_{ij} evaluated at R_{ij}^{-1}(t) only (nuisance removed but H, Psi remain)
  # X(G_{i,j}^{-1}(t)) with G^{-1}_{ij} = Psi_j^{-1}(H_i^{-1}(...)) — here simplified to rmat_ij only
  make_X_Ginv_array <- function(A_ij, omega0, rmat_ij, obsGrid) {
    d <- nrow(omega0)
    arr <- array(0, c(d,d,length(obsGrid)))
    for (k in seq_along(obsGrid)) {
      # rmat_ij[k] = R_{ij}(t_k) plays the role of the warped time argument
      arr[,,k] <- A_ij * ( generate_spd_matrix_3d(rmat_ij[k], axis_variation) - omega0 ) + omega0
    }
    arr
  }

  # Build both oracle-aligned arrays: component-first [[j]][[i]]
  X_Ginv <- vector("list", p)
  X_aligned <- vector("list", p)
  for (j in seq_len(p)) {
    X_Ginv[[j]] <- vector("list", n)
    X_aligned[[j]] <- vector("list", n)
    for (i in seq_len(n)) {
      X_Ginv[[j]][[i]] <- make_X_Ginv_array(amat[i,j], omega0, rmat[[i]][[j]], obsGrid)   # d×d×T
      X_aligned[[j]][[i]] <- make_X_aligned_array(amat[i,j], omega0, obsGrid)              # d×d×T
    }
  }

  # -----------------------------------------------------------------------------
  # Add random perturbation noise to get observed data Y
  # -----------------------------------------------------------------------------
  U <- rand_orth(3)   # shared random orthogonal matrix for all noise realizations
  simdata_Y <- vector("list", p)
  for (j in 1:p) {
    simdata_Y[[j]] <- vector("list", n)
    for (i in 1:n) {
      Tn <- dim(simdata_X[[j]][[i]])[3]
      arrY <- array(0, c(d,d,Tn))
      for (k in 1:Tn)
        arrY[,,k] <- makeY(simdata_X[[j]][[i]][,,k], U, sigma_pert)  # add noise at each time point
      simdata_Y[[j]][[i]] <- arrY
    }
  }
  
  ## ------------------------------
  ## Return: all necessary components
  ## (All SPD trajectories are 3×3×T arrays nested as [[i]][[j]])
  ## ------------------------------
  list(
    # --- Grids and template ---
    obsGrid     = obsGrid,       # Observation time grid (length m)
    tau  = tau,    # True template trajectory {τ(t)} on obsGrid
    Leta       = Leta,         # List of true component tempo trajectories        
    
    # --- Warping components ---
    psi        = psi,          # Component-level forward maps {Ψ_j} on t_Grid
    psi_inv    = psi_inv,      # Component-level inverse maps {Ψ_j^{-1}}
    hmat       = hmat,         # n×m subject-level warping H_i(t) on obsGrid
    hinvmat    = hinvmat,      # n×m subject-level inverse warping H_i^{-1}(t)
    rmat       = rmat,         # Forward nuisance warping {R_{ij}(t)}
    rinv       = rinv,         # Inverse nuisance warping {R_{ij}^{-1}(t)}
    
    # --- Parameters and constants ---
    amat       = amat,         # n×p matrix of amplitude scalars a_{ij}
    omega0    = omega0,           # Reference (origin) SPD matrix (here: zero matrix)
    U_fixed    = U,            # Fixed orthogonal matrix used for perturbation noise
    
    # --- Generated SPD trajectories ---
    X          = simdata_X,            # Noiseless trajectories X_{ij}(t)
    Y          = simdata_Y,            # Contaminated SPD trajectories Y_{ij}(t)
    
    # --- Aligned versions (for analysis) ---
    X_Ginv  = X_Ginv,  # Partially aligned trajectories (G-only)
    X_aligned  = X_aligned   # Fully aligned (oracle) trajectories
  )
}