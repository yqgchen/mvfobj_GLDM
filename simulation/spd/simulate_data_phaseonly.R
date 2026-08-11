# =============================================================================
# File: simulate_data_phaseonly.R
# Purpose: Phase-only variant of simulate_data.R — generates synthetic SPD
#          matrix-valued functional data without amplitude variation (A_{ij} = 1).
#          All other warping and noise components are identical to simulate_data.R.
# Dependencies: base R only (no external packages required)
# =============================================================================

#' Simulate SPD matrix-valued functional data — phase-only case (no amplitude variation)
#'
#' @description
#' Identical to simulate_data() except that all amplitude scalars A_{ij} = 1,
#' so only phase (warping) variation and perturbation noise are present.
#' Used for the "phaseonly" simulation setting in the GLDM paper.
#'
#' @param n            Integer; number of subjects.
#' @param sigma_warp   Numeric; SD of z_i ~ N(0, sigma_warp) controlling H_i.
#' @param sigma_dist   Numeric; SD of r ~ N(0, sigma_dist) controlling R_{ij}.
#' @param sigma_pert   Numeric; scale for log-normal multiplicative SPD noise.
#' @param axis_variation Numeric; eigenvalue spread in the SPD template.
#' @param obsGrid      Numeric vector; common time grid in [0,1].
#' @param seed         Integer or NULL; passed to set.seed() if non-NULL.
#'
#' @return Same named list structure as simulate_data(). The amat field will
#'   be an n x p matrix of ones (no amplitude variation).
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
  # Same construction as simulate_data.R: orthogonal frame Q with time-varying
  # columns v2(t), v3(t) and fixed v1; eigenvalues grow/shrink with axis_variation.
  generate_spd_matrix_3d <- function(tt, axis_variation) {
    v2_fun <- function(s) c(cos(pi/4 + pi*s) / sqrt(2),
                            cos(pi/4 + pi*s) / sqrt(2),
                            sin(pi/4 + pi*s))
    v3_fun <- function(s) c(cos(3*pi/4 + pi*s) / sqrt(2),
                            cos(3*pi/4 + pi*s) / sqrt(2),
                            sin(3*pi/4 + pi*s))
    v2 <- normalize_vec(v2_fun(tt))
    v3 <- normalize_vec(v3_fun(tt))
    v1 <- normalize_vec(c(1, -1, 0))  # fixed eigenvector direction
    Q <- cbind(v1, v2, v3)            # orthogonal frame
    # Q <- qr.Q(qr(cbind(v1, v2, v3)))  # ensure orthonormal frame
    # if (det(Q) < 0) Q[, 1] <- -Q[, 1]
    D <- diag(c(
      1+axis_variation*tt,
      1/(1+axis_variation*tt),
      1/(1+axis_variation*tt)
    ))
    Q %*% D %*% t(Q)
  }

  # Linear interpolation / composition of warping maps on a grid
  compose_map <- function(f_xy, x_grid, xout) approx(x_grid, f_xy, xout = xout, rule = 2)$y

  ## Calculate a SPD matrix to a power
  # Symmetrizes A then applies spectral decomposition; pmax clamps negative eigenvalues.
  mat_sqrt <- function(A, pow = 0.5) {
    eig <- eigen((A + t(A))/2, symmetric = TRUE)
    V <- eig$vectors; d <- pmax(eig$values, 0)
    V %*% diag(d^pow) %*% t(V)
  }

  ## Random orthogonal matrix (Haar-uniform)
  # QR decomposition of a Gaussian matrix; sign flip ensures det = +1.
  rand_orth <- function(d = 3) {
    M <- matrix(rnorm(d*d), d, d)
    Q <- qr.Q(qr(M))
    if (det(Q) < 0) Q[, 1] <- -Q[, 1]
    Q
  }

  ## Multiplicative SPD-preserving noise (mean-preserving)
  # Y = X^{1/2} S X^{1/2}, S = U diag(xi) U^T, xi ~ LogNormal(mu, sigma).
  # Mean-log mu = -sigma^2/2 ensures E[xi] = 1 (mean-preserving property).
  makeY <- function(X, U, sigma_pert = 0.2) {
    d <- 3
    sigma <- sigma_pert/10          # log-scale SD (small perturbation)
    mu <- - sigma^2 / 2             # mean-log ensuring E[xi_k] = 1
    xi <- rlnorm(d, meanlog = mu, sdlog = sigma)
    S  <- U %*% diag(xi) %*% t(U)  # random SPD noise matrix
    Xh <- mat_sqrt(X, 0.5)         # X^{1/2}
    Y  <- Xh %*% S %*% Xh          # congruence transform
    (Y + t(Y)) / 2                  # symmetrize
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
  tau <- lapply(obsGrid, function(tt) generate_spd_matrix_3d(tt, axis_variation))
  tau <- simplify2array(tau)  # converts list of 3x3 matrices to a 3x3xT array
  omega0 <- diag(0, 3)        # reference matrix (zero matrix)

  # -----------------------------------------------------------------------------
  # Warping maps H, Psi, R
  # -----------------------------------------------------------------------------

  # Subject-level H_i, H_i^{-1}: exponential warps H_i^{-1}(t) = (e^{z_i t}-1)/(e^{z_i}-1)
  m <- length(obsGrid)
  hinvmat <- hmat <- matrix(0, n, m)
  z_subj <- rnorm(n, 0, sigma_warp)
  for (i in 1:n) {
    hinvmat[i, ] <- if (abs(z_subj[i]) < .Machine$double.eps)
      obsGrid else (exp(obsGrid * z_subj[i]) - 1)/(exp(z_subj[i]) - 1)
    hmat[i, ] <- approx(hinvmat[i, ], obsGrid, obsGrid, rule = 2)$y
  }

  # Component-level Psi_j, Psi_j^{-1}: beta-CDF mixture warps
  p <- 4
  alpha <- c(2, 1)    # beta shape parameters
  beta <- c(2, 1/2)   # beta rate parameters
  lambda <- 0.5       # mixing weight
  psi <- psi_inv <- vector("list", p)
  for (j in 1:(p/2)) {
    psi[[j]] <- lambda * pbeta(obsGrid, alpha[j], beta[j]) + (1 - lambda)*obsGrid
    psi_inv[[j]] <- approx(psi[[j]], obsGrid, obsGrid, rule = 2)$y
    # Reflected: Psi_{j+2}^{-1}(t) = 2t - Psi_j^{-1}(t)
    psi_inv[[j+2]] <- 2*obsGrid - psi_inv[[j]]
    psi[[j+2]] <- approx(psi_inv[[j+2]], obsGrid, obsGrid, rule = 2)$y
  }

  # R_{i,j} and R_{i,j}^{-1}: nuisance within-component warps (independent per (i,j))
  rmat <- rinv <- vector("list", n)
  for (i in 1:n) {
    rmat[[i]] <- vector("list", p)
    rinv[[i]] <- vector("list", p)
    for (j in 1:p) {
      r <- rnorm(1, 0, sigma_dist)
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
  # NOTE: Phase-only case — all amplitudes fixed at 1 (no amplitude variation).
  # -----------------------------------------------------------------------------
  amat <- matrix(1, n, p)  # all A_{ij} = 1 (constant; no amplitude variation)
  
  # -----------------------------------------------------------------------------
  # Generate noiseless SPD trajectories X_{i,j}(t) for t in obsGrid
  # -----------------------------------------------------------------------------

  # Build a single noiseless SPD trajectory: X_{ij}(t) = omega0 + A_ij*(tau(s(t)) - omega0)
  # where s(t) = R_{ij}(Psi_j(H_i(t))) is the composed warp.
  make_X <- function(A_ij, omega0, rmat_ij, h_row, psi_j, obsGrid) {
    d <- nrow(omega0)
    arr <- array(0, c(d,d,length(obsGrid)))
    # s = (R_{i,j} ∘ Psi_j ∘ H_i)(t) for t in obsGrid
    s <- compose_map(
      rmat_ij, obsGrid,
      compose_map(psi_j, obsGrid, h_row)
    )
    for (k in seq_along(obsGrid)) {
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

  # Fully aligned (oracle): X_{ij}(G_{ij}^{-1}(t)) = omega0 + A_ij*(tau(t) - omega0)
  make_X_aligned_array <- function(A_ij, omega0, obsGrid) {
    d <- nrow(omega0)
    arr <- array(0, c(d,d,length(obsGrid)))
    for (k in seq_along(obsGrid)) {
      arr[,,k] <- A_ij * ( generate_spd_matrix_3d(obsGrid[k], axis_variation) - omega0 ) + omega0
    }
    arr
  }

  # Partially aligned (nuisance removed only): X_{ij} at R_{ij}^{-1}(t)
  make_X_Ginv_array <- function(A_ij, omega0, rmat_ij, obsGrid) {
    d <- nrow(omega0)
    arr <- array(0, c(d,d,length(obsGrid)))
    for (k in seq_along(obsGrid)) {
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
  U <- rand_orth(3)   # shared orthogonal matrix for all noise realizations
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