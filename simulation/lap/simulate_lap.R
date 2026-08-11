# =============================================================================
# File: simulate_lap.R
# Purpose: Generate Laplacian-valued functional data for GLDM simulation study
#          with both amplitude and phase variation (the "bothvar" setting).
# Dependencies: none (self-contained helpers defined inline)
# =============================================================================

# -----------------------------------------------------------------------------
# Template generation: tau(t) as a Laplacian trajectory on a fixed node set
# -----------------------------------------------------------------------------

#' Generate the latent template Laplacian trajectory tau(t)
#'
#' @description
#' Builds a d x d x T array representing the template Laplacian trajectory
#' tau(t) on the given obsGrid. The edge-weight function Xi_{uv}(t) is
#' determined by a deterministic rule based on node-index parity, and only
#' active edges (mask = 1) contribute. The mask M is drawn once from
#' Bernoulli(p_edge) on the upper triangle and symmetrized.
#'
#' @param obsGrid Numeric vector; time grid on [0, 1] (length T).
#' @param nnode   Integer; number of graph nodes (d = nnode).
#' @param p_edge  Numeric in (0, 1); Bernoulli edge probability for the mask.
#' @param seed    Integer or NULL; optional RNG seed for reproducibility.
#'
#' @return A list with:
#'   - tau: d x d x T array of Laplacian matrices on obsGrid
#'   - obsGrid: the input time grid (returned for convenience)
#'   - mask: d x d symmetric binary edge-presence matrix
make_tau <- function(
    obsGrid = seq(0, 1, 0.05),  # time grid on [0,1]
    nnode   = 10,               # number of nodes (fixed across runs)
    p_edge  = 0.5,              # Bernoulli prob. for the symmetric mask
    seed    = NULL              # optional seed for reproducibility
) {
  if (!is.null(seed)) set.seed(seed)
  Tn <- length(obsGrid)
  d  <- nnode
  
  # --- (1) Symmetric binary mask M ~ Bernoulli(p_edge) on upper triangle (u<v) ---
  M <- matrix(0L, d, d)
  idx <- which(upper.tri(M), arr.ind = TRUE)
  M[idx] <- rbinom(n = nrow(idx), size = 1, prob = p_edge)
  M <- M + t(M)  # symmetrize; diagonal remains 0
  
  # --- (2) Edge-weight trajectory Xi(t) per t, then Laplacian L(t) = diag(Xi 1) - Xi ---
  # Build Xi(t) given t in [0,1] and mask M
  Xi_eval <- function(u, v, t, mask) {
    
    if (mask[u, v] == 0) {
      return(0)
    }
    
    if (u > v) {
      return(Xi_eval(v, u, t, mask))
    }
    
    # case 1: (u+v) is odd
    if ((u + v) %% 2 == 1) {
      return(0.1 + 0.9 * t)
    }
    
    # case 2: both are odd
    if ((u %% 2 == 1) && (v %% 2 == 1)) {
      return(1 / (1 + 9 * t))
    }
    
    # case 3: both are even
    if ((u %% 2 == 0) && (v %% 2 == 0)) {
      x <- t  
      return(0.4 * sin(pi * (x - 0.5)) + 0.6)
    }
    return(0)
  }
  
  build_Xi_at_t <- function(t, M) {
    d  <- nrow(M)
    Xi <- matrix(0, d, d)
    
    # Compute only for upper-triangular entries
    for (u in 1:(d - 1)) {
      for (v in (u + 1):d) {
        
        # Evaluate Xi 
        val <- Xi_eval(u, v, t, M)
        Xi[u, v] <- val
        Xi[v, u] <- val   # enforce symmetry
      }
    }
    
    Xi
  }
  # Laplacian map: L = diag(Xi * 1) - Xi
  laplacian_map <- function(Xi) {
    deg <- rowSums(Xi)
    diag(deg, d) - Xi
  }
  
  # --- (3) Build tau(t) on obsGrid and stack to a 3D array (d x d x T) ---
  tau_list <- vector("list", Tn)
  for (k in seq_len(Tn)) {
    Xi_t     <- build_Xi_at_t(obsGrid[k],M)
    tau_list[[k]] <- laplacian_map(Xi_t)
  }
  tau_arr <- simplify2array(tau_list)  # d x d x T
  
  # Return everything needed downstream
  list(
    tau     = tau_arr,   # d x d x T Laplacian trajectory
    obsGrid = obsGrid,   # time grid
    mask    = M          # fixed symmetric edge mask used to build Xi(t)
  )
}

# -----------------------------------------------------------------------------
# Top-level helper functions (also used by simulate_from_tau)
# -----------------------------------------------------------------------------

#' Evaluate edge-weight Xi_{uv}(t) for a given node pair and time point
#'
#' @description
#' Returns the scalar edge weight at time t for the (u, v) pair according to
#' a parity-based deterministic rule. The function enforces u <= v (symmetric),
#' and returns 0 for masked-out edges (mask[u,v] == 0).
#'
#' @param u    Integer; first node index.
#' @param v    Integer; second node index.
#' @param t    Numeric scalar; time point in [0, 1].
#' @param mask d x d binary matrix; 1 if edge (u, v) is active, 0 otherwise.
#'
#' @return Numeric scalar; edge weight Xi_{uv}(t).
Xi_eval <- function(u, v, t, mask) {

  if (mask[u, v] == 0) {
    return(0)  # edge absent in this graph realization
  }

  if (u > v) {
    return(Xi_eval(v, u, t, mask))  # enforce canonical u < v order
  }

  # case 1: (u+v) is odd — linearly increasing edge weight
  if ((u + v) %% 2 == 1) {
    return(0.1 + 0.9 * t)
  }

  # case 2: both u and v are odd — hyperbolically decreasing edge weight
  if ((u %% 2 == 1) && (v %% 2 == 1)) {
    return(1 / (1 + 9 * t))
  }

  # case 3: both u and v are even — sinusoidal edge weight
  if ((u %% 2 == 0) && (v %% 2 == 0)) {
    x <- t
    return(0.4 * sin(pi * (x - 0.5)) + 0.6)
  }
  return(0)
}

#' Build the edge-weight matrix Xi at a single time point t
#'
#' @description
#' Evaluates Xi_eval for all active upper-triangle pairs and symmetrizes.
#'
#' @param t Numeric scalar; time point in [0, 1].
#' @param M d x d binary symmetric edge mask.
#'
#' @return d x d numeric matrix of edge weights at time t (symmetric, diagonal 0).
build_Xi_at_t <- function(t, M) {
  d  <- nrow(M)
  Xi <- matrix(0, d, d)

  # Compute only for upper-triangular entries; lower triangle is set by symmetry
  for (u in 1:(d - 1)) {
    for (v in (u + 1):d) {

      # Evaluate Xi
      val <- Xi_eval(u, v, t, M)
      Xi[u, v] <- val
      Xi[v, u] <- val   # enforce symmetry
    }
  }

  Xi
}

#' Convert an edge-weight matrix to a graph Laplacian
#'
#' @description
#' Computes L = diag(A * 1) - A for a symmetric nonneg-edge matrix A.
#'
#' @param A d x d symmetric matrix of edge weights (nonneg, zero diagonal).
#'
#' @return d x d graph Laplacian matrix (positive semidefinite, zero row sums).
laplacian_map <- function(A) {
  deg <- rowSums(A)  # degree of each node
  diag(deg, nrow(A)) - A
}

#' Evaluate tau(t) exactly (no grid interpolation)
#'
#' @description
#' Returns the d x d Laplacian matrix at time t by evaluating the edge-weight
#' function analytically and applying the Laplacian map.
#'
#' @param t Numeric scalar; time point in [0, 1].
#' @param M d x d binary symmetric edge mask.
#'
#' @return d x d Laplacian matrix.
# Exact tau(s) evaluator (NO interpolation)
tau_at_exact <- function(t, M) {
  laplacian_map(build_Xi_at_t(t, M))
}

# -----------------------------------------------------------------------------
# Data simulation: generate observed Laplacian-valued trajectories Y_{ij}(t)
# -----------------------------------------------------------------------------

#' Simulate Laplacian-valued functional data from a fixed template tau
#'
#' @description
#' Generates n subjects x p component Laplacian-valued trajectories following
#' the GLDM generative model: X_{ij}(t) = A_{ij} * (tau(G_{ij}(t)) - omega0)
#' + omega0, where G_{ij} = R_{ij} o Psi_j o H_i is the global deformation
#' (nuisance R_{ij}, component Psi_j, subject H_i). Amplitude factors A_{ij}
#' are drawn from Uniform(1, 2). Observed Y_{ij}(t) adds edge-wise Beta(2,2)
#' perturbation noise (controlled by sigma_pert).
#'
#' @param tau       d x d x T array; the fixed Laplacian template on obsGrid.
#' @param obsGrid   Numeric vector of length T; time grid where tau is defined.
#' @param mask      d x d binary symmetric matrix; edge mask from make_tau().
#' @param n         Integer; number of subjects.
#' @param p         Integer; number of components (must be even).
#' @param sigma_warp Numeric; SD of subject-level warp parameter z_i ~ N(0, sigma_warp).
#' @param sigma_dist Numeric; SD of nuisance warp parameter r_{ij} ~ N(0, sigma_dist).
#' @param sigma_pert Numeric; amplitude of Beta(2,2) edge-weight perturbation.
#' @param seed       Integer or NULL; optional RNG seed for reproducibility.
#'
#' @return A list with all ground-truth components needed for evaluation:
#'   - obsGrid, tau, Leta: template and component tempo trajectories
#'   - psi, psi_inv: component-level warps and inverses (length p)
#'   - hmat, hinvmat: subject-level forward and inverse warps (n x T matrices)
#'   - rmat, rinv: nuisance warps per (i, j) (lists of length n, each length p)
#'   - amat: n x p amplitude factors
#'   - X: noiseless Laplacian trajectories (list of p lists of n arrays)
#'   - Y: noisy observed trajectories (list of p lists of n arrays)
#'   - X_Ginv: X evaluated at G_{ij}^{-1}(t) (list p x n)
#'   - X_aligned: fully aligned X (list p x n)
simulate_from_tau <- function(
    tau,                  # d x d x T array: fixed template trajectory on obsGrid
    obsGrid,              # vector: grid where 'tau' is defined
    mask,                 # binary mask used in tau generation
    n = 15,               # number of subjects
    p = 4,                # number of components (must be even for Psi construction)
    sigma_warp = 0.5,     # subject-level deformation parameter
    sigma_dist = 0.5,     # nuisance time distortion parameter
    sigma_pert = 0.5,     # edge-weight perturbation noise level
    seed = NULL           # optional RNG seed for reproducibility
) {
  if (!is.null(seed)) set.seed(seed)
  deg = dim(tau[,,1])[1]
  omega0 = diag(0,deg) # reference origin (zero matrix of size d x d)
  d <- nrow(omega0)
  m <- length(obsGrid)

  # --- Subject-level H_i and inverse H_i^{-1} ---
  # H_i is an exponential warp: H_i^{-1}(t) = (exp(z_i * t) - 1) / (exp(z_i) - 1)
  hmat <- hinvmat <- matrix(0, n, m)
  z_subj <- rnorm(n, 0, sigma_warp)
  for (i in 1:n) {
    hinvmat[i, ] <- if (abs(z_subj[i]) < .Machine$double.eps)
      obsGrid else (exp(obsGrid * z_subj[i]) - 1) / (exp(z_subj[i]) - 1)
    hmat[i, ] <- approx(hinvmat[i, ], obsGrid, obsGrid, rule = 2)$y  # numeric inverse
  }

  # --- Component-level Psi_j and inverse Psi_j^{-1} ---
  # Psi_j is a Beta-CDF mixture: lambda * F_{Beta}(t) + (1-lambda) * t
  # For j in {3,4}: Psi_j^{-1} = 2t - Psi_{j-2}^{-1}(t) (reflective symmetry)
  alpha <- c(2, 1)
  beta  <- c(2, 1/2)
  lambda <- 0.5
  psi <- psi_inv <- vector("list", p)
  for (j in 1:(p/2)) {
    psi[[j]]     <- lambda * pbeta(obsGrid, alpha[j], beta[j]) + (1 - lambda) * obsGrid
    psi_inv[[j]] <- approx(psi[[j]], obsGrid, obsGrid, rule = 2)$y
    psi_inv[[j+2]] <- 2*obsGrid - psi_inv[[j]]  # reflected inverse for partner component
    psi[[j+2]]   <- approx(psi_inv[[j+2]], obsGrid, obsGrid, rule = 2)$y
  }

  # --- Nuisance R_{ij} and inverse ---
  # R_{ij} is an exponential warp with random parameter r_{ij} ~ N(0, sigma_dist)
  rmat <- rinv <- vector("list", n)
  for (i in 1:n) {
    rmat[[i]] <- rinv[[i]] <- vector("list", p)
    for (j in 1:p) {
      r <- rnorm(1, 0, sigma_dist)
      rinv_ij <- if (abs(r) < .Machine$double.eps)
        obsGrid else (exp(obsGrid * r) - 1) / (exp(r) - 1)
      rinv[[i]][[j]] <- rinv_ij
      rmat[[i]][[j]] <- approx(rinv_ij, obsGrid, obsGrid, rule = 2)$y  # numeric inverse
    }
  }

  # --- Component tempo trajectories: Leta_j(t) = tau(Psi_j(t)) ---
  # Stored for evaluation purposes; represents the true component-level tempo
  Leta <- vector("list", p)
  for (j in seq_len(p)) {
    Tn <- length(obsGrid)
    arr <- array(0, c(nrow(tau[,,1]), nrow(tau[,,1]), Tn))
    for (k in seq_len(Tn)) {
      s <- psi[[j]][k]                # s = Psi_j(t_k)
      arr[,,k] <- tau_at_exact(s, mask)
    }
    Leta[[j]] <- arr
  }

  # --- Amplitudes a_{ij} ~ Uniform(1, 2) (amplitude + phase setting) ---
  amat <- matrix(runif(n * p, 1, 2), n, p)
  
  # --- Generate noiseless trajectories X_{ij}(t) ---
  # Internal helper: build X_{ij}(t) = A_ij*(tau(G_{ij}(t)) - omega0) + omega0
  # where G_{ij} = R_{ij} o Psi_j o H_i
  make_X<- function(A_ij, omega0, rmat_ij, h_row, psi_j, obsGrid, mask) {
    d   <- nrow(omega0)
    Tn  <- length(obsGrid)
    arr <- array(0, c(d, d, Tn))
    # s(t) = (R_{ij} ∘ Psi_j ∘ H_i)(t): compose all three warps at each time point
    s_vals <- approx(obsGrid, rmat_ij, xout = approx(obsGrid, psi_j, xout = h_row, rule = 2)$y, rule = 2)$y
    for (k in seq_len(Tn)) {
      Tau_s   <- tau_at_exact(s_vals[k], mask)  # evaluate template at warped time
      arr[,,k] <- A_ij * (Tau_s - omega0) + omega0  # amplitude scaling
    }
    arr
  }
  
  simdata_X <- vector("list", p)
  for (j in 1:p) {
    simdata_X[[j]] <- vector("list", n)
    for (i in 1:n) {
      simdata_X[[j]][[i]] <- make_X(amat[i, j], omega0, rmat[[i]][[j]], hmat[i, ], psi[[j]], obsGrid, mask)
    }
  }
  
  # --- "True aligned" variants for analysis ---
  # X_aligned: X(G_{ij}^{-1}(R_{ij}^{-1}(t))) = a_ij * (tau(t) - omega0) + omega0
  make_X_aligned <- function(A_ij, omega0, obsGrid, mask) {
    d <- nrow(omega0); Tn <- length(obsGrid)
    arr <- array(0, c(d, d, Tn))
    for (k in seq_len(Tn)) {
      Tau_t    <- tau_at_exact(obsGrid[k], mask)
      arr[,,k] <- A_ij * (Tau_t - omega0) + omega0
    }
    arr
  }
  
  # X_Ginv: X(G_{ij}^{-1}(t)) = a_ij * (tau(R_{ij}(t)) - omega0) + omega0
  make_X_Ginv <- function(A_ij, omega0, rmat_ij, obsGrid, mask) {
    d <- nrow(omega0); Tn <- length(obsGrid)
    arr <- array(0, c(d, d, Tn))
    for (k in seq_len(Tn)) {
      Tau_r    <- tau_at_exact(rmat_ij[k], mask)
      arr[,,k] <- A_ij * (Tau_r - omega0) + omega0
    }
    arr
  }
  
  X_Ginv   <- vector("list", p)
  X_aligned <- vector("list", p)
  for (j in seq_len(p)) {
    X_Ginv[[j]]    <- vector("list", n)
    X_aligned[[j]] <- vector("list", n)
    for (i in seq_len(n)) {
      X_Ginv[[j]][[i]]    <- make_X_Ginv(amat[i, j], omega0, rmat[[i]][[j]], obsGrid, mask)
      X_aligned[[j]][[i]] <- make_X_aligned(amat[i, j], omega0, obsGrid, mask)
    }
  }
  
  # ================================
  # Laplacian <-> Edge space helpers and perturbation map P (local scope)
  # ================================
  # Local redefinition ensures these helpers are available in simulate_from_tau's
  # lexical scope even when called in isolation.
  laplacian_map <- function(A) {
    # L = diag(A 1) - A, for symmetric A with zero diagonal and nonnegative edges
    deg <- rowSums(A)
    diag(deg, nrow(A)) - A
  }
  laplacian_inv <- function(L) {
    # Recover edge weights from Laplacian: A_{uv} = -L_{uv} (u!=v), A_{uu}=0
    A <- -L
    diag(A) <- 0
    # Ensure symmetry and nonnegativity (numerical safety)
    A <- (A + t(A)) / 2
    A[A < 0] <- 0
    A
  }
  apply_P_edgewise <- function(L, sigma_pert) {
    # If sigma_pert = 0 -> identity map
    if (sigma_pert <= 0) return(L)
    A <- laplacian_inv(L)
    # Add symmetric Beta(2,2)-centered jitter ONLY on existing edges (A != 0)
    for (u in 1:(nrow(A) - 1)) {
      for (v in (u + 1):ncol(A)) {
        if (A[u, v] != 0) {
          eps <- rbeta(1, 2, 2)         # in (0,1); mean = 0.5
          delta <- 0.1 * sigma_pert * (eps - 0.5)  # centered in ~[-0.05*sigma_pert, 0.05*sigma_pert]
          val <- A[u, v] + delta
          # Keep valid (nonnegative) edge weights
          val <- max(val, 0)
          A[u, v] <- A[v, u] <- val
        }
      }
    }
    # Rebuild Laplacian
    laplacian_map(A)
  }
  
  # ================================
  # Add perturbation to get Y via P_{i,j,k} in Laplacian/edge space
  # ================================
  simdata_Y <- vector("list", p)
  for (j in 1:p) {
    simdata_Y[[j]] <- vector("list", n)
    for (i in 1:n) {
      Tn <- dim(simdata_X[[j]][[i]])[3]
      arrY <- array(0, c(d, d, Tn))
      for (k in 1:Tn) {
        # X[[j]][[i]][,,k] is a Laplacian already (by construction from tau)
        Lx <- simdata_X[[j]][[i]][,,k]
        arrY[,,k] <- apply_P_edgewise(Lx, sigma_pert)
      }
      simdata_Y[[j]][[i]] <- arrY
    }
  }
  
  
  # --- Return package ---
  list(
    # Grids and template
    obsGrid   = obsGrid,
    tau       = tau,
    Leta      = Leta,
    # Warping maps
    psi       = psi,
    psi_inv   = psi_inv,
    hmat      = hmat,
    hinvmat   = hinvmat,
    rmat      = rmat,
    rinv      = rinv,
    # Parameters
    amat      = amat,
    # Trajectories
    X         = simdata_X,
    Y         = simdata_Y,
    # Aligned versions
    X_Ginv    = X_Ginv,
    X_aligned = X_aligned
  )
}
