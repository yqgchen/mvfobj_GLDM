# =============================================================================
# File: fctns_estimation.R
# Purpose: Core estimation functions for the GLDM applied to SPD matrix-valued
#          functional data. Implements amplitude normalization, subject-level
#          and component-level warping via Frobenius pairwise warping, latent
#          template estimation, and the full GLDM fitting pipeline.
# Dependencies: frechet (LocCovReg, CovFMean), minqa (bobyqa), SMFilter (FDist2),
#               pracma (trapz), stats (approx, optim)
# =============================================================================

# -----------------------------------------------------------------------------
# Low-level accessor and utility helpers
# -----------------------------------------------------------------------------

#' Extract a single d x d matrix slice from a component-first SPD list
#'
#' @param X Component-first list: X[[j]][[i]] is a d x d x T array.
#' @param i Integer; subject index.
#' @param j Integer; component index.
#' @param t Integer; time index (third dimension).
#' @return A d x d matrix.
get_mat <- function(X, i, j, t) X[[j]][[i]][,,t, drop = FALSE][,,1]

#' Set a single d x d matrix slice in a component-first SPD list
#'
#' @param X Component-first list to modify.
#' @param i Integer; subject index.
#' @param j Integer; component index.
#' @param t Integer; time index.
#' @param M d x d matrix to insert.
#' @return The modified list X.
set_mat <- function(X, i, j, t, M) { X[[j]][[i]][,,t] <- M; X }

#' Get the number of subjects from a component-first SPD list
#'
#' @param X Component-first list: X[[j]][[i]] is a d x d x T array.
#' @return Integer; n = length(X[[1]]).
n_subjects <- function(X) length(X[[1]])

#' Get the number of components from a component-first SPD list
#'
#' @param X Component-first list: X[[j]][[i]] is a d x d x T array.
#' @return Integer; p = length(X).
n_components <- function(X) length(X)

#' Get the (d, d, T) dimensions of the arrays in a component-first SPD list
#'
#' @param X Component-first list; X[[1]][[1]] must be a d x d x T array.
#' @return Integer vector of length 3: c(d, d, T).
dims_of <- function(X) c(dim(X[[1]][[1]])[1:2], # d, d (e.g., matrix dimensions)
                         dim(X[[1]][[1]])[3])    # T (e.g., number of time points)

#' Frobenius distance from a d x d matrix to a reference matrix
#'
#' @description
#' Computes ||X - omega0||_F using base::norm(..., "F").
#'
#' @param X      d x d matrix.
#' @param omega0 d x d reference matrix; defaults to zero matrix if NULL.
#' @return Non-negative scalar: Frobenius distance.
frob_dist <- function(X, omega0 = NULL) {
  stopifnot(is.matrix(X))
  if (is.null(omega0)) omega0 <- matrix(0, nrow(X), ncol(X))
  base::norm(X - omega0, type = "F")
}

#' Compute sup_t ||Y_{ij}(t) - omega0||_F for one (i,j) trajectory
#'
#' @description
#' Given a d x d x T array Y, returns the maximum Frobenius distance
#' from omega0 over all time indices t = 1..T.
#'
#' @param Y      d x d x T array (SPD trajectory for one subject-component pair).
#' @param omega0 d x d reference matrix; defaults to zero matrix if NULL.
#' @return Non-negative scalar: max_t ||Y[,,t] - omega0||_F.
normA_one <- function(Y, omega0 = NULL) {
  stopifnot(length(dim(Y)) == 3, dim(Y)[1] == dim(Y)[2])
  d <- dim(Y)[1]; Tn <- dim(Y)[3]
  if (is.null(omega0)) omega0 <- matrix(0, d, d)
  dists <- vapply(seq_len(Tn), function(tt) frob_dist(Y[,,tt, drop = FALSE][,,1], omega0), numeric(1))
  max(dists, na.rm = TRUE)
}

#' Compute sup-norm amplitudes for all subjects in one component
#'
#' @description
#' For component j, applies normA_one to each subject i = 1..n and
#' returns a vector of length n.
#'
#' @param X      Component-first list: X[[j]][[i]] is a d x d x T array.
#' @param j      Integer; component index.
#' @param omega0 d x d reference matrix; defaults to zero matrix if NULL.
#' @return Numeric vector of length n: sup-norm amplitudes for all subjects.
normA_component <- function(X, j, omega0 = NULL) {
  n <- n_subjects(X)
  vapply(seq_len(n), function(i) normA_one(X[[j]][[i]], omega0), numeric(1))
}


# -----------------------------------------------------------------------------
# Amplitude estimation and normalization
# -----------------------------------------------------------------------------

#' Compute amplitude factors A_{i,j} for SPD/PSD trajectories
#'
#' @description
#' Given a component-first SPD/PSD trajectory object \code{X} such that
#' \code{X[[j]][[i]]} is a \eqn{d×d×T} matrix-valued trajectory for subject \eqn{i}
#' and component \eqn{j}, compute amplitude factors \eqn{A_{i,j}} measuring
#' the maximal deviation (e.g., Frobenius distance) from a reference matrix \eqn{w_0},
#' normalized by a global constant \code{Const}.
#'
#' The resulting amplitude matrix \eqn{A} is often used to form amplitude-normalized
#' trajectories:
#' \deqn{X^{*}_{i,j}(t) = w_0 + (X_{i,j}(t) - w_0) / A_{i,j}.}
#'
#' @param X  List of length \eqn{p} (components);
#'           each \code{X[[j]]} is a list of length \eqn{n} (subjects);
#'           each \code{X[[j]][[i]]} is a \eqn{d×d×T} SPD/PSD array.
#' @param omega0 Reference \eqn{d×d} matrix (the "center" of scaling).
#'        If \code{NULL}, defaults to a zero matrix.
#'
#' @return A numeric matrix \eqn{A} of size \eqn{n×p},
#'         where entry \eqn{A_{i,j}} is the amplitude factor for subject \eqn{i}
#'         and component \eqn{j}.
#'
#' @examples
#' # A <- get_amp(X, omega0 = diag(0, d))
#' 
#' @export
get_amp <- function(X, omega0 = NULL) {
  n <- n_subjects(X); p <- n_components(X)
  d <- dim(X[[1]][[1]])[1]
  if (is.null(omega0)) omega0 <- matrix(0, d, d)
  Apre <- matrix(NA_real_, nrow = n, ncol = p)
  for (j in seq_len(p)) Apre[, j] <- normA_component(X, j, omega0)
  c0 <- min(Apre) # c0hat = min_{i,j} sup_t ||Y_ij(t) - omega0||_F
  Apre/c0
}

#' Normalize SPD/PSD trajectories by removing amplitude factors
#'
#' @description
#' Given a list of d×d×T arrays Y[[j]][[i]] (subjects i = 1..n, components j = 1..p), reference matrix omega0, 
#' and an amplitude matrix A (n×p), return trajectories normalized by
#'   Y_norm = omega0 + (Y - omega0) / A[i, j].
#' This removes the per-(i,j) amplitude factors.
#'
#' @param Y  List: length p; each Y[[j]] is a list of length n;
#'           each Y[[j]][[i]] is a d×d×T array (SPD/PSD trajectory).
#' @param A  Numeric matrix (n×p) of amplitudes (e.g., from get_amp()).
#' @param omega0 Reference d×d matrix; if NULL, defaults to a zero matrix of size d.
#'
#' @return A list with the same structure as Y, containing amplitude-normalized arrays.
#' @examples
#' # Y_norm <- normalize(Y = Y, A = A, omega0 = diag(0, d))
#' @export
normalize <- function(Y, A, omega0 = NULL) {
  # ---- validate structure ----
  if (!is.list(Y) || length(Y) == 0L)
    stop("Y must be a non-empty list: Y[[j]][[i]] is a d×d×T array.")
  p <- length(Y)
  if (!is.list(Y[[1]]) || length(Y[[1]]) == 0L)
    stop("Each Y[[j]] must be a non-empty list of subjects.")
  
  n <- length(Y[[1]])
  if (!is.matrix(A) || any(dim(A) != c(n, p)))
    stop(sprintf("A must be an n×p matrix with dims (%d, %d).", n, p))
  
  ref <- Y[[1]][[1]]
  
  if (!is.array(ref) || length(dim(ref)) != 3L)
    stop("Each Y[[j]][[i]] must be a 3D array (d×d×T).")
  
  d    <- dim(ref)[1]
  Tlen <- dim(ref)[3]
  if (dim(ref)[1] != dim(ref)[2])
    stop("Arrays must be square (d×d×T).")
  
  # check all arrays share dims
  for (j in seq_len(p)) {
    if (length(Y[[j]]) != n)
      stop("All components must have the same number of subject n.")
    for (i in seq_len(n)) {
      a <- Y[[j]][[i]]
      if (!is.array(a) || length(dim(a)) != 3L ||
          any(dim(a)[1:2] != d) || dim(a)[3] != Tlen) {
        stop("All Y[[j]][[i]] must be d×d×T with identical dimensions.")
      }
    }
  }
  # check A
  if ( any( !is.finite(A) ) ) {
    stop("A includes infinite values.")
  }
  if ( any( !(A>0) ) ) {
    stop("A includes nonpositive values.")
  }
  
  # default omega0 = 0_dxd
  if (is.null(omega0)) omega0 <- matrix(0, d, d)
  if (!is.matrix(omega0) || any(dim(omega0) != c(d, d)))
    stop("omega0 must be a d×d matrix matching Y.")
  
  # ---- normalize by amplitude per (i,j) ----
  Y_norm <- vector("list", p)
  for (j in seq_len(p)) {
    Y_norm[[j]] <- vector("list", n)
    for (i in seq_len(n)) {
      a_ij <- A[i, j]
      
      out <- array(0, dim = c(d, d, Tlen))
      for (t in seq_len(Tlen)) {
        # center-preserving scaling: omega0 + a_ij*(Y - omega0)
        M <- omega0 + (1/a_ij) * (Y[[j]][[i]][, , t] - omega0)
        # symmetrize for numerical stability (keeps PSD under mild roundoff)
        #out[, , t] <- (M + t(M)) / 2
        out[, , t] <- M
      }
      Y_norm[[j]][[i]] <- out
    }
  }
  
  Y_norm
}


# -----------------------------------------------------------------------------
# Frechet mean and lambda utilities
# -----------------------------------------------------------------------------

#' Fréchet mean over a list at each time (Frobenius metric)
#'
#' @description
#' Given a list of d×d×T arrays (all with identical dimensions),
#' compute the Fréchet mean (with Frobenius metric) at each time index k,
#' returning a single d×d×T array as the mean trajectory.
#'
#' @param Y_list A list whose elements are d×d×T arrays (SPD trajectories).
#'               All arrays must share the same dimensions.
#' @return A d×d×T array: the Fréchet mean (Frobenius) at each time.
#' @importFrom frechet CovFMean
#' @examples
#' # Y_list <- list(array1, array2, array3)
#' # mean_traj <- Fmean_list(Y_list)
Fmean_list <- function(Y_list) {
  if (!is.list(Y_list) || length(Y_list) < 1)
    stop("Y_list must be a non-empty list of d×d×T arrays.")
  
  ref <- Y_list[[1]]
  if (!is.array(ref) || length(dim(ref)) != 3)
    stop("Each element of Y_list must be a 3D array (d×d×T).")
  
  d <- dim(ref)[1]
  Tlen <- dim(ref)[3]
  if (dim(ref)[1] != dim(ref)[2])
    stop("Arrays must be square (d×d×T).")
  
  for (i in seq_along(Y_list)) {
    a <- Y_list[[i]]
    if (!is.array(a) || length(dim(a)) != 3 || any(dim(a)[1:2] != d) || dim(a)[3] != Tlen)
      stop("All arrays must be d×d×T and share the same dimensions.")
  }
  
  L <- length(Y_list)
  mean_arr <- array(0, dim = c(d, d, Tlen))
  
  for (k in seq_len(Tlen)) {
    # Stack slices into a d×d×L array
    M_k <- array(0, dim = c(d, d, L))
    for (i in seq_len(L)) M_k[, , i] <- Y_list[[i]][, , k]
    Fmean <- frechet::CovFMean(M = M_k, optns = list(metric = "frobenius"))$Mout
    mean_arr[, , k] <- Fmean[[1]]
  }
  
  mean_arr
}

#' Generate a default lambda via Integrated MSE (Frobenius) against the Fréchet mean
#'
#' @description
#' Compute IMSE = mean_i ∫ ||Y_i(t) - \bar{Y}(t)||_F^2 dt,
#' then return lambda = scale * IMSE. 
#'
#' @param Y_list A list of d×d×T arrays (same dims).
#' @param workGrid Numeric vector of length T (integration grid).
#' @param scale Numeric scalar to scale the IMSE (default 1e-4).
#' @return A numeric scalar lambda.
#' @examples
#' # lambda <- get_lambda(Y_list, workGrid, scale = 1e-4)
get_lambda <- function(Y_list, workGrid, scale = 1e-4) {
  if (!is.list(Y_list) || length(Y_list) < 1)
    stop("Y_list must be a non-empty list of d×d×T arrays.")
  ref <- Y_list[[1]]
  if (!is.array(ref) || length(dim(ref)) != 3)
    stop("Each element of Y_list must be a d×d×T array.")
  d <- dim(ref)[1]
  Tlen <- dim(ref)[3]
  if (dim(ref)[1] != dim(ref)[2])
    stop("Each array must be square (d×d×T).")
  if (length(workGrid) != Tlen)
    stop("length(workGrid) must equal T (third dimension of arrays).")
  
  # Prefer SMFilter::FDist2 if available; otherwise fallback
  FDist2_local <- if (exists("FDist2", where = asNamespace("SMFilter"), inherits = FALSE)) {
    SMFilter::FDist2
  } else {
    function(A, B) sum((A - B)^2)
  }
  
  mean_traj <- Fmean_list(Y_list)
  
  imse_vals <- vapply(seq_along(Y_list), function(i) {
    difft <- vapply(seq_len(Tlen), function(k) {
      FDist2_local(Y_list[[i]][, , k], mean_traj[, , k])
    }, numeric(1))
    pracma::trapz(workGrid, difft)
  }, numeric(1))
  
  imse <- mean(imse_vals)
  scale * imse
}


# -----------------------------------------------------------------------------
# Pairwise Frobenius warping (sparse and dense versions)
# -----------------------------------------------------------------------------

#' PSD Pairwise Warping
#'
#' This function performs pairwise warping for PSD matrix-valued random processes
#' based on the Frobenius metric.
#'
#' @param Lt A list of numeric vectors, each containing time points for the corresponding subject.
#' @param Ly A list of 3D arrays, where each array contains observed matrices for each subject at corresponding time points.
#' @param optns A list of control parameters for the function.
#'
#' @return A list containing:
#' \item{h}{A matrix of warping functions evaluated on \code{workGrid}.}
#' \item{hInv}{A matrix of inverse warping functions evaluated on \code{workGrid}.}
#' \item{yhAligned}{A list of matrices representing aligned y for each subject.}
#' \item{yh}{Estimated y from local Fréchet regression.}
#' \item{workGrid}{A vector containing the grid points used for evaluation.}
#' \item{optns}{The list of control options used.}
#'
#' @importFrom frechet LocCovReg CovFMean
#' @importFrom pracma trapz
#' @importFrom stats approx optim quantile weighted.mean
#' @importFrom minqa bobyqa
#' @export
#'

FrobPW <- function ( Lt=NULL, Ly=NULL, optns = list() ){
  
  # check input
  if ( !is.list (Lt) | any( !(sapply(Lt, is.numeric) & sapply(Lt, is.vector)) ) ) {
    stop("Missing input or incorrect format of Lt.")
  }
  nsubj <- length(Lt)
  ntime_per_subj <- sapply( Lt, length )
  
  if ( !is.null(Ly) ) {
    if ( !is.list(Ly) ) {
      stop( "Ly should be a list." )
    } else if ( any( !sapply(Ly, is.array) ) ) {
      stop( "Ly[[i]] should be a array for all plausible i.")
    } else if ( any( abs( sapply(Ly, function(x) dim(x)[3])- ntime_per_subj ) > 0 ) ) {
      stop( "Mismatched numbers of time points per subject between Lt and Ly." )
    } else if ( length(Ly) != nsubj ) {
      stop( "Mismatched numbers of trajectories/subjects between Lt and Ly. ")
    } else if ( !all( unlist( lapply( Ly, sapply, is.vector ) ) ) | !all( unlist( lapply( Ly, sapply, is.numeric ) ) ) ) {
      stop("Incorrect format of Ly.")
    }
  } else {
    stop("Ly should be input.")
  }
  
  # set up some default options
  if ( is.null(optns$ngrid) ) optns$ngrid <- 51
  if (is.null(optns$kernelReg)) optns$kernelReg <- 'epan'
  optns$metric <- 'frobenius'
  tin <- unique( sort( unlist( Lt ) ) )
  trange <- range(tin)
  tin <- lapply( Lt, function(t) ( t - trange[1] ) / diff(trange) )
  M <- optns$ngrid
  workGrid = seq( 0, 1, length.out = M ) # a grid on [0,1]
  
  ## pre-smoothing
  timingPrsm <- Sys.time()
  if ( !is.null(Ly) ) {
    yhat <- lapply(1:nsubj, function(i) {
      LocCovReg( x = matrix(tin[[i]], ncol = 1), M = Ly[[i]], xout = matrix(workGrid, ncol = 1), optns = optns )
    })
  }
  yhat <- lapply(1:nsubj, function(i) {
    abind::abind(yhat[[i]]$Mout, along = 3)  # list of matrices → 3D array
  })
  
  timingPrsm <- Sys.time() - timingPrsm
  
  # time warping
  
  res <- FrobPWdense( tVec = workGrid, yhat = yhat, optns = optns )
  
  return(list(
    h = res$h * diff(trange) + trange[1],
    hInv = res$hInv * diff(trange) + trange[1],
    yhatAligned = res$yhatAligned,
    yh = yhat,
    workGrid = res$workGrid * diff(trange) + trange[1],
    optns = res$optns,
    costs = res$costs
  ))
}

#' @title Frobenius pairwise warping for completely/densely observed PSD matrix-valued random processes
#'
#' @description Pairwise warping for PSD matrix-valued random processes
#' \eqn{\{X_i(\cdot)\}_{i=1}^{n}} endowed with the Frobenius metric,
#' where each \eqn{X_i(\cdot)} is a PSD matrix-valued random process densely observed
#' on the same grid across \eqn{i} (and hence no local Fréchet regression is involved).
#'
#' @param tVec A vector holding the common equidistant time grid on which all \eqn{\{X_i(\cdot)\}_{i=1}^{n}} are observed.
#' @param yhat A list of 3D arrays with \code{yhat[[i]]} holding the \eqn{i}-th PSD matrix-valued
#' trajectory \eqn{X_i(\cdot)} on \code{tVec}, of which each slice \code{yhat[[i]][,,j]} corresponds to a time point in \code{tVec}.
#' @param optns A list of control parameters specified by \code{list(name=value)},
#' including options for optimization and smoothing.
#'
#' @return A list of the following:
#' \item{h}{A matrix holding the warping functions evaluated on \code{workGrid}, each row corresponding to a subject.}
#' \item{hInv}{A matrix holding the inverse warping functions evaluated on \code{workGrid}, each row corresponding to a subject.}
#' \item{yhAligned}{A list of 3D arrays. The \eqn{i}-th array holds the \eqn{i}-th aligned
#' PSD matrix-valued trajectory on \code{workGrid}, where each slice \code{yhAligned[[i]][,,j]} corresponds to one time point in \code{workGrid}.}
#' \item{workGrid}{A copy of \code{tVec}. An equidistant grid on which
#' (inverse) warping functions in \code{h} and \code{hInv}
#' and PSD matrix-valued trajectories in \code{yhat} and \code{yhAligned} are evaluated.}
#' \item{optns}{Control options used.}
#' \item{costs}{The mean cost associated with each trajectory.}
#' \item{timingWarp}{The time cost of the warping process.}
#'
#' @importFrom pracma trapz
#' @importFrom stats approx optim weighted.mean
#' @importFrom utils installed.packages
#' @importFrom minqa bobyqa
#' @importFrom SMFilter FDist2
#' @export
FrobPWdense <- function(tVec, yhat, optns = list()) {
  ## ------------------------------------------------------------
  ## 0. Set default options
  ## ------------------------------------------------------------
  if (is.null(optns$nknots))   optns$nknots  <- 2
  if (is.null(optns$choice))   optns$choice  <- "truncated"  
  if (!optns$choice %in% c("unweighted", "weighted", "truncated")) {
    stop("optns$choice must be one of 'unweighted', 'weighted', or 'truncated'.")
  }
  if (is.null(optns$isPWL))    optns$isPWL   <- TRUE
  if (is.null(optns$seed))     optns$seed    <- 666
  if (is.null(optns$verbose))  optns$verbose <- FALSE
  
  ## lambda: use SPD-based scaling (user-defined get_lambda)
  if (is.null(optns$lambda)) {
    optns$lambda <- get_lambda(yhat, tVec)
  }
  lambda <- optns$lambda
  
  ## ------------------------------------------------------------
  ## 1. Normalize the time grid to [0,1]
  ## ------------------------------------------------------------
  trange   <- range(tVec)
  workGrid <- (tVec - trange[1]) / diff(trange)  
  M        <- length(tVec)
  nsubj    <- length(yhat)
  numOfKcurves <- nsubj - 1
  
  ## Storage containers
  gijMat      <- array(dim = c(numOfKcurves, M, nsubj))  # g_{j→i}
  distMat     <- matrix(nrow = nsubj, ncol = numOfKcurves)
  hMat        <- array(dim = c(nsubj, M))
  hInvMat     <- array(dim = c(nsubj, M))
  yhatAligned <- vector("list", nsubj)
  
  ## ------------------------------------------------------------
  ## 2. Helper functions
  ## ------------------------------------------------------------
  
  # Map target grid tj back to the closest points on the workGrid
  get_tJ <- function(tj = workGrid) {
    workGrid[round(tj * (M - 1)) + 1]
  }
  
  # Extract the matrix trajectory Y_j(t) evaluated at tj grid
  getXtJ <- function(Y, j, tj = workGrid) {
    Y[[j]][,, round(tj * (M - 1)) + 1]
  }
  
  # Convert knot parameters into a full warping function on [0,1]
  getSol <- function(res, tGrid) {
    stats::approx(
      x   = seq(0, 1, length.out = (2 + optns$nknots)),
      y   = c(0, sort(res), 1),
      xout = tGrid
    )$y
  }
  
  # Cost function for optimization
  # Integrates Frobenius distance + time regularization penalty
  theCostOptim <- function(x, i, j, lambda, ti) {
    tj <- getSol(x, ti)
    FD <- numeric(length(ti))
    for (k in seq_along(ti)) {
      FD[k] <- FDist2(
        getXtJ(yhat, j, tj[k]),
        getXtJ(yhat, i, ti[k])
      )
    }
    pracma::trapz(ti, FD + lambda * (get_tJ(tj) - ti)^2)
  }
  
  # Optimize the knot positions for g_{j→i}
  getGijOptim <- function(i, j, lambda, minqaAvail) {
    s0 <- seq(0, 1, length.out = (2 + optns$nknots))[2:(1 + optns$nknots)]
    
    if (!minqaAvail) {
      optimRes <- stats::optim(
        par    = s0,
        fn     = theCostOptim,
        method = "L-BFGS-B",
        lower  = rep(1e-6, optns$nknots),
        upper  = rep(1 - 1e-6, optns$nknots),
        i = i, j = j, lambda = lambda, ti = workGrid
      )
    } else {
      optimRes <- minqa::bobyqa(
        par   = s0,
        fn    = theCostOptim,
        lower = rep(1e-6, optns$nknots),
        upper = rep(1 - 1e-6, optns$nknots),
        i = i, j = j, lambda = lambda, ti = workGrid
      )
    }
    getSol(optimRes$par, seq(0, 1, length.out = M))
  }
  
  # Frobenius distance vector over the grid
  getFD <- function(Y, ti, tj, i, j) {
    FD <- numeric(length(ti))
    for (k in seq_along(ti)) {
      FD[k] <- FDist2(
        getXtJ(Y, j, tj[k]),
        getXtJ(Y, i, ti[k])
      )
    }
    FD
  }
  
  ## ------------------------------------------------------------
  ## 3. Check availability of the 'minqa' optimization backend
  ## ------------------------------------------------------------
  if (!"minqa" %in% utils::installed.packages()[, 1] && optns$isPWL) {
    warning("Package 'minqa' not installed. Using 'optim' (L-BFGS-B) instead.")
    minqaAvail <- FALSE
  } else {
    minqaAvail <- TRUE
  }
  ## ------------------------------------------------------------
  ## 4. Perform pairwise warping
  ## ------------------------------------------------------------
  timingWarp <- Sys.time()
  
  for (i in seq_len(nsubj)) {
    if (optns$verbose) {
      cat("Computing pairwise warping for trajectory #", i, "out of", nsubj, "\n")
    }
    
    set.seed(i + optns$seed)
    candidateKcurves <- seq_len(nsubj)[-i]
    
    for (j_idx in seq_len(numOfKcurves)) {
      j <- candidateKcurves[j_idx]
      
      # Optimize g_{j→i}
      gijMat[j_idx, , i] <- getGijOptim(i, j, lambda, minqaAvail)
      
      # Integrated Frobenius distance after alignment
      FD_vec <- getFD(yhat, workGrid, gijMat[j_idx, , i], i, j)
      distMat[i, j_idx] <- sqrt(pracma::trapz(workGrid, FD_vec))
    }
    
    ## ---- Compute inverse warping h^{-1} depending on the averaging rule ----
    if (optns$choice == "weighted") {
      w <- 1 / distMat[i, ]
      w[!is.finite(w)] <- 0
      hInvMat[i, ] <- apply(gijMat[, , i, drop = FALSE], 2, weighted.mean, w)
      
    } else if (optns$choice == "truncated") {
      keep_idx <- (distMat[i, ] <= quantile(distMat[i, ], 0.90))
      hInvMat[i, ] <- apply(gijMat[keep_idx, , i, drop = FALSE], 2, mean)
      
    } else {  
      hInvMat[i, ] <- apply(gijMat[, , i, drop = FALSE], 2, mean)
    }
    
    ## ---- Compute forward warping h by inversion ----
    hMat[i, ] <- stats::approx(
      y = workGrid,
      x = hInvMat[i, ],
      xout = workGrid
    )$y
    
    ## ---- Align trajectory:  X_i(h_i(t)) ----
    yhatAligned[[i]] <- getXtJ(yhat, j = i, tj = hMat[i, ])
  }
  
  timingWarp <- Sys.time() - timingWarp
  
  ## ------------------------------------------------------------
  ## 5. Rescale all warping functions back to the original time range
  ## ------------------------------------------------------------
  for (i in seq_len(nsubj)) {
    for (j_idx in seq_len(numOfKcurves)) {
      gijMat[j_idx, , i] <- gijMat[j_idx, , i] * diff(trange) + trange[1]
    }
  }
  
  h_out    <- hMat    * diff(trange) + trange[1]
  hInv_out <- hInvMat * diff(trange) + trange[1]
  work_out <- workGrid * diff(trange) + trange[1]
  
  ## ------------------------------------------------------------
  ## 6. Return results
  ## ------------------------------------------------------------
  return(list(
    h          = h_out,
    hInv       = hInv_out,
    workGrid   = work_out,
    yhatAligned = yhatAligned,
    optns      = optns,
    costs      = rowMeans(distMat),
    timingWarp = timingWarp,
    g          = gijMat
  ))
}



# -----------------------------------------------------------------------------
# Subject-level warping estimation
# -----------------------------------------------------------------------------

#' @title Estimate subject-level warping functions
#' @description
#' For each component, fit pairwise time-warping on normalized SPD trajectories
#' (d×d×T arrays) using the Frobenius metric
#' to get the subject-level estimate H. The inverse Hinv is obtained approx().
#'
#' @param workGrid Numeric vector of time points (length T).
#' @param Lspd_normd List of length n_comp; each Lspd_normd[[j]] is a list of length n_subj.
#'   Each Lspd_normd[[j]][[i]] is a d×d×T array (SPD trajectory) on workGrid for subject i, component j.
#' @param nknots Integer; number of knots for warping estimation (default: 4).
#' @param lambda Optional numeric vector of length n_comp or a single scalar. If missing,
#'   a default per-component value is computed as lambda_j = 1e-4 × IMSE_j, where
#'   IMSE_j = mean_i ∫ ||Y_{ij}(t) − \bar Y_j(t)||_F^2 dt and \bar Y_j is the Fréchet (Frobenius) mean trajectory.
#'
#' @return A list with:
#' \item{H}{n_subj × T matrix; estimated subject-level warping functions H_i(t).}
#' \item{Hinv}{n_subj × T matrix; numerical inverses H_i^{-1}(t).}
#' \item{workGrid}{The time grid used for evaluation.}
#' \item{Loptns}{List of control options actually used per component.}
#' \item{LtimingWarp}{List of timing information per component.}
#'
#' @importFrom pracma trapz
#' @export
get_subjwf_spd <- function(
    workGrid,
    Lspd_normd,
    nknots = 4,
    lambda
) {
  ## ---- checks ----
  if (!is.numeric(workGrid) || length(workGrid) < 2)
    stop("workGrid must be a numeric vector of length >= 2.")
  if (!is.list(Lspd_normd) || length(Lspd_normd) < 1)
    stop("Lspd_normd must be a non-empty list over components.")
  n_comp <- length(Lspd_normd)
  n_subj_vec <- sapply(Lspd_normd, length)
  if (any(abs(diff(n_subj_vec)) > 0))
    stop("Numbers of subjects differ across components.")
  n_subj <- n_subj_vec[1]
  
  ## ---- lambda handling ----
  missing_lambda <- missing(lambda)
  if (!missing_lambda) {
    if (length(lambda) == 1L) {
      lambda <- rep(lambda, n_comp)
      message("Only one lambda provided — applied to all components.")
    } else if (length(lambda) != n_comp) {
      stop("Length of 'lambda' must be 1 or equal to the number of components.")
    }
  }
  
  ## ---- per-component warping (matrix-valued, Frobenius) ----
  Lres <- lapply(seq_len(n_comp), function(j) {
    Yj <- Lspd_normd[[j]]  # list over subjects; each is d×d×T array
    
    lambda_j <- if (missing_lambda) {
      get_lambda(Yj, workGrid)   # default: IMSE × 1e-4
    } else {
      lambda[j]
    }
    
    FrobPWdense(
      workGrid,
      Yj,
      optns = list(
        lambda = lambda_j,
        nknots = nknots
      )
    )
  })
  names(Lres) <- names(Lspd_normd)
  
  ## ---- average forward warps across components ----
  # In our matrix pipeline, FrobPWdense(...)$hInv is the forward warp H_i(t) (n_subj × T).
  # We average H_i estimates across components to reduce component-specific bias.
  H <- Reduce(`+`, lapply(Lres, function(res) res$hInv)) / n_comp

  ## ---- numeric inverse via approx ----
  # Obtain H_i^{-1} by inverting H_i numerically: swap x and y in approx().
  Hinv <- t(apply(H, 1L, function(h) {
    approx(x = h, y = workGrid, xout = workGrid, rule = 2)$y
  }))
  
  list(
    H           = H,        # n_subj × T forward warps
    Hinv        = Hinv,     # n_subj × T inverse warps
    workGrid    = workGrid,
    Loptns      = lapply(Lres, `[[`, "optns"),
    LtimingWarp = lapply(Lres, `[[`, "timingWarp")
  )
}

# -----------------------------------------------------------------------------
# Grid interpolation helpers (match_times and interp_mat_at)
# -----------------------------------------------------------------------------

#' Find nearest-neighbor indices mapping tout onto tin
#'
#' @description
#' For each element of tout, find the index of the closest element in tin
#' (ties broken toward the left neighbor). Used for "nearest" interpolation
#' of SPD array slices when exact time matching is sufficient.
#'
#' @param tin  Numeric vector; the reference grid (sorted).
#' @param tout Numeric vector; the query time points.
#' @return Integer vector of length length(tout); indices into tin.
match_times <- function (
    tin, # vector of input time points
    tout # vector of desired output time points
) {
  idx <- findInterval(tout, tin)
  # Compare which side is closer, clamping to valid range [1, length(tin)]
  idx <- ifelse(
    idx == 0, 1,
    ifelse(idx == length(tin), length(tin),
           ifelse(abs(tout - tin[idx]) < abs(tout - tin[idx + 1]), idx, idx + 1))
  )
  idx
}

#' Entry-wise linear interpolation of a d x d x T array at a scalar time s
#'
#' @description
#' For each entry (a, b) of the d x d matrix, interpolates the time series
#' arr[a, b, ] (sampled on workGrid) at a single query time s using approx().
#' Then symmetrizes the result to preserve symmetry.
#'
#' @param arr d x d x T array; matrix-valued trajectory sampled on workGrid.
#' @param s   Scalar; query time point.
#' @return A d x d symmetric matrix: interpolated value at time s.
#' @note Uses workGrid from the calling environment (not passed as argument).
interp_mat_at <- function(arr, s) {
  d <- dim(arr)[1]
  out <- matrix(0, d, d)
  for (a in seq_len(d)) for (b in seq_len(d)) {
    out[a, b] <- approx(x = workGrid, y = arr[a, b, ], xout = s, rule = 2)$y
  }
  (out + t(out)) / 2  # symmetrize for numerical stability
}

# -----------------------------------------------------------------------------
# Component alignment and global template estimation
# -----------------------------------------------------------------------------

#' Align component-wise SPD trajectories using subject-level inverse warps
#'
#' @description
#' Given normalized SPD trajectories (component-first: Lspd_normd[[j]][[i]] is d×d×T on workGrid)
#' and subject-level inverse warps Hinv (n×T), align time per subject to obtain (primitive) component tempo trajectories. 
#'
#' @param workGrid Numeric vector of length T (evaluation time grid).
#' @param Lspd_normd List of length p (components). Each Lspd_normd[[j]] is a list of length n (subjects),
#'   and each Lspd_normd[[j]][[i]] is a d×d×T array (already normalized).
#' @param Hinv Numeric matrix n×T: subject-level inverse warps H_i^{-1}(t) evaluated on workGrid.
#' @param method "nearest" (use match_times) or "interp" (entry-wise linear interpolation). Default "nearest".
#' @param return_mean Logical; if TRUE, also return component-wise Fréchet mean trajectories (Etahat).
#'
#' @return A list with
#'   - aligned: list length p; each [[j]] is a list length n of d×d×T aligned arrays
#'   - Etahat:  (optional) list length p; each [[j]] is d×d×T (component-wise mean)
#' @examples
#' # res <- align_comp_spd(workGrid, Xhat, Hinv = res_H$Hinv, method = "nearest", return_mean = TRUE)
#' @export
align_comp_spd <- function(workGrid, Lspd_normd, Hinv,
                           method = c("nearest", "interp"),
                           return_mean = TRUE) {
  method <- match.arg(method)
  
  # ---- basic checks ----
  if (!is.numeric(workGrid) || length(workGrid) < 2)
    stop("workGrid must be a numeric vector (length >= 2).")
  if (!is.list(Lspd_normd) || length(Lspd_normd) < 1)
    stop("Lspd_normd must be a non-empty list over components [[j]][[i]].")
  p <- length(Lspd_normd)
  n_vec <- sapply(Lspd_normd, length)
  if (any(abs(diff(n_vec)) > 0)) stop("Number of subjects differs across components.")
  n <- n_vec[1]
  
  ref <- Lspd_normd[[1]][[1]]
  if (!is.array(ref) || length(dim(ref)) != 3L) stop("Each trajectory must be a d×d×T array.")
  d <- dim(ref)[1]; Tlen <- dim(ref)[3]
  if (dim(ref)[1] != dim(ref)[2]) stop("Arrays must be square (d×d×T).")
  if (length(workGrid) != Tlen) stop("length(workGrid) must equal T (third dim of arrays).")
  if (!is.matrix(Hinv) || nrow(Hinv) != n || ncol(Hinv) != Tlen)
    stop("Hinv must be an n×T matrix matching subjects and workGrid length.")
  
  
  # ---- alignment per component ----
  aligned <- vector("list", p)
  for (j in seq_len(p)) {
    aligned[[j]] <- vector("list", n)
    for (i in seq_len(n)) {
      arr <- array(0, dim = c(d, d, Tlen))
      if (method == "nearest") {
        idx <- match_times(tin = workGrid, tout = Hinv[i, ])
        for (k in seq_len(Tlen)) {
          M <- Lspd_normd[[j]][[i]][, , idx[k]]
          arr[, , k] <- (M + t(M)) / 2
        }
      } else { # "interp"
        for (k in seq_len(Tlen)) {
          s <- Hinv[i, k]
          arr[, , k] <- interp_mat_at(Lspd_normd[[j]][[i]], s)
        }
      }
      aligned[[j]][[i]] <- arr
    }
  }
  
  # ---- optional: component-wise Fréchet mean (Frobenius) at each time ----
  out <- list(aligned = aligned)
  if (isTRUE(return_mean)) {
    Etahat <- lapply(aligned, Fmean_list)  # uses frechet::CovFMean internally
    out$Etahat <- Etahat
  }
  out
}

#' Preliminary global deformation for SPD trajectories (Frobenius metric)
#'
#' @description
#' Pools all component-subject SPD trajectories, applies global pairwise warping
#' via `FrobPWdense`, and restructures outputs back into component-wise lists.
#'
#' @param workGrid Numeric vector of length T (evaluation time grid).
#' @param Lspd_normd List of length p (components). Each Lspd_normd[[j]] is a
#'   list of length n (subjects), each element a d×d×T SPD array.
#' @param nknots Integer, number of warping knots (default 4).
#' @param lambda Optional scalar. If missing, determined automatically by `get_lambda()`.
#'
#' @return A list containing:
#' \itemize{
#'   \item{H}{ list length p; each [[j]] is n×T matrix of global deformation functions.}
#'   \item{Hinv}{ list length p; each [[j]] is n×T matrix of inverse deformation functions.}
#'   \item{Lspd_normd_aligned}{ list length p; each [[j]] is a list of d×d×T aligned SPD trajectories.}
#'   \item{workGrid}{ The time grid.}
#'   \item{optns}{ The options passed to `FrobPWdense`.}
#'   \item{timingWarp}{ Runtime information.}
#' }
#' @export
align_glob_spd <- function(workGrid, Lspd_normd, nknots = 4, lambda) {
  
  # --- basic checks ---
  if (!is.numeric(workGrid)) stop("workGrid must be numeric.")
  p <- length(Lspd_normd)
  n_vec <- sapply(Lspd_normd, length)
  if (any(n_vec != n_vec[1])) stop("Different numbers of subjects across components.")
  n <- n_vec[1]
  
  ref <- Lspd_normd[[1]][[1]]
  if (!is.array(ref) || length(dim(ref)) != 3L)
    stop("Each element of Lspd_normd[[j]][[i]] must be a d×d×T array.")
  d <- dim(ref)[1]; Tlen <- dim(ref)[3]
  if (d != dim(ref)[2]) stop("SPD arrays must be square.")
  if (length(workGrid) != Tlen) stop("workGrid length must equal T.")
  
  # --- pool all subjects across components ---
  # Concatenates all n subjects from all p components into one list of length p*n.
  # FrobPWdense treats each element as an independent curve for global pairwise warping.
  yhat <- do.call(c, Lspd_normd)   # list length = p * n

  # --- set lambda if missing ---
  if (missing(lambda)) {
    lambda <- get_lambda(yhat, workGrid, scale = 1e-4)
  } else if (length(lambda) != 1) {
    warning("lambda has length > 1; using the first element.")
    lambda <- lambda[1]
  }
  
  # --- call pairwise Frobenius warping ---
  res <- FrobPWdense(
    tVec = workGrid,
    yhat = yhat,
    optns = list(
      lambda = lambda,
      nknots = nknots
    )
  )
  
  # --- index map for component blocks ---
  # The pooled list has subjects ordered as: comp1 subjects, comp2 subjects, ...
  # idx_by_comp[[j]] gives the row indices in the (p*n)-row warp matrices for component j.
  idx_by_comp <- lapply(seq_len(p), function(j) {
    base <- (j - 1L) * n
    base + seq_len(n)
  })

  # According to FrobPWdense return semantics:
  # hInv = forward warp H (what we call "H"); h = inverse warp H^{-1} (what we call "Hinv").
  H    <- lapply(seq_len(p), function(j) res$hInv[idx_by_comp[[j]], , drop=FALSE])
  Hinv <- lapply(seq_len(p), function(j) res$h   [idx_by_comp[[j]], , drop=FALSE])
  
  # --- aligned trajectories ---
  Lspd_normd_aligned <- lapply(seq_len(p), function(j) {
    out_j <- res$yhatAligned[idx_by_comp[[j]]]
    names(out_j) <- names(Lspd_normd[[j]])
    out_j
  })
  names(Lspd_normd_aligned) <- names(H) <- names(Hinv) <- names(Lspd_normd)
  
  # --- return result ---
  list(
    H = H,
    Hinv = Hinv,
    Lspd_normd_aligned = Lspd_normd_aligned,
    workrid = workGrid,
    optns = res$optns,
    timingWarp = res$timingWarp
  )
}

#' Compute the latent template (SPD) 
#'
#' @description
#' Pools all globally aligned SPD trajectories across components and subjects,
#' then computes the Fréchet mean (Frobenius metric).
#'
#' @param workGrid Numeric vector of length T (evaluation grid).
#' @param Lspd_normd_aligned List of length p (components). Each [[j]] is a list
#'   of length n (subjects), and each [[j]][[i]] is a d x d x T array (globally aligned).
#'
#' @return A list with:
#' \itemize{
#'   \item{Tau}{ d x d x T array; latent template trajectory (Fréchet mean).}
#'   \item{workGrid}{ the input time grid (returned for convenience).}
#' }
#' @examples
#' # Tau_res <- get_template_spd(workGrid, Lspd_normd_aligned)
get_template_spd <- function(workGrid, Lspd_normd_aligned) {
  # ---- basic checks ----
  if (!is.numeric(workGrid) || length(workGrid) < 2)
    stop("workGrid must be a numeric vector (length >= 2).")
  if (!is.list(Lspd_normd_aligned) || length(Lspd_normd_aligned) < 1)
    stop("Lspd_normd_aligned must be a non-empty list over components [[j]][[i]].")
  
  p <- length(Lspd_normd_aligned)
  n_vec <- sapply(Lspd_normd_aligned, length)
  if (any(n_vec != n_vec[1])) stop("Numbers of subjects differ across components.")
  n <- n_vec[1]
  
  ref <- Lspd_normd_aligned[[1]][[1]]
  if (!is.array(ref) || length(dim(ref)) != 3L)
    stop("Each trajectory must be a d x d x T array.")
  d <- dim(ref)[1]; Tlen <- dim(ref)[3]
  if (d != dim(ref)[2]) stop("Arrays must be square (d x d x T).")
  if (length(workGrid) != Tlen) stop("length(workGrid) must equal T.")
  
  # ---- pool all (j,i) and compute time-wise Fréchet mean (Frobenius) ----
  pooled <- do.call(c, Lspd_normd_aligned)     # length = p * n; elements are d x d x T
  tau <- Fmean_list(pooled)                    # returns d x d x T
  
  list(
    tau = tau,
    workGrid = workGrid
  )
}





# -----------------------------------------------------------------------------
# CM template estimation (get_temp_spd: pairwise warping on representative Z_i)
# -----------------------------------------------------------------------------

#' @title Randomly sample representative processes \eqn{Z^*_{i}}
#' @description
#' Given SPD trajectories 
#' (\code{Lspd[[j]][[i]]} is a d×d×T array for component j and subject i),
#' randomly pick one component per subject to form a list
#' (\code{spd_z[[i]]}).
#'
#' @param Lspd A list of length p (components). Each \code{Lspd[[j]]} is a list
#' of length n (subjects), and each \code{Lspd[[j]][[i]]} is a d×d×T SPD array.
#'
#' @return A list of length n where each \code{spd_z[[i]]} is the selected SPD
#' trajectory for subject i.
#'
#' @examples
#' # Example: component-first SPD data (Lspd[[j]][[i]])
#' # spd_z <- get_repr_spd(spd_list)
#'
#' @export
get_repr_spd <- function(Lspd, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  # Basic validation
  if (!is.list(Lspd) || length(Lspd) < 1)
    stop("Lspd must be a non-empty list of components (Lspd[[j]][[i]] is d×d×T).")
  
  p <- length(Lspd)          # number of components
  n <- length(Lspd[[1]])     # number of subjects
  
  # Randomly assign a component index to each subject
  idx <- sample(seq_len(p), size = n, replace = TRUE)
  
  # Construct subject-first list spd_z[[i]]
  spd_z <- vector("list", n)
  for (i in seq_len(n)) {
    spd_z[[i]] <- Lspd[[ idx[i] ]][[ i ]]
  }
  
  spd_z
}


#' @title Estimate latent template for SPD trajectories.
#' @description
#' Given normalized SPD trajectories
#' (Lspd_normd[[j]][[i]] is d×d×T on workGrid), this function estimates latent template for SPD trajectories.
#'
#' @param workGrid Numeric vector of time points (length T).
#' @param Lspd_normd List of length p (components). Each Lspd_normd[[j]] is a list of length n (subjects),
#'   and each Lspd_normd[[j]][[i]] is a d×d×T SPD/PSD array (already normalized) on workGrid.
#' @param nknots Integer; number of knots for warping estimation (default: 4).
#' @param lambda Optional scalar regularization. If missing, set to 1e-4 × IMSE(Z, workGrid).
#' @param method Alignment for building the template: "nearest" (index matching) or "interp" (entry-wise linear interpolation).
#'
#' @return A list with:
#' \item{tau}{d×d×T array: the estimated latent template trajectory.}
#' \item{W}{n×T matrix: forward global warps (res$hInv).}
#' \item{Winv}{n×T matrix: inverse global warps (res$h).}
#' \item{lambda}{Scalar lambda actually used.}
#' \item{Z}{Subject-first representative processes used (list length n).}
#' \item{aligned_Z}{Subject-first aligned Z (list length n).}
#'
#' @examples
#' # Tau <- get_template_spd(workGrid, Lspd_normd = Xhat, nknots = 4)
#'
#' @export
get_temp_spd <- function(workGrid, Lspd_normd, nknots = 4, lambda, method = c("nearest","interp"),seed=NULL) {
  method <- match.arg(method)
  
  ## ---- basic checks ----
  if (!is.numeric(workGrid) || length(workGrid) < 2)
    stop("workGrid must be a numeric vector of length >= 2.")
  if (!is.list(Lspd_normd) || length(Lspd_normd) < 1)
    stop("Lspd_normd must be a non-empty component-first list [[j]][[i]].")
  p <- length(Lspd_normd)
  n_vec <- sapply(Lspd_normd, length)
  if (any(abs(diff(n_vec)) > 0)) stop("Numbers of subjects differ across components.")
  n <- n_vec[1]
  
  ref <- Lspd_normd[[1]][[1]]
  if (!is.array(ref) || length(dim(ref)) != 3L)
    stop("Each Lspd_normd[[j]][[i]] must be a 3D array (d×d×T).")
  d <- dim(ref)[1]; Tlen <- dim(ref)[3]
  if (dim(ref)[1] != dim(ref)[2]) stop("Arrays must be square (d×d×T).")
  if (length(workGrid) != Tlen) stop("length(workGrid) must equal T (third dim).")
  
  ## ---- helper: Frobenius distance (uses SMFilter::FDist2 if available) ----
  FDist2_local <- if (requireNamespace("SMFilter", quietly = TRUE) &&
                      exists("FDist2", where = asNamespace("SMFilter"), inherits = FALSE)) {
    SMFilter::FDist2
  } else {
    function(A, B) sum((A - B)^2)
  }
  
  
  ## ---- step 1: representative processes Z_i* (subject-first) ----
  Z <- get_repr_spd(Lspd_normd, seed=seed)  # returns list length n with d×d×T arrays
  
  ## ---- step 2: lambda (if missing) ----
  missing_lambda <- missing(lambda)
  if (missing_lambda) {
    lambda <- get_lambda(Z, workGrid)
  } else {
    if (!is.numeric(lambda) || length(lambda) != 1L)
      stop("lambda must be a scalar.")
  }
  
  ## ---- step 3: global warps W via FrobPWdense on subject-first Z ----
  # FrobPWdense(workGrid, Yj) expects a list over subjects; use Z directly.
  res <- FrobPWdense(
    workGrid,
    Z,
    optns = list(lambda = lambda, nknots = nknots)
  )
  W    <- res$hInv   # forward
  Winv <- res$h      # inverse
  
  ## ---- step 4: align Z by Winv and form template tau ----
  # alignment: nearest index or interpolation
  
  aligned_Z <- vector("list", n)
  if (method == "nearest") {
    for (i in seq_len(n)) {
      out <- array(0, dim = c(d, d, Tlen))
      idx <- match_times(workGrid, Winv[i, ])
      for (k in seq_len(Tlen)) {
        M <- Z[[i]][, , idx[k]]
        out[, , k] <- (M + t(M)) / 2
      }
      aligned_Z[[i]] <- out
    }
  } else { # "interp"
    for (i in seq_len(n)) {
      out <- array(0, dim = c(d, d, Tlen))
      for (k in seq_len(Tlen)) out[, , k] <- interp_mat_at(Z[[i]], Winv[i, k])
      aligned_Z[[i]] <- out
    }
  }
  
  tau <- Fmean_list(aligned_Z)
  
  list(
    tau    = tau,
    W         = W,
    Winv      = Winv,
    lambda    = lambda,
    Z         = Z,
    aligned_Z = aligned_Z
  )
}


# -----------------------------------------------------------------------------
# Component-level warping functions (Psi) and tempo trajectories
# -----------------------------------------------------------------------------

#' @title Estimate component-level warping functions.
#' @description
#' Get the component-level warping function Psi using the latent template
#' \code{tau} (d×d×T) and the component-specific tempo trajectory \code{Eta_j} (d×d×T). 
#' The inverse of warping function is obtained approx().
#' @param workGrid Numeric vector of time points (length T).
#' @param tau   d×d×T array: latent template trajectory on \code{workGrid}.
#' @param Leta     List of length p (components). Each \code{Leta[[j]]} is a d×d×T SPD/PSD array
#'                 giving the component-specific aligned mean trajectory (tempo) on \code{workGrid}.
#' @param nknots   Integer; number of knots for warping estimation (default: 4).
#' @param lambda   Optional scalar regularization. If missing, set to
#'                 \eqn{10^{-4} \times \frac{1}{p}\sum_j \int \| \Eta_j(t) - \tau(t) \|_F^2\,dt}.
#'
#' @return A list with:
#' \item{H}{list of length p; each element is a numeric vector (length T) holding the estimates of component-level warping functions on time points in \code{workGrid}.}
#' \item{Hinv}{list of length p; each element is a numeric vector (length T) holding the inverse of component-level warping functions on time points in \code{workGrid}.}
#' \item{workGrid}{The time grid used for evaluation.}
#' \item{Loptns}{List of control options used for each component.}
#' \item{LtimingWarp}{List of timing information per component.}
#'
#' @export
get_compwf_spd <- function(workGrid, tau, Leta, nknots = 4, lambda) {
  ## ---- checks ----
  if (!is.numeric(workGrid) || length(workGrid) < 2)
    stop("workGrid must be a numeric vector of length >= 2.")
  if (!is.array(tau) || length(dim(tau)) != 3L)
    stop("tau must be a 3D array (d×d×T).")
  d <- dim(tau)[1]; Tlen <- dim(tau)[3]
  if (dim(tau)[1] != dim(tau)[2]) stop("tau must be square (d×d×T).")
  if (length(workGrid) != Tlen) stop("length(workGrid) must equal T (third dim of tau).")
  
  if (!is.list(Leta) || length(Leta) < 1)
    stop("Leta must be a non-empty list over components; each Leta[[j]] is d×d×T.")
  p <- length(Leta)
  for (j in seq_len(p)) {
    A <- Leta[[j]]
    if (!is.array(A) || length(dim(A)) != 3L ||
        any(dim(A)[1:2] != c(d, d)) || dim(A)[3] != Tlen) {
      stop("All Leta[[j]] must be d×d×T with the same dims as tau.")
    }
  }
  
  ## ---- default lambda if missing ----
  missing_lambda <- missing(lambda)
  if (missing_lambda) {
    # Prefer SMFilter::FDist2 if available
    FDist2_local <- if (requireNamespace("SMFilter", quietly = TRUE) &&
                        exists("FDist2", where = asNamespace("SMFilter"), inherits = FALSE)) {
      SMFilter::FDist2
    } else {
      function(A, B) sum((A - B)^2)
    }
    # Average IMSE over components
    imse_j <- vapply(seq_len(p), function(j) {
      difft <- vapply(seq_len(Tlen), function(k) {
        FDist2_local(Leta[[j]][, , k], tau[, , k])
      }, numeric(1))
      pracma::trapz(workGrid, difft)
    }, numeric(1))
    lambda <- 1e-4 * mean(imse_j)
  } else if (!is.numeric(lambda) || length(lambda) != 1L) {
    stop("lambda must be a numeric scalar.")
  }
  
  ## ---- per-component pairwise warping against tau ----
  Lres <- lapply(seq_len(p), function(j) {
    # Two-curve warping: subject list = {tau, Eta_j}
    res2 <- FrobPWdense(
      workGrid,
      list(tau, Leta[[j]]),
      optns = list(lambda = lambda, nknots = nknots)
    )
    # Follow the distributional convention:
    # take the warp associated with the 2nd curve (component j)
    list(
      hInv       = res2$hInv[2, ],  # forward warp H_j(t)
      h         = res2$h[2, ],     # inverse warp H_j^{-1}(t)
      optns      = res2$optns,
      timingWarp = res2$timingWarp
    )
  })
  
  list(
    H           = lapply(Lres, `[[`, "hInv"),
    Hinv        = lapply(Lres, `[[`, "h"),
    workGrid    = workGrid,
    Loptns      = lapply(Lres, `[[`, "optns"),
    LtimingWarp = lapply(Lres, `[[`, "timingWarp")
  )
}

#' @title Component tempos by warping the template trajectory. 
#' @description
#' Warp the latent template trajectory tau (d×d×T) by each component-level
#' warping function Psi_j (length T on workGrid) to obtain component-specific
#' tempo trajectories.
#'
#' @param workGrid Numeric vector of time points (length T).
#' @param tau   3D array d×d×T: latent template trajectory on workGrid.
#' @param Psi      List of length p; each Psi[[j]] is a numeric vector (length T)
#'                 giving the component-level forward warp evaluated on workGrid.
#' @param method   "nearest" (index matching via match_times) or
#'                 "interp" (entry-wise linear interpolation). Default "nearest".
#'
#' @return A list of length p; each element is a d×d×T array, i.e., the resolved
#'         component tempo aligned by Psi_j.
#'
#' @examples
#' # Leta <- get_comptempo_spd(workGrid, tau, Psi, method = "nearest")
#'
#' @export
get_comptempo_spd <- function(workGrid, tau, Psi, method = c("nearest", "interp")) {
  method <- match.arg(method)
  
  # ---- basic checks ----
  if (!is.numeric(workGrid) || length(workGrid) < 2)
    stop("workGrid must be a numeric vector (length >= 2).")
  if (!is.array(tau) || length(dim(tau)) != 3L)
    stop("tau must be a 3D array (d×d×T).")
  d <- dim(tau)[1]; Tlen <- dim(tau)[3]
  if (dim(tau)[1] != dim(tau)[2]) stop("tau must be square (d×d×T).")
  if (length(workGrid) != Tlen) stop("length(workGrid) must equal T (third dim of tau).")
  if (!is.list(Psi) || length(Psi) < 1) stop("Psi must be a non-empty list of component warps.")
  p <- length(Psi)
  for (j in seq_len(p)) {
    if (!is.numeric(Psi[[j]]) || length(Psi[[j]]) != Tlen)
      stop("Each Psi[[j]] must be a numeric vector of length T (workGrid).")
  }
  
  # ---- build component tempos ----
  Leta <- vector("list", p)
  
  if (method == "nearest") {
    # require match_times(tin, tout) to be defined in the environment
    for (j in seq_len(p)) {
      idx <- match_times(tin = workGrid, tout = Psi[[j]])
      arr <- array(0, dim = c(d, d, Tlen))
      for (k in seq_len(Tlen)) {
        M <- tau[, , idx[k]]
        arr[, , k] <- (M + t(M)) / 2
      }
      Leta[[j]] <- arr
    }
  } else { # "interp"
    for (j in seq_len(p)) {
      arr <- array(0, dim = c(d, d, Tlen))
      for (k in seq_len(Tlen)) {
        s <- Psi[[j]][k]
        arr[, , k] <- interp_mat_at(tau, s)
      }
      Leta[[j]] <- arr
    }
  }
  
  Leta
}

# -----------------------------------------------------------------------------
# Global deformation composition and alignment
# -----------------------------------------------------------------------------

# ---- helper: compose two warps on a common grid --------------------------------
# Given two monotone warps f1 and f2 sampled on the same grid 'grid',
# return the composition f1(f2(x)) evaluated at x = 'grid'; numeric vector of length T.
composite_two <- function(f1_vec, f2_vec, grid) {
  # f1, f2 are numeric vectors of length T on 'grid'
  # Step 1) evaluate u(x) = f2(x) on grid (already sample-aligned)
  u <- f2_vec
  # Step 2) evaluate f1(u) by interpolation back on 'grid'
  as.numeric(approx(x = grid, y = f1_vec, xout = u, rule = 2)$y)
}

# ---- main: final global deformation functions (SPD version) --------------------
#' Final global deformation functions (SPD; Frobenius metric)
#'
#' @description
#' Given subject-level warps (from `get_subjwf_spd`) and component-level warps
#' (from `get_compwf_spd`), compose them to obtain final global deformation
#' functions per component and subject on a target `workGrid`.
#'
#' @param res_H   Output of `get_subjwf_spd`, with elements:
#'                H (n x T_H), Hinv (n x T_H), workGrid.
#' @param res_Psi Output of `get_compwf_spd`, with elements:
#'                H (list length p with n x T_P), Hinv (list length p with n x T_P), workGrid.
#' @param workGrid Numeric vector: target grid on which to report final warps.
#'
#' @return A list with
#' \itemize{
#'   \item{H}{ list length p; each [[j]] is n x T matrix of forward global deformations \(G_{i,j}(t)\).}
#'   \item{Hinv}{ list length p; each [[j]] is n x T matrix of inverse deformations \(G^{-1}_{i,j}(t)\).}
#'   \item{workGrid}{ the target grid used for evaluation.}
#' }
#' @export
get_globwf_spd <- function(res_H, res_Psi, workGrid) {
  # ---- basic checks ----
  stopifnot(is.numeric(workGrid), length(workGrid) >= 2)
  
  # subject-level warps (matrices n x T_H)
  H_subj    <- res_H$H
  Hinv_subj <- res_H$Hinv
  grid_H    <- res_H$workGrid
  if (!is.matrix(H_subj) || !is.matrix(Hinv_subj))
    stop("res_H$H and res_H$Hinv must be n x T matrices.")
  if (nrow(H_subj) != nrow(Hinv_subj))
    stop("Row mismatch between res_H$H and res_H$Hinv.")
  n <- nrow(H_subj)
  
  # component-level warps (vectors length T_P per component)
  Psi_list    <- res_Psi$H
  PsiInv_list <- res_Psi$Hinv
  grid_Psi    <- res_Psi$workGrid
  if (!is.list(Psi_list) || !is.list(PsiInv_list) || length(Psi_list) < 1L)
    stop("res_Psi$H and res_Psi$Hinv must be non-empty lists of vectors.")
  p <- length(Psi_list)
  if (length(PsiInv_list) != p)
    stop("res_Psi$H and res_Psi$Hinv must have the same length.")
  
  # ---- helpers: resample to target grid and compose f1∘f2 on 'grid' ----
  to_grid_mat <- function(mat, x_old, x_new) {
    if (isTRUE(all.equal(x_old, x_new))) return(mat)
    t(apply(mat, 1, function(v) approx(x = x_old, y = v, xout = x_new, rule = 2)$y))
  }
  to_grid_vec <- function(vec, x_old, x_new) {
    if (isTRUE(all.equal(x_old, x_new))) return(vec)
    approx(x = x_old, y = vec, xout = x_new, rule = 2)$y
  }
  # f1_vec(f2_vec(t_k)): evaluate f1 at inputs f2(t_k)
  compose_on_grid <- function(f1_vec, f2_vec, grid) {
    # both f1_vec and f2_vec are sampled on 'grid'
    approx(x = grid, y = f1_vec, xout = f2_vec, rule = 2)$y
  }
  
  # ---- resample everything onto 'workGrid' ----
  H_subj_g    <- to_grid_mat(H_subj,    grid_H,  workGrid)
  Hinv_subj_g <- to_grid_mat(Hinv_subj, grid_H,  workGrid)
  Psi_g       <- lapply(Psi_list,    to_grid_vec, x_old = grid_Psi, x_new = workGrid)
  PsiInv_g    <- lapply(PsiInv_list, to_grid_vec, x_old = grid_Psi, x_new = workGrid)
  
  # ---- compose per component j and subject i ----
  H_out    <- vector("list", p)   # forward:  G_{i,j} = Psi_j ∘ H_i
  Hinv_out <- vector("list", p)   # inverse:  G^{-1}_{i,j} = H^{-1}_i ∘ Psi^{-1}_j
  
  for (j in seq_len(p)) {
    if (!is.numeric(Psi_g[[j]]) || length(Psi_g[[j]]) != length(workGrid))
      stop("Each res_Psi$H[[j]] must be a numeric vector of length T after resampling.")
    if (!is.numeric(PsiInv_g[[j]]) || length(PsiInv_g[[j]]) != length(workGrid))
      stop("Each res_Psi$Hinv[[j]] must be a numeric vector of length T after resampling.")
    
    Gj_mat    <- matrix(NA_real_, nrow = n, ncol = length(workGrid))
    GjInv_mat <- matrix(NA_real_, nrow = n, ncol = length(workGrid))
    
    for (i in seq_len(n)) {
      # forward composition
      Gj_mat[i, ]    <- compose_on_grid(f1_vec = Psi_g[[j]],
                                        f2_vec = H_subj_g[i, ],
                                        grid   = workGrid)
      # inverse composition
      GjInv_mat[i, ] <- compose_on_grid(f1_vec = Hinv_subj_g[i, ],
                                        f2_vec = PsiInv_g[[j]],
                                        grid   = workGrid)
    }
    H_out[[j]]    <- Gj_mat
    Hinv_out[[j]] <- GjInv_mat
  }
  
  # optional: carry names over if present
  names(H_out)    <- names(Psi_list)
  names(Hinv_out) <- names(PsiInv_list)
  
  list(
    H        = H_out,     # list length p; each n x T forward global warps
    Hinv     = Hinv_out,  # list length p; each n x T inverse global warps
    workGrid = workGrid
  )
}



#' Align Xhat using final inverse global warps G^{-1}_{i,j} (SPD, component-first)
#'
#' @param workGrid Numeric vector of length T.
#' @param Xhat Component-first SPD trajectories; Xhat[[j]][[i]] is d×d×T.
#' @param res_G Output of `get_globwf_spd()`, where
#'        res_G$Hinv[[j]] is an n×T matrix giving G^{-1}_{i,j}(t).
#' @param method "nearest" (default) or "interp".
#'
#' @return Same shape as Xhat: out[[j]][[i]] is aligned d×d×T.
#' @export
align_Xhat_Ginv_spd <- function(workGrid, Xhat, res_G, 
                                method = c("nearest","interp")) {
  method <- match.arg(method)
  
  # ---- basic checks ----
  if (!is.numeric(workGrid) || length(workGrid) < 2)
    stop("workGrid must be numeric.")
  if (!is.list(Xhat) || length(Xhat) < 1)
    stop("Xhat must be a non-empty component-first list.")
  
  p <- length(Xhat)
  n <- length(Xhat[[1]])
  Tlen <- length(workGrid)
  
  if (!is.list(res_G$Hinv) || length(res_G$Hinv) != p)
    stop("res_G$Hinv must be a list of length p.")
  if (any(sapply(res_G$Hinv, nrow) != n))
    stop("Each res_G$Hinv[[j]] must have n rows.")
  if (any(sapply(res_G$Hinv, ncol) != Tlen))
    stop("Each res_G$Hinv[[j]] must have length T columns equal to workGrid length.")
  
  ref <- Xhat[[1]][[1]]
  if (!is.array(ref) || length(dim(ref)) != 3L)
    stop("Each Xhat[[j]][[i]] must be d×d×T.")
  d <- dim(ref)[1]
  if (dim(ref)[1] != dim(ref)[2] || dim(ref)[3] != Tlen)
    stop("Mismatch: Xhat must be d×d×T with T = length(workGrid).")
  
  # ---- align ----
  Xout <- vector("list", p)
  for (j in seq_len(p)) {
    Hinv_j <- res_G$Hinv[[j]]   # n×T
    Xout[[j]] <- vector("list", n)
    for (i in seq_len(n)) {
      
      # final inverse warp trajectory s(t) = G^{-1}_{i,j}(t)
      s_vec <- Hinv_j[i, ]
      
      out <- array(0, dim = c(d, d, Tlen))
      if (method == "nearest") {
        idx <- match_times(tin = workGrid, tout = s_vec)
        for (k in seq_len(Tlen)) {
          M <- Xhat[[j]][[i]][,, idx[k]]
          out[,,k] <- (M + t(M)) / 2  # symmetrize only
        }
      } else {  # "interp"
        for (k in seq_len(Tlen)) {
          out[,,k] <- interp_mat_at(Xhat[[j]][[i]], s_vec[k])
        }
      }
      
      Xout[[j]][[i]] <- out
    }
  }
  
  Xout
}

# -----------------------------------------------------------------------------
# Top-level GLDM fitting functions
# -----------------------------------------------------------------------------

#' Fit GLDM (Geodesic latent deformation model) for SPD Matrix-Valued Trajectories
#'
#' @description
#' Perform GLDM to SPD-valued trajectories (matrix processes).
#' @param obsGrid   numeric vector: time grid at which X processes are observed (length T1).
#' @param workGrid  numeric vector: time grid where estimated warping functions are evaluated (length T2).
#' @param Y      list of observed SPD-valued trajectories with potential contamination per component per subject.
#'                  Y[[j]][[i]] is a 3D array (d × d × T1) for subject i, component j.
#' @param nknots    integer; number of spline knots used for warping estimation.
#' @param lambdaH  vector of length p; regularization parameters for subject-level warp estimation.
#' @param lambdaW  scalar; regularization parameter for template estimation.
#' @param lambdaPsi scalar; regularization parameter for component-level warp estimation.
#' @param method    character; interpolation or alignment method (default "nearest").
#' @param sigma_pert  scalar indicating whether the observed processes Y is contaminated or not.
#'
#' @return A list containing all intermediate and final alignment results:
#' \itemize{
#'   \item \code{Ahat}: n x p matrix of amplitude factors.
#'   \item \code{norm_data}: amplitude-normalized SPD trajectories (X*).
#'   \item \code{Xhat}: smoothed SPD trajectories (from LocCovReg when doLFR=TRUE).
#'   \item \code{res_H}: subject-level warping functions (H_i).
#'   \item \code{res_eta_pre}: component-level aligned trajectories (eta_j, primitive).
#'   \item \code{res_G_pre}: primitive global deformation functions.
#'   \item \code{res_G_ref}: refined global deformation functions (GLDM).
#'   \item \code{res_G_old}: global deformation functions (CM).
#'   \item \code{res_tau}: global latent template (GLDM).
#'   \item \code{res_tau_old}: global latent template (CM).
#'   \item \code{res_Psi}: component-level warping functions (GLDM).
#'   \item \code{res_Psi_old}: component-level warping functions (CM).
#'   \item \code{res_eta}: resolved component tempos.
#'   \item \code{Xhat_ginv_pre}: globally aligned trajectories via primitive G.
#'   \item \code{Xhat_ginv_ref}: globally aligned trajectories via refined G (GLDM).
#'   \item \code{Xhat_ginv_old}: globally aligned trajectories via CM G.
#'   \item \code{workGrid}: the time grid used.
#' }
#'
#' @details
#' Each helper function (`get_amp`, `normalize`, `get_subjwf_spd`, etc.)
#' must already be defined. They should handle SPD-valued data using
#' Log-Euclidean or affine-invariant distances.
#'
#' @export
gldm_spd_once <- function(obsGrid, Y, nknots = 4, lambdaH,lambdaW,lambdaG, lambdaPsi,method = "nearest",doLFR=TRUE, seed=NULL) {

  if (!is.list(Y) || length(Y) == 0L)
    stop("Y must be a non-empty list: Y[[j]][[i]] is a d×d×T1 array.")
  p <- length(Y)
  if (!is.list(Y[[1]]) || length(Y[[1]]) == 0L)
    stop("Each Y[[j]] must be a non-empty list of subjects.")

  n <- length(Y[[1]])
  d <- dim(Y[[1]][[1]])[1]

  # ---- (0) Xhat estimation ----
  m <- length(obsGrid)
  x_mat  <- matrix(obsGrid, nrow = length(obsGrid))
  Xhat_array <- vector("list", p)
  # Xhat (estimated from noisy observations Y)
  for (j in 1:p) {
    Xhat_array[[j]] <- vector("list", n)
    for (i in 1:n) {
      Xhat_array[[j]][[i]] <- array(0, dim = c(d, d, m))
      if (!doLFR) {
        Xhat_array[[j]][[i]] <- Y[[j]][[i]]
      } else {
        y <- LocCovReg(x = x_mat, M = Y[[j]][[i]],
                       xout = x_mat,
                       optns = list(metric = "frobenius", kernel = "epan"))
        for (kk in 1:m) Xhat_array[[j]][[i]][, , kk] <- y$Mout[[kk]]
      }
    }
  }
  
  
  # ---- (1) Amplitude and normalization ----
  
  Ahat <- get_amp(Xhat_array)
  
  norm_data <- normalize(Xhat_array, Ahat)
  
  # ---- (2) Subject-level warping ----
  res_H <- get_subjwf_spd(
    workGrid = obsGrid,
    Lspd_normd = norm_data,
    nknots = nknots,
    lambda = lambdaH
  )
  
  # ---- (3) Component-level alignment ----
  res_eta_pre <- align_comp_spd(
    workGrid = obsGrid,
    Lspd_normd = norm_data,
    Hinv = res_H$Hinv,
    method = method
  )
  
  res_G_pre <- align_glob_spd(
    workGrid = obsGrid,
    Lspd_normd = norm_data,
    nknots = nknots,
    lambda = lambdaG
  )

  # ---- (4) Global latent template ----
  res_tau <- get_template_spd(
    workGrid = obsGrid,
    Lspd_normd = res_G_pre$Lspd_normd_aligned
  )
  
  # ---- (5) Component-level warping functions ----
  res_Psi <- get_compwf_spd(
    workGrid = obsGrid,
    tau = res_tau$tau,
    Leta = res_eta_pre$Etahat,
    nknots = nknots,
    lambda = lambdaPsi
  )
  
  # ---- (4) Global latent template ----
  res_tau_old <- get_temp_spd(
    workGrid = obsGrid,
    Lspd_normd = norm_data,
    nknots = nknots,
    lambda = lambdaW,
    method = method,
    seed = seed
  )
  
  # ---- (5) Component-level warping functions ----
  res_Psi_old <- get_compwf_spd(
    workGrid = obsGrid,
    tau = res_tau_old$tau,
    Leta = res_eta_pre$Etahat,
    nknots = nknots,
    lambda = lambdaPsi
  )
  
  # ---- (6) Resolved component tempos ----
  res_eta_resolved <- get_comptempo_spd(
    workGrid = obsGrid,
    tau = res_tau$tau,
    Psi = res_Psi$H
  )
  
  res_G_old <- get_globwf_spd(
    res_H = res_H, 
    res_Psi = res_Psi_old, 
    workGrid = obsGrid
  )
  
  res_G_ref <- get_globwf_spd(
    res_H = res_H,
    res_Psi = res_Psi,
    workGrid = obsGrid
  )
  # ---- (7) Global alignment of trajectories ----
  Xhat_ginv_pre <- align_Xhat_Ginv_spd(
    workGrid = obsGrid,
    Xhat = Xhat_array,
    res_G= res_G_pre,
    method = method
  )

  Xhat_ginv_ref <- align_Xhat_Ginv_spd(
    workGrid = obsGrid,
    Xhat = Xhat_array,
    res_G= res_G_ref,
    method = method
  )
  
  Xhat_ginv_old <- align_Xhat_Ginv_spd(
    workGrid = obsGrid,
    Xhat = Xhat_array,
    res_G= res_G_old,
    method = method
  )
  # ---- (8) Return all results ----
  list(
    Ahat = Ahat,
    norm_data = norm_data,
    Xhat = Xhat_array,
    res_H = res_H,
    res_eta_pre = res_eta_pre,
    res_G_pre = res_G_pre,
    res_G_ref = res_G_ref,
    res_G_old = res_G_old,
    res_tau = res_tau,
    res_tau_old = res_tau_old,
    res_Psi = res_Psi,
    res_Psi_old = res_Psi_old,
    res_eta = res_eta_resolved,
    Xhat_ginv_pre = Xhat_ginv_pre,
    Xhat_ginv_ref = Xhat_ginv_ref,
    Xhat_ginv_old = Xhat_ginv_old,
    workGrid = obsGrid
  )
}



#' Fit GLDM for SPD trajectories — phase-only case (no amplitude variation)
#'
#' @description
#' Same as \code{gldm_spd_once} but skips amplitude estimation and normalization.
#' The smoothed trajectories Xhat are used directly as Xstar (norm_data).
#' Use this when A_{ij} = 1 for all (i,j), i.e., no amplitude variation.
#'
#' @inheritParams gldm_spd_once
#'
#' @return A list containing alignment results (same structure as gldm_spd_once
#'   but without amplitude estimation):
#' \itemize{
#'   \item \code{norm_data}: smoothed trajectories used as X* (= Xhat, no normalization).
#'   \item \code{Xhat}: smoothed trajectories (same as norm_data).
#'   \item \code{res_H}, \code{res_eta_pre}, \code{res_G_pre}, \code{res_G_ref},
#'         \code{res_G_old}, \code{res_tau}, \code{res_tau_old}, \code{res_Psi},
#'         \code{res_Psi_old}, \code{res_eta}: warping and template estimates.
#'   \item \code{Xhat_ginv_pre}, \code{Xhat_ginv_ref}, \code{Xhat_ginv_old}:
#'         globally aligned trajectories.
#'   \item \code{workGrid}: the time grid used.
#' }
#'
#' @export
gldm_spd_phaseonly <- function(obsGrid, Y, nknots = 4, lambdaH, lambdaW, lambdaG, lambdaPsi, method = "nearest", doLFR = TRUE, seed = NULL) {

  if (!is.list(Y) || length(Y) == 0L)
    stop("Y must be a non-empty list: Y[[j]][[i]] is a d×d×T1 array.")
  p <- length(Y)
  if (!is.list(Y[[1]]) || length(Y[[1]]) == 0L)
    stop("Each Y[[j]] must be a non-empty list of subjects.")

  n <- length(Y[[1]])
  d <- dim(Y[[1]][[1]])[1]

  # ---- (0) Xhat estimation (smoothing) ----
  m <- length(obsGrid)
  x_mat <- matrix(obsGrid, nrow = length(obsGrid))
  Xhat_array <- vector("list", p)
  for (j in 1:p) {
    Xhat_array[[j]] <- vector("list", n)
    for (i in 1:n) {
      Xhat_array[[j]][[i]] <- array(0, dim = c(d, d, m))
      if (!doLFR) {
        Xhat_array[[j]][[i]] <- Y[[j]][[i]]
      } else {
        y <- LocCovReg(x = x_mat, M = Y[[j]][[i]],
                       xout = x_mat,
                       optns = list(metric = "frobenius", kernel = "epan"))
        for (kk in 1:m) Xhat_array[[j]][[i]][, , kk] <- y$Mout[[kk]]
      }
    }
  }

  # ---- (1) Skip amplitude — use Xhat directly as Xstar ----
  norm_data <- Xhat_array

  # ---- (2) Subject-level warping ----
  res_H <- get_subjwf_spd(
    workGrid = obsGrid, Lspd_normd = norm_data,
    nknots = nknots, lambda = lambdaH
  )

  # ---- (3) Component-level alignment ----
  res_eta_pre <- align_comp_spd(
    workGrid = obsGrid, Lspd_normd = norm_data,
    Hinv = res_H$Hinv, method = method
  )

  res_G_pre <- align_glob_spd(
    workGrid = obsGrid, Lspd_normd = norm_data,
    nknots = nknots, lambda = lambdaG
  )

  # ---- (4) Global latent template (GLDM) ----
  res_tau <- get_template_spd(
    workGrid = obsGrid, Lspd_normd = res_G_pre$Lspd_normd_aligned
  )

  # ---- (5) Component-level warping functions (GLDM) ----
  res_Psi <- get_compwf_spd(
    workGrid = obsGrid, tau = res_tau$tau,
    Leta = res_eta_pre$Etahat, nknots = nknots, lambda = lambdaPsi
  )

  # ---- (4') Global latent template (CM) ----
  res_tau_old <- get_temp_spd(
    workGrid = obsGrid, Lspd_normd = norm_data,
    nknots = nknots, lambda = lambdaW, method = method, seed = seed
  )

  # ---- (5') Component-level warping functions (CM) ----
  res_Psi_old <- get_compwf_spd(
    workGrid = obsGrid, tau = res_tau_old$tau,
    Leta = res_eta_pre$Etahat, nknots = nknots, lambda = lambdaPsi
  )

  # ---- (6) Resolved component tempos ----
  res_eta_resolved <- get_comptempo_spd(
    workGrid = obsGrid, tau = res_tau$tau, Psi = res_Psi$H
  )

  res_G_old <- get_globwf_spd(res_H = res_H, res_Psi = res_Psi_old, workGrid = obsGrid)
  res_G_ref <- get_globwf_spd(res_H = res_H, res_Psi = res_Psi, workGrid = obsGrid)

  # ---- (7) Global alignment of trajectories ----
  Xhat_ginv_pre <- align_Xhat_Ginv_spd(
    workGrid = obsGrid, Xhat = Xhat_array, res_G = res_G_pre, method = method
  )
  Xhat_ginv_ref <- align_Xhat_Ginv_spd(
    workGrid = obsGrid, Xhat = Xhat_array, res_G = res_G_ref, method = method
  )
  Xhat_ginv_old <- align_Xhat_Ginv_spd(
    workGrid = obsGrid, Xhat = Xhat_array, res_G = res_G_old, method = method
  )

  # ---- (8) Return all results ----
  list(
    norm_data = norm_data,
    Xhat = Xhat_array,
    res_H = res_H,
    res_eta_pre = res_eta_pre,
    res_G_pre = res_G_pre,
    res_G_ref = res_G_ref,
    res_G_old = res_G_old,
    res_tau = res_tau,
    res_tau_old = res_tau_old,
    res_Psi = res_Psi,
    res_Psi_old = res_Psi_old,
    res_eta = res_eta_resolved,
    Xhat_ginv_pre = Xhat_ginv_pre,
    Xhat_ginv_ref = Xhat_ginv_ref,
    Xhat_ginv_old = Xhat_ginv_old,
    workGrid = obsGrid
  )
}

