# =============================================================================
# File: fctns_evaluation.R
# Purpose: Evaluation metrics (ISE/MISE) for GLDM and SRVF-based alignment
#          methods on SPD matrix-valued trajectories. Computes errors at all
#          levels: subject warp (H), template (tau), component warp (Psi), and
#          globally aligned matrix trajectories (X).
# Dependencies: pracma (trapz), SMFilter (FDist2)
#               Assumes FDist2() is defined in the global environment (via SMFilter).
# =============================================================================

# -----------------------------------------------------------------------------
# Individual error component functions
# -----------------------------------------------------------------------------

#' Compute per-subject ISE (Integrated squared error) and overall MISE (Mean of integrated squared error) for warping functions
#'
#' @param Hhat      n x T2 matrix: estimated forward warps evaluated on workGrid.
#' @param Htrue     n x T1 matrix: ground-truth forward warps evaluated on obsGrid.
#' @param obsGrid numeric vector of length T1: grid for Htrue (e.g., seq(0,1,by=0.05)).
#' @param workGrid  numeric vector of length T2: grid for Hhat (and integration).
#' @param scale     numeric scalar to rescale the final MISE (default 1).
#'
#' @return list with
#'   - Hise: length-n vector of ISEs
#'   - Hmise: scalar, mean(Hise) * scale
#'
#' @examples
#' # res <- compute_H_mise(Hhat = res_H$H, Htrue = simdata$hmat,
#' #                       obsGrid = seq(0,1,by=0.05), workGrid = workGrid, scale = 100)
#' # res$Hmise
compute_H_mise <- function(Hhat, Htrue, obsGrid, workGrid, scale = 1) {
  # ---- checks ----
  if (!is.matrix(Hhat) || !is.matrix(Htrue)) stop("Hhat and Htrue must be matrices.")
  if (nrow(Hhat) != nrow(Htrue)) stop("Hhat and Htrue must have the same number of subjects (rows).")
  if (!is.numeric(obsGrid) || !is.numeric(workGrid)) stop("obsGrid and workGrid must be numeric vectors.")
  if (ncol(Hhat) != length(workGrid)) stop("ncol(Hhat) must equal length(workGrid).")
  n <- nrow(Hhat)
  T2 <- length(workGrid)
  
  Hise <- numeric(n)
  for (i in seq_len(n)) {
    # interpolate true H to workGrid
    Htrue_dense <- approx(x = obsGrid, y = Htrue[i, ], xout = workGrid, rule = 2)$y
    diff_sq <- (Hhat[i, ] - Htrue_dense)^2
    Hise[i] <- pracma::trapz(workGrid, diff_sq)
  }
  list(Hise = Hise, Hmise = mean(Hise) * scale)
}

#' @title Compute ISE (Integrated squared error) between two matrix-valued trajectories 
#' @description
#' Compute the integrated squared Frobenius distance (ISE) between
#' estimated and true matrix-valued trajectories on a shared time grid:
#' \deqn{ISE = ∫ ||tau_true(t) - tau_hat(t)||_F^2 dt.}
#'
#' @param tau_hat 3D array (d × d × T): estimated trajectory.
#' @param tau_true 3D array (d × d × T): ground-truth trajectory.
#' @param workGrid numeric vector of length T: integration grid.
#' @param scale numeric scalar; optional scaling factor (default = 1).
#'
#' @return numeric scalar: ISE × scale.
#' @details
#' Requires an existing function `FDist2(A, B)` that computes
#' squared Frobenius distance between two matrices.
#'
#' @export
compute_T_ise <- function(tau_hat, tau_true, workGrid, scale = 1) {
  # ---- prerequisites ----
  if (!exists("FDist2", mode = "function"))
    stop("FDist2() not found. Please define FDist2(A, B) before calling this function.")
  
  # ---- dimension checks ----
  if (!is.array(tau_hat) || length(dim(tau_hat)) != 3L)
    stop("tau_hat must be a 3D array of size d × d × T.")
  if (!is.array(tau_true) || length(dim(tau_true)) != 3L)
    stop("tau_true must be a 3D array of size d × d × T.")
  
  d1 <- dim(tau_hat)[1]
  d2 <- dim(tau_hat)[2]
  Tlen <- dim(tau_hat)[3]
  if (d1 != d2)
    stop("tau_hat must be square in its first two dimensions (d × d × T).")
  if (any(dim(tau_true) != c(d1, d2, Tlen)))
    stop("tau_true must have the same dimensions as tau_hat (d × d × T).")
  
  if (!is.numeric(workGrid))
    stop("workGrid must be numeric.")
  if (length(workGrid) != Tlen)
    stop("length(workGrid) must equal T (third dimension of tau_hat).")
  
  # ---- compute ISE ----
  tdiff <- numeric(Tlen)
  for (j in seq_len(Tlen)) {
    tdiff[j] <- FDist2(tau_true[, , j], tau_hat[, , j])
  }
  ise <- pracma::trapz(workGrid, tdiff)
  
  as.numeric(ise) * scale
}

#' Compute mean ISE across components for SPD trajectories (Frobenius)
#'
#' @description
#' Given two component-first lists of SPD trajectories (each element is a d×d×T array),
#' compute the Integrated Squared Frobenius Error (ISE) per component on a shared grid
#' and return their average. Optionally returns the per-component ISE vector as well.
#'
#' @param Leta_hat list length p; each [[j]] is a d × d × T array (estimate).
#' @param Leta_true list length p; each [[j]] is a d × d × T array (ground truth).
#' @param workGrid numeric vector of length T (integration grid).
#' @param scale numeric; optional scaling factor applied to outputs (default 1).
#' @param return_per_comp logical; if TRUE, also return per-component ISE vector.
#'
#' @return If return_per_comp = FALSE (default), a scalar: mean ISE × scale.
#'         If TRUE, a list with:
#'           - mean_ise: scalar mean ISE × scale
#'           - ise_per_comp: numeric vector length p (each × scale)
#'
#' @details
#' Requires an existing function FDist2(A, B) that returns squared Frobenius distance.
#'
#' @export
compute_Leta_ise <- function(Leta_hat, Leta_true, workGrid, scale = 1, return_per_comp = FALSE) {
  # ---- prerequisites ----
  if (!exists("FDist2", mode = "function"))
    stop("FDist2() not found. Please define FDist2(A, B) before calling this function.")
  
  # ---- basic checks ----
  if (!is.list(Leta_hat) || length(Leta_hat) < 1)
    stop("Leta_hat must be a non-empty list of d×d×T arrays.")
  if (!is.list(Leta_true) || length(Leta_true) != length(Leta_hat))
    stop("Leta_true must be a list with the same length as Leta_hat.")
  
  p <- length(Leta_hat)
  
  # check shapes against the first element
  ref_hat  <- Leta_hat[[1]]
  ref_true <- Leta_true[[1]]
  if (!is.array(ref_hat)  || length(dim(ref_hat))  != 3L)
    stop("Each Leta_hat[[j]] must be a 3D array (d × d × T).")
  if (!is.array(ref_true) || length(dim(ref_true)) != 3L)
    stop("Each Leta_true[[j]] must be a 3D array (d × d × T).")
  
  d1 <- dim(ref_hat)[1]; d2 <- dim(ref_hat)[2]; Tlen <- dim(ref_hat)[3]
  if (d1 != d2) stop("Arrays must be square in the first two dimensions (d × d × T).")
  if (any(dim(ref_true) != c(d1, d2, Tlen)))
    stop("Leta_true[[1]] must match Leta_hat[[1]] in dimensions.")
  if (!is.numeric(workGrid) || length(workGrid) != Tlen)
    stop("workGrid must be numeric with length equal to T (third dimension).")
  
  # verify all components share the same shape
  for (j in seq_len(p)) {
    Ah <- Leta_hat[[j]]; At <- Leta_true[[j]]
    if (!is.array(Ah) || !is.array(At) || length(dim(Ah)) != 3L || length(dim(At)) != 3L)
      stop("All elements must be 3D arrays (d × d × T).")
    if (any(dim(Ah) != c(d1, d2, Tlen)) || any(dim(At) != c(d1, d2, Tlen)))
      stop("All arrays must share the same d × d × T dimensions across components.")
  }
  
  # ---- per-component ISE ----
  ise_per_comp <- numeric(p)
  for (j in seq_len(p)) {
    td <- numeric(Tlen)
    Ah <- Leta_hat[[j]]
    At <- Leta_true[[j]]
    for (k in seq_len(Tlen)) {
      td[k] <- FDist2(At[, , k], Ah[, , k])
    }
    ise_per_comp[j] <- pracma::trapz(workGrid, td)
  }
  
  mean_ise <- mean(ise_per_comp)
  
  if (isTRUE(return_per_comp)) {
    return(list(
      mean_ise     = as.numeric(mean_ise) * scale,
      ise_per_comp = as.numeric(ise_per_comp) * scale
    ))
  } else {
    return(as.numeric(mean_ise) * scale)
  }
}


#' Compute per-component ISE (Integrated squared error) and overall MISE (Mean of integrated squared error) for component-wise warps (Psi)
#'
#' @description
#' For each component k = 1..p, interpolate the ground-truth curve `Ptrue`
#' from `obsGrid` to `workGrid`, compute pointwise squared error against
#' the estimate `Phat`, and integrate with the trapezoidal rule to get ISE.
#' Returns the vector of per-component ISEs and the (scaled) mean across components (MISE).
#'
#' @param Phat       list (length p) or matrix (p x T2): estimated warps on `workGrid`.
#' @param Ptrue      list (length p) or matrix (p x T1): ground-truth warps on `obsGrid`.
#' @param obsGrid  numeric vector (length T1): grid for `Ptrue` (e.g., `seq(0,1,by=0.05)`).
#' @param workGrid   numeric vector (length T2): grid for `Phat` and for integration.
#' @param scale      numeric scalar to rescale the final MISE (default 1).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{Pise}: numeric vector (length p) of per-component ISEs.
#'   \item \code{Pmise}: numeric scalar, \code{mean(Pise) * scale}.
#' }
#'
#' @details
#' Integration is performed via \code{pracma::trapz}. Ensure the \code{pracma} package is available.
#'
#' @examples
#' \dontrun{
#' p <- 4
#' obsGrid <- seq(0, 1, by = 0.05)
#' workGrid  <- seq(0, 1, length.out = 201)
#' set.seed(1)
#' Ptrue <- lapply(seq_len(p), function(k) obsGrid + 0.05 * sin(2*pi*obsGrid + k))
#' Phat  <- lapply(seq_len(p), function(k) {
#'   approx(obsGrid, Ptrue[[k]], xout = workGrid)$y + rnorm(length(workGrid), sd = 0.01)
#' })
#' res <- compute_P_mise(Phat, Ptrue, obsGrid, workGrid, scale = 100)
#' res$Pmise
#' }
compute_P_mise <- function(Phat, Ptrue, obsGrid, workGrid, scale = 1) {
  # ---- normalize inputs to lists of numeric vectors (length p) ----
  to_list_by_row <- function(M) {
    if (is.matrix(M)) {
      split(M, row(M)) |> lapply(as.numeric)
    } else if (is.list(M)) {
      M
    } else stop("Phat/Ptrue must be either a list or a matrix.")
  }
  PhatL  <- to_list_by_row(Phat)
  PtrueL <- to_list_by_row(Ptrue)
  
  # basic checks
  if (!is.numeric(obsGrid) || !is.numeric(workGrid))
    stop("obsGrid and workGrid must be numeric vectors.")
  if (length(PhatL) != length(PtrueL))
    stop("Phat and Ptrue must have the same number of components (p).")
  
  p <- length(PhatL)
  # ensure Phat components match workGrid length
  for (k in seq_len(p)) {
    if (!is.numeric(PhatL[[k]]))
      stop("Each element of Phat must be a numeric vector.")
    if (length(PhatL[[k]]) != length(workGrid))
      stop("Each Phat[[k]] must have length equal to length(workGrid).")
    if (!is.numeric(PtrueL[[k]]))
      stop("Each element of Ptrue must be a numeric vector.")
  }
  
  # ---- compute ISE per component ----
  Pise <- numeric(p)
  for (k in seq_len(p)) {
    Ptrue_dense <- approx(x = obsGrid, y = PtrueL[[k]], xout = workGrid, rule = 2)$y
    diff_sq <- (PhatL[[k]] - Ptrue_dense)^2
    Pise[k] <- pracma::trapz(workGrid, diff_sq)
  }
  
  # ---- return ----
  list(Pise = Pise, Pmise = mean(Pise) * scale)
}


#' Compute ISE (Integrated squared error) and overall MISE (Mean of integrated squared error) for multivariate matrix-valued trajectories
#'
#' @description
#' For each component j = 1..p and subject i = 1..n, this function integrates
#' over `workGrid` the squared distance between two matrix-valued trajectories
#' `XGinv_hat[[j]][[i]][,,k]` and `XGinv_true[[j]][[i]][,,k]` using `FDist2`.
#' It returns the per-(i,j) ISE matrix and overall MISE (mean of all entries),
#' optionally rescaled by `scale`.
#'
#' @param XGinv_hat  list of length p; each element is a list of length n;
#'                   innermost elements are 3D arrays (d × d × T2) evaluated on `workGrid`.
#' @param XGinv_true same structure as `XGinv_hat`.
#' @param workGrid   numeric vector of length T2: integration grid.
#' @param scale      numeric scalar to rescale the final MISE (default 1).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{Xise}: numeric matrix (n × p), entry (i, j) is the ISE for subject i, component j.
#'   \item \code{Xmise}: numeric scalar, \code{mean(Xise) * scale}.
#' }
#'
#' @details
#' Requires an existing function \code{FDist2(A, B)} returning a nonnegative scalar distance
#' between two d × d matrices. Integration uses \code{pracma::trapz}.
compute_X_mise <- function(XGinv_hat, XGinv_true, workGrid, scale = 1) {
  if (!exists("FDist2", mode = "function"))
    stop("FDist2() not found. Please define FDist2(A, B) before calling this function.")
  if (!is.list(XGinv_hat) || !is.list(XGinv_true))
    stop("XGinv_hat and XGinv_true must be lists of length p (components).")
  if (!is.numeric(workGrid))
    stop("workGrid must be a numeric vector.")
  
  p <- length(XGinv_hat)
  if (length(XGinv_true) != p)
    stop("XGinv_hat and XGinv_true must have the same number of components (p).")
  
  n <- length(XGinv_hat[[1]])
  if (length(XGinv_true[[1]]) != n)
    stop("Each component must contain the same number of subjects (n).")
  
  T2 <- length(workGrid)
  Xise <- matrix(NA_real_, nrow = n, ncol = p)
  
  for (j in seq_len(p)) {
    for (i in seq_len(n)) {
      Xhat_ij  <- XGinv_hat[[j]][[i]]
      Xtrue_ij <- XGinv_true[[j]][[i]]
      
      if (!is.array(Xhat_ij) || !is.array(Xtrue_ij) ||
          length(dim(Xhat_ij)) != 3L || length(dim(Xtrue_ij)) != 3L)
        stop(sprintf("Elements must be 3D arrays (d x d x T2) at (j=%d, i=%d).", j, i))
      
      if (dim(Xhat_ij)[3] != T2 || dim(Xtrue_ij)[3] != T2)
        stop(sprintf("Time length mismatch at (j=%d, i=%d).", j, i))
      
      diff_vec <- numeric(T2)
      for (k in seq_len(T2))
        diff_vec[k] <- FDist2(Xhat_ij[,,k], Xtrue_ij[,,k])
      
      Xise[i, j] <- pracma::trapz(workGrid, diff_vec)
    }
  }
  
  list(Xise = Xise, Xmise = mean(Xise, na.rm = TRUE) * scale)
}

# -----------------------------------------------------------------------------
# Global warp MISE functions (gWMISE)
# -----------------------------------------------------------------------------

#' Compute gWMISE with primitive G estimate (includes R_{ij})
#'
#' @param G_hat_list List of p matrices (n x M): estimated G_{ij} on obsGrid.
#' @param simdata    Simulation data with hmat, psi, rmat.
#' @param obsGrid    Numeric vector (time grid).
#' @param scale      Rescaling factor (default 100).
#' @return Numeric scalar: mean ISE across all (i,j), times scale.
compute_gWMISE_p <- function(G_hat_list, simdata, obsGrid, scale = 100) {
  n <- nrow(simdata$hmat); p <- length(simdata$psi)
  trapz_int <- function(x, y) {
    nn <- length(x); sum(diff(x) * (y[-nn] + y[-1]) / 2)
  }
  total <- 0
  for (j in 1:p) for (i in 1:n) {
    s1 <- simdata$hmat[i, ]
    s2 <- approx(obsGrid, simdata$psi[[j]], xout = s1, rule = 2)$y
    s3 <- approx(obsGrid, simdata$rmat[[i]][[j]], xout = s2, rule = 2)$y
    total <- total + trapz_int(obsGrid, (G_hat_list[[j]][i, ] - s3)^2)
  }
  total / (n * p) * scale
}

#' Compute gWMISE with refined G estimate (excludes R_{ij})
#'
#' @param G_hat_list List of p matrices (n x M): estimated G_{ij} on obsGrid.
#' @param simdata    Simulation data with hmat, psi.
#' @param obsGrid    Numeric vector (time grid).
#' @param scale      Rescaling factor (default 100).
#' @return Numeric scalar: mean ISE across all (i,j), times scale.
compute_gWMISE_r <- function(G_hat_list, simdata, obsGrid, scale = 100) {
  n <- nrow(simdata$hmat); p <- length(simdata$psi)
  trapz_int <- function(x, y) {
    nn <- length(x); sum(diff(x) * (y[-nn] + y[-1]) / 2)
  }
  total <- 0
  for (j in 1:p) for (i in 1:n) {
    s1 <- simdata$hmat[i, ]
    s2 <- approx(obsGrid, simdata$psi[[j]], xout = s1, rule = 2)$y
    total <- total + trapz_int(obsGrid, (G_hat_list[[j]][i, ] - s2)^2)
  }
  total / (n * p) * scale
}

# -----------------------------------------------------------------------------
# High-level evaluation wrappers
# -----------------------------------------------------------------------------

#' Evaluate GLDM-SPD fit by computing ISE/MISE metrics
#'
#' @description
#' Given the fitted GLDM-SPD results (`gldmfit`) and the ground-truth simulation data (`simdata`),
#' this function computes alignment and warping errors at all levels:
#' subject-level (H), template-level (Tau), component-level (Psi), and
#' globally aligned matrix trajectories (X).
#'
#' @param ldmfit   list; output from \code{gldm_spd()} containing estimated objects.
#' @param simdata  list; ground-truth simulation data containing \code{hmat}, \code{tau_dense},
#'                 \code{psi}, and \code{X_aligned_Ginv}.
#' @param obsGrid numeric vector; grid used for ground-truth warps (e.g., seq(0,1,by=0.05)).
#' @param workGrid  numeric vector; working grid used in the estimation.
#'
#' @return A named list of metrics:
#' \itemize{
#'   \item \code{Hmise}: mean integrated squared error for subject-level warps.
#'   \item \code{Tise}: integrated squared error for the latent template.
#'   \item \code{Pmise}: mean integrated squared error for component-level warps.
#'   \item \code{Xmise}: mean integrated squared error for globally aligned matrices.
#' }
#'
#' @export
evaluate_gldm_spd <- function(ldmfit, simdata, obsGrid, workGrid) {
  # ---- compute subject-level warping MISE ----
  Hmise <- compute_H_mise(
    Hhat      = ldmfit$res_H$H,
    Htrue     = simdata$hmat,
    obsGrid = obsGrid,
    workGrid  = workGrid,
    scale     = 100
  )$Hmise
  
  # ---- compute template ISE ----
  Tise <- compute_T_ise(
    tau_hat  = ldmfit$res_tau$Tauhat,
    tau_true = simdata$tau_dense,
    workGrid = workGrid
  )
  
  # ---- compute component-level warping MISE ----
  Pmise <- compute_P_mise(
    Phat      = ldmfit$res_Psi$H,
    Ptrue     = simdata$psi,
    obsGrid = obsGrid,
    workGrid  = workGrid,
    scale     = 100
  )$Pmise
  
  # ---- compute global matrix alignment MISE ----
  Xmise <- compute_X_mise(
    XGinv_hat  = ldmfit$Xhat_ginv,
    XGinv_true = simdata$X_aligned_Ginv,
    workGrid   = workGrid,
    scale      = 1
  )$Pmise
  
  # ---- return all metrics ----
  list(
    Hmise = Hmise,
    Tise  = Tise,
    Pmise = Pmise,
    Xmise = Xmise
  )
}


#' Evaluate a single GLDM-SPD fit (full metric suite for one MC replicate)
#'
#' @description
#' Computes all ISE/MISE metrics comparing both GLDM and CM methods against
#' ground truth from one simulated dataset. Covers subject-level warping (H),
#' template (tau, new and old), component warps (Psi, new and old), component
#' tempos (eta), global warps (gWMISE), and globally aligned trajectories (X).
#'
#' @param gldmfit  List; output of \code{gldm_spd_once()}. Must contain
#'                 res_H, res_tau, res_tau_old, res_Psi, res_Psi_old,
#'                 res_eta_pre, res_eta, res_G_pre, res_G_ref, res_G_old,
#'                 Xhat_ginv_pre, Xhat_ginv_ref, Xhat_ginv_old.
#' @param simdata  List; output of \code{simulate_data()}. Must contain
#'                 hmat, tau, psi, Leta, X_Ginv, X_aligned.
#' @param obsGrid  Numeric vector; the common observation time grid.
#'
#' @return A named list of scalar metrics:
#'   Hmise, Tise_old, Tise_new, Pmise_old, Pmise_new, Emise_pre, Emise_final,
#'   gWMISE_cm, gWMISE_gldm_p, gWMISE_gldm_r,
#'   Xmise_Ginv_old, Xmise_Ginv_pre, Xmise_Ginv_final,
#'   Xmise_aligned_old, Xmise_aligned_p, Xmise_aligned_r.
evaluate_gldm_spd_once <- function(gldmfit, simdata, obsGrid) {
  # ---- compute subject-level warping MISE ----
  Hmise <- compute_H_mise(
    Hhat      = gldmfit$res_H$H,
    Htrue     = simdata$hmat,
    obsGrid = obsGrid,
    workGrid  = obsGrid,
    scale     = 100
  )$Hmise
  
  # ---- compute template ISE ----
  Tise_old <- compute_T_ise(
    tau_hat  = gldmfit$res_tau_old$tau,
    tau_true = simdata$tau,
    workGrid = obsGrid
  )
  Tise_new <- compute_T_ise(
    tau_hat  = gldmfit$res_tau$tau,
    tau_true = simdata$tau,
    workGrid = obsGrid
  )
  # ---- compute component-level warping MISE ----
  Pmise_old <- compute_P_mise(
    Phat      = gldmfit$res_Psi_old$H,
    Ptrue     = simdata$psi,
    obsGrid = obsGrid,
    workGrid  = obsGrid,
    scale     = 100
  )$Pmise
  
  Pmise_new <- compute_P_mise(
    Phat      = gldmfit$res_Psi$H,
    Ptrue     = simdata$psi,
    obsGrid = obsGrid,
    workGrid  = obsGrid,
    scale     = 100
  )$Pmise
  
  Emise_pre <- compute_Leta_ise(
    Leta_hat = gldmfit$res_eta_pre$Etahat,
    Leta_true = simdata$Leta,
    workGrid = obsGrid
  )
  
  Emise_final <- compute_Leta_ise(
    Leta_hat = gldmfit$res_eta,
    Leta_true = simdata$Leta,
    workGrid = obsGrid
  )
  
  # ---- compute global matrix alignment MISE ----
  Xmise_Ginv_old <- compute_X_mise(
    XGinv_hat  = gldmfit$Xhat_ginv_old,
    XGinv_true = simdata$X_Ginv,
    workGrid   = obsGrid,
    scale      = 1
  )$Xmise
  
  Xmise_Ginv_pre <- compute_X_mise(
    XGinv_hat  = gldmfit$Xhat_ginv_pre,
    XGinv_true = simdata$X_Ginv,
    workGrid   = obsGrid,
    scale      = 1
  )$Xmise
  
  Xmise_Ginv_final <- compute_X_mise(
    XGinv_hat  = gldmfit$Xhat_ginv_ref,
    XGinv_true = simdata$X_Ginv,
    workGrid   = obsGrid,
    scale      = 1
  )$Xmise
  
  Xmise_aligned_old <- compute_X_mise(
    XGinv_hat  = gldmfit$Xhat_ginv_old,
    XGinv_true = simdata$X_aligned,
    workGrid   = obsGrid,
    scale      = 1
  )$Xmise
  
  Xmise_aligned_p <- compute_X_mise(
    XGinv_hat  = gldmfit$Xhat_ginv_pre,
    XGinv_true = simdata$X_aligned,
    workGrid   = obsGrid,
    scale      = 1
  )$Xmise
  
  Xmise_aligned_r <- compute_X_mise(
    XGinv_hat  = gldmfit$Xhat_ginv_ref,
    XGinv_true = simdata$X_aligned,
    workGrid   = obsGrid,
    scale      = 1
  )$Xmise
  # ---- compute global warping MISE (gWMISE) ----
  gWMISE_cm         <- compute_gWMISE_p(gldmfit$res_G_old$H, simdata, obsGrid, scale = 100)
  gWMISE_gldm_p   <- compute_gWMISE_p(gldmfit$res_G_pre$H, simdata, obsGrid, scale = 100)
  gWMISE_gldm_r <- compute_gWMISE_r(gldmfit$res_G_ref$H, simdata, obsGrid, scale = 100)

  # ---- return all metrics ----
  list(
    Hmise = Hmise,
    Tise_old  = Tise_old,
    Tise_new  = Tise_new,
    Pmise_old = Pmise_old,
    Pmise_new = Pmise_new,
    Emise_pre = Emise_pre,
    Emise_final = Emise_final,
    gWMISE_cm = gWMISE_cm,
    gWMISE_gldm_p = gWMISE_gldm_p,
    gWMISE_gldm_r = gWMISE_gldm_r,
    Xmise_Ginv_old = Xmise_Ginv_old,
    Xmise_Ginv_pre = Xmise_Ginv_pre,
    Xmise_Ginv_final = Xmise_Ginv_final,
    Xmise_aligned_old = Xmise_aligned_old,
    Xmise_aligned_p = Xmise_aligned_p,
    Xmise_aligned_r = Xmise_aligned_r
  )
}