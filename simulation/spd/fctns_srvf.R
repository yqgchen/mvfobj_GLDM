# =============================================================================
# File: fctns_srvf.R
# Purpose: SRVF-based baseline methods for SPD functional data alignment.
#          Implements two Step-1 variants (Karcher-mean and pairwise alignment)
#          and a shared Steps 2-4 pipeline (template, component warps, G_ref).
#          Also provides an evaluation wrapper for SRVF fits.
# Dependencies: fdasrvf (multivariate_karcher_mean, multiple_align_multivariate,
#                        curve_pair_align)
#               fctns_estimation.R (align_Xhat_Ginv_spd, align_comp_spd,
#                                   get_template_spd, get_globwf_spd)
#               fctns_evaluation.R (compute_H_mise, compute_P_mise,
#                                   compute_T_ise, compute_gWMISE_p,
#                                   compute_gWMISE_r, compute_X_mise)
# =============================================================================

# -----------------------------------------------------------------------------
# SPD/PSD vectorization utility
# -----------------------------------------------------------------------------

#' Vectorize a symmetric matrix via its upper triangle (row-major order)
#'
#' @param mat A d x d symmetric matrix.
#' @return A numeric vector of length d*(d+1)/2 containing the upper triangular
#'   entries (including diagonal) in row-major order.
spd_to_vec <- function(mat) {
  d <- nrow(mat)
  idx <- which(upper.tri(mat, diag = TRUE), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  mat[idx]
}

# ==================================================================
# SRVF Step 1: H_i and G_ij_pre estimation
# ==================================================================

#' SRVF Step 1 via Karcher mean alignment
#'
#' @description
#' Estimate subject-level warps H_i and global warps G_{ij} (pre-refinement)
#' using multivariate Karcher mean and multiple_align_multivariate from fdasrvf.
#'
#' @param Xstar List of p components; each a list of n d x d x M arrays
#'              (amplitude-normalized trajectories).
#' @param obsGrid Numeric vector of length M (time grid).
#'
#' @return A list with res_H (subject warps) and res_G_pre (global warps).
srvf_step1_km <- function(Xstar, obsGrid) {
  p <- length(Xstar); n <- length(Xstar[[1]]); M <- length(obsGrid)
  d <- dim(Xstar[[1]][[1]])[1]
  vec_dim <- d * (d + 1) / 2

  # H_i: (vec_dim * p)-dim multivariate curve per subject
  beta_H <- array(0, c(vec_dim * p, M, n))
  for (i in 1:n) for (j in 1:p) for (k in 1:M) {
    offset <- (j - 1) * vec_dim
    beta_H[(offset + 1):(offset + vec_dim), k, i] <- spd_to_vec(Xstar[[j]][[i]][,,k])
  }
  # Karcher mean of all n curves in SRVF geometry (mode="O" for open curves);
  # then align every curve to the mean — out_H$gam[, i] is H_i^{-1} on obsGrid
  mu_H <- multivariate_karcher_mean(beta_H, mode = "O")
  out_H <- multiple_align_multivariate(beta_H, mu = mu_H$betamean,
             mode = "O", rotation = FALSE, scale = FALSE, lambda = 0, verbose = FALSE)

  H_hat <- Hinv_hat <- matrix(0, n, M)
  for (i in 1:n) {
    Hinv_hat[i, ] <- out_H$gam[, i]
    H_hat[i, ] <- approx(out_H$gam[, i], obsGrid, xout = obsGrid, rule = 2)$y
  }

  # G_ij: vec_dim-dim curves, n*p total
  N <- n * p
  beta_G <- array(0, c(vec_dim, M, N))
  ij_map <- matrix(0, N, 2); idx <- 0
  for (i in 1:n) for (j in 1:p) {
    idx <- idx + 1; ij_map[idx, ] <- c(i, j)
    for (k in 1:M) beta_G[, k, idx] <- spd_to_vec(Xstar[[j]][[i]][,,k])
  }
  # Same Karcher-mean alignment on all n*p (i,j)-pair curves;
  # out_G$gam[, idx] is the warp G_{ij}^{-1} on obsGrid for the pair encoded at idx
  mu_G <- multivariate_karcher_mean(beta_G, mode = "O")
  out_G <- multiple_align_multivariate(beta_G, mu = mu_G$betamean,
             mode = "O", rotation = FALSE, scale = FALSE, lambda = 0, verbose = FALSE)

  res_G_pre <- list(H = vector("list", p), Hinv = vector("list", p), workGrid = obsGrid)
  for (j in 1:p) {
    res_G_pre$H[[j]] <- matrix(0, n, M)
    res_G_pre$Hinv[[j]] <- matrix(0, n, M)
  }
  for (idx2 in 1:N) {
    i <- ij_map[idx2, 1]; j <- ij_map[idx2, 2]
    res_G_pre$Hinv[[j]][i, ] <- out_G$gam[, idx2]
    res_G_pre$H[[j]][i, ] <- approx(out_G$gam[, idx2], obsGrid, xout = obsGrid, rule = 2)$y
  }

  list(
    res_H = list(H = H_hat, Hinv = Hinv_hat, workGrid = obsGrid),
    res_G_pre = res_G_pre
  )
}

#' SRVF Step 1 via pairwise alignment
#'
#' @description
#' Estimate subject-level warps H_i and global warps G_{ij} (pre-refinement)
#' using pairwise curve_pair_align from fdasrvf, averaged across all pairs.
#'
#' @param Xstar List of p components; each a list of n d x d x M arrays
#'              (amplitude-normalized trajectories).
#' @param obsGrid Numeric vector of length M (time grid).
#'
#' @return A list with res_H (subject warps) and res_G_pre (global warps).
srvf_step1_pair <- function(Xstar, obsGrid) {
  p <- length(Xstar); n <- length(Xstar[[1]]); M <- length(obsGrid)
  d <- dim(Xstar[[1]][[1]])[1]
  vec_dim <- d * (d + 1) / 2

  # H_i: pairwise alignment of (vec_dim * p)-dim curves
  beta_H <- array(0, c(vec_dim * p, M, n))
  for (i in 1:n) for (j in 1:p) for (k in 1:M) {
    offset <- (j - 1) * vec_dim
    beta_H[(offset + 1):(offset + vec_dim), k, i] <- spd_to_vec(Xstar[[j]][[i]][,,k])
  }

  H_hat <- Hinv_hat <- matrix(0, n, M)
  for (i in 1:n) {
    gam_list <- matrix(0, n - 1, M)
    cnt <- 0
    for (r in 1:n) {
      if (r == i) next
      cnt <- cnt + 1
      # Align subject r's curve to subject i's; gam = warping that maps r's time to i's time
      out_pair <- curve_pair_align(beta_H[,,r], beta_H[,,i], mode = "O",
                                  rotation = FALSE, scale = FALSE)
      gam_list[cnt, ] <- out_pair$gam
    }
    avg_gam <- colMeans(gam_list)
    avg_gam <- (avg_gam - min(avg_gam)) / (max(avg_gam) - min(avg_gam))
    Hinv_hat[i, ] <- avg_gam
    H_hat[i, ] <- approx(avg_gam, obsGrid, xout = obsGrid, rule = 2)$y
  }

  # G_ij: pairwise alignment of vec_dim-dim curves
  N <- n * p
  beta_G <- array(0, c(vec_dim, M, N))
  ij_map <- matrix(0, N, 2); idx <- 0
  for (i in 1:n) for (j in 1:p) {
    idx <- idx + 1; ij_map[idx, ] <- c(i, j)
    for (k in 1:M) beta_G[, k, idx] <- spd_to_vec(Xstar[[j]][[i]][,,k])
  }

  G_fwd <- G_inv <- matrix(0, N, M)
  for (idx2 in 1:N) {
    gam_list <- matrix(0, N - 1, M)
    cnt <- 0
    for (ridx in 1:N) {
      if (ridx == idx2) next
      cnt <- cnt + 1
      # Align curve ridx to curve idx2; gam = warping of ridx to idx2's time axis
      out_pair <- curve_pair_align(beta_G[,,ridx], beta_G[,,idx2], mode = "O",
                                  rotation = FALSE, scale = FALSE)
      gam_list[cnt, ] <- out_pair$gam
    }
    avg_gam <- colMeans(gam_list)
    avg_gam <- (avg_gam - min(avg_gam)) / (max(avg_gam) - min(avg_gam))
    G_inv[idx2, ] <- avg_gam
    G_fwd[idx2, ] <- approx(avg_gam, obsGrid, xout = obsGrid, rule = 2)$y
  }

  res_G_pre <- list(H = vector("list", p), Hinv = vector("list", p), workGrid = obsGrid)
  for (j in 1:p) {
    res_G_pre$H[[j]] <- matrix(0, n, M)
    res_G_pre$Hinv[[j]] <- matrix(0, n, M)
  }
  for (idx2 in 1:N) {
    i <- ij_map[idx2, 1]; j <- ij_map[idx2, 2]
    res_G_pre$H[[j]][i, ] <- G_fwd[idx2, ]
    res_G_pre$Hinv[[j]][i, ] <- G_inv[idx2, ]
  }

  list(
    res_H = list(H = H_hat, Hinv = Hinv_hat, workGrid = obsGrid),
    res_G_pre = res_G_pre
  )
}

# ==================================================================
# SRVF Steps 2-4: Shared pipeline (tau, eta, Psi, G_ref)
# ==================================================================
#' SRVF pipeline: template, component warps, and refined global warps
#'
#' @description
#' Given Step 1 outputs (H_i and G_pre), estimate the latent template tau,
#' component-level warps Psi_j, refined global warps G_ref = Psi_j o H_i,
#' and align Xhat accordingly.
#'
#' @param Xstar  Amplitude-normalized trajectories (list of p lists of n arrays).
#' @param Xhat   Amplitude-scaled trajectories (list of p lists of n arrays).
#' @param obsGrid Numeric vector (time grid).
#' @param res_H  Subject-level warping result from Step 1.
#' @param res_G_pre Global warping result from Step 1.
#' @param method Interpolation method (default "nearest").
#'
#' @return A list containing res_H, res_G_pre, res_G_ref, res_tau, res_Psi,
#'         eta_pre, Xhat_ginv_pre, Xhat_ginv_ref.
srvf_pipeline <- function(Xstar, Xhat, obsGrid, res_H, res_G_pre,
                          method = "nearest") {
  p <- length(Xstar); n <- length(Xstar[[1]]); M <- length(obsGrid)
  d <- dim(Xstar[[1]][[1]])[1]
  vec_dim <- d * (d + 1) / 2

  # Step 2: tau = Frechet mean of X* aligned by G_pre^{-1}
  Xstar_aligned <- align_Xhat_Ginv_spd(
    workGrid = obsGrid, Xhat = Xstar, res_G = res_G_pre, method = method
  )
  tau_hat <- get_template_spd(
    workGrid = obsGrid, Lspd_normd_aligned = Xstar_aligned
  )$tau

  # Step 3a: eta_j = component-wise Frechet mean after H^{-1} alignment
  eta_pre <- align_comp_spd(
    workGrid = obsGrid, Lspd_normd = Xstar,
    Hinv = res_H$Hinv, method = method
  )

  # Step 3b: Psi_j via curve_pair_align(tau, eta_j)
  # Warping direction: gam = Psi^{-1}, so Psi = gam^{-1}
  Psi_hat <- PsiInv_hat <- vector("list", p)
  for (j in 1:p) {
    tau_vec <- eta_j_vec <- matrix(0, vec_dim, M)
    for (k in 1:M) {
      tau_vec[, k] <- spd_to_vec(tau_hat[,,k])
      eta_j_vec[, k] <- spd_to_vec(eta_pre$Etahat[[j]][,,k])
    }
    out_pair <- curve_pair_align(tau_vec, eta_j_vec, mode = "O",
                                rotation = FALSE, scale = FALSE)
    gam <- out_pair$gam
    PsiInv_hat[[j]] <- gam
    Psi_hat[[j]] <- approx(gam, obsGrid, xout = obsGrid, rule = 2)$y
  }
  res_Psi <- list(H = Psi_hat, Hinv = PsiInv_hat, workGrid = obsGrid)

  # Step 4: G_ref = Psi_j o H_i
  res_G_ref <- get_globwf_spd(
    res_H = res_H, res_Psi = res_Psi, workGrid = obsGrid
  )

  # Align Xhat (amplitude-scaled) for PMISE
  Xhat_ginv_pre <- align_Xhat_Ginv_spd(
    workGrid = obsGrid, Xhat = Xhat, res_G = res_G_pre, method = method
  )
  Xhat_ginv_ref <- align_Xhat_Ginv_spd(
    workGrid = obsGrid, Xhat = Xhat, res_G = res_G_ref, method = method
  )

  list(
    res_H = res_H,
    res_G_pre = res_G_pre,
    res_G_ref = res_G_ref,
    res_tau = list(tau = tau_hat),
    res_Psi = res_Psi,
    eta_pre = eta_pre,
    Xhat_ginv_pre = Xhat_ginv_pre,
    Xhat_ginv_ref = Xhat_ginv_ref
  )
}

# ==================================================================
# Evaluate one SRVF method
# ==================================================================
#' Evaluate SRVF fit against ground truth
#'
#' @param srvf_fit Output from srvf_pipeline().
#' @param simdata  Ground-truth simulation data.
#' @param obsGrid  Numeric vector (time grid).
#'
#' @return A named list of metrics: sWMISE, cWMISE, TISE,
#'         gWMISE_p, gWMISE_r, PMISE_pre, PMISE_ref.
evaluate_srvf <- function(srvf_fit, simdata, obsGrid) {
  sWMISE <- compute_H_mise(
    Hhat = srvf_fit$res_H$H, Htrue = simdata$hmat,
    obsGrid = obsGrid, workGrid = obsGrid, scale = 100
  )$Hmise

  cWMISE <- compute_P_mise(
    Phat = srvf_fit$res_Psi$H, Ptrue = simdata$psi,
    obsGrid = obsGrid, workGrid = obsGrid, scale = 100
  )$Pmise

  TISE <- compute_T_ise(
    tau_hat = srvf_fit$res_tau$tau, tau_true = simdata$tau,
    workGrid = obsGrid
  )

  gWMISE_p <- compute_gWMISE_p(srvf_fit$res_G_pre$H, simdata, obsGrid, scale = 100)
  gWMISE_r <- compute_gWMISE_r(srvf_fit$res_G_ref$H, simdata, obsGrid, scale = 100)

  PMISE_pre <- compute_X_mise(
    XGinv_hat = srvf_fit$Xhat_ginv_pre, XGinv_true = simdata$X_aligned,
    workGrid = obsGrid, scale = 1
  )$Xmise
  PMISE_ref <- compute_X_mise(
    XGinv_hat = srvf_fit$Xhat_ginv_ref, XGinv_true = simdata$X_aligned,
    workGrid = obsGrid, scale = 1
  )$Xmise

  list(sWMISE = sWMISE, cWMISE = cWMISE, TISE = TISE,
       gWMISE_p = gWMISE_p, gWMISE_r = gWMISE_r,
       PMISE_pre = PMISE_pre, PMISE_ref = PMISE_ref)
}
