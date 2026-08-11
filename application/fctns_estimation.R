# =============================================================================
# File: fctns_estimation.R
# Purpose: Core estimation functions for the GLDM mortality application,
#          including model fitting (ldm_distn), cross-sectional alignment
#          (get_pxcalign, get_sxcalign), and composite warping (get_compwf).
# Dependencies: purrr, frechet, fdapace, pracma, OPW, minqa
# =============================================================================

# Helper: atomically save an RDS file to path. -----------
safe_save_rds <- function(obj, path) {
  # Write to a temporary file in the same directory
  dir <- dirname(path)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- file.path(dir, paste0(".", basename(path), ".tmp"))  # hidden temp file
  
  saveRDS(obj, file = tmp, compress = "xz")  # write fully to tmp
  
  # On Windows, renaming over an existing file can fail → remove if exists
  if (file.exists(path)) file.remove(path)
  
  # Atomic-ish publish step (POSIX rename is atomic if same filesystem)
  ok <- file.rename(tmp, path)
  if (!ok) {
    # Fallback if rename fails (e.g., across filesystems)
    ok_copy <- file.copy(tmp, path, overwrite = TRUE)
    file.remove(tmp)
    if (!ok_copy) stop("Failed to publish RDS to: ", path)
  }
}

# compute the squared L2 distance between two functions on a common grid ----
get_sqL2 <- function(
    y1, # numeric vector; first function's values on \code{x}
    y2, # numeric vector; second function's values on \code{x}
    x   # numeric vector; common support grid
) {
  pracma::trapz(x, (y1 - y2)^2)
  # returns \int (y1 - y2)^2 dx approximated via the trapezoid rule
}

# resolve a character option: tries exact prefix match (pmatch), then fuzzy match (agrep);
# returns the matched element of valid_options. ----
resolve_option <- function(
    option,        # character string; the option supplied by the user
    valid_options  # character vector; the allowed option names
) {
  idx <- pmatch(option, valid_options)
  if (is.na(idx)) {
    idx <- agrep(option, valid_options, max.distance = 0.2)
  }
  if (length(idx) == 0) stop("Unknown option: ", option)
  if (length(idx) > 1) stop("Ambiguous option: ", option)
  valid_options[idx]
}

# get the reference point omega0 ----
# returns a list with qout (quantile function of omega0 on qSup),
# qSup, dout (density function of omega0), and dSup.
get_omega0 <- function (
  Lqf, # list of quantile functions per component per subject, 
  # e.g., \code{Lqf[[j]][[i]]} is a matrix with rows corresponding to time points and
  # columns corresponding to quantile values on \code{qSup} for subject \eqn{i} and component \eqn{j}.
  qSup, # vector of support points on \eqn{[0,1]} of quantile functions in \code{Lqf}.
  ref_source = "pooled" # character specifying the principal geodesics of which sets of objects (here, distributions) are used to determine the reference point: \code{"start"} (using objects at the first time point from each component, i.e., \eqn{p} sets of objects), \code{"end"} (using objects at the last time point from each component, i.e., \eqn{p} sets of objects), \code{"both"} (using objects at the first and last time points from each component, i.e., \eqn{2p} sets of objects), \code{"pooled"} (using objects at the first and last time points pooling all components together, i.e., 2 sets of objects
) {
  
  if ( any( !sapply( do.call(c,Lqf), is.matrix ) ) ) {
    stop ('Not all elements in Lqf are matrices.')
  }
  if ( length( unique( lapply( do.call(c,Lqf), dim ) ) ) > 1 ) {
    stop ('The dimensions of matrices in Lqf do not match.')
  }
  nqSup_Lqf <- ncol(Lqf[[1]][[1]])
  if ( missing(qSup) ) {
    qSup <- seq( 0, 1, length.out = nqSup_Lqf )
  } else {
    if ( length(qSup) != nqSup_Lqf ) {
      stop ("The length of qSup does not match the number of columns of Lqf[[j]][[i]].")
    }
  }
  
  ## prepare the list of matrices of quantile functions used to determine the reference point ----
  # Later, the first principal geodesic will be estimated for quantile functions in each matrix.
  valid_ref_sources <- c( "start", "end", "both", "pooled" )
  ref_source_resolved <- resolve_option(option = ref_source, valid_options = valid_ref_sources )
  Lqf4omega0 <- switch(
    ref_source_resolved,
    "start" = Lqf %>% llply( function ( res_comp ) {
      res_comp %>% laply( function (res_comp_subj) {
        head(res_comp_subj,1)
      })
    }),
    "end" = Lqf %>% llply( function ( res_comp ) {
      res_comp %>% laply( function (res_comp_subj) {
        tail(res_comp_subj,1)
      })
    }),
    "both" = c(
      Lqf %>% llply( function ( res_comp ) {
        res_comp %>% laply( function (res_comp_subj) {
          head(res_comp_subj,1)
        })
      }),
      Lqf %>% llply( function ( res_comp ) {
        res_comp %>% laply( function (res_comp_subj) {
          tail(res_comp_subj,1)
        })
      })
    ),
    "pooled" = list(
      Lqf %>% llply( function ( res_comp ) {
        res_comp %>% laply( function (res_comp_subj) {
          head(res_comp_subj,1)
        })
      }) %>% do.call( what = rbind ),
      Lqf %>% llply( function ( res_comp ) {
        res_comp %>% laply( function (res_comp_subj) {
          tail(res_comp_subj,1)
        })
      }) %>% do.call( what = rbind )
    )
  )
  
  ## estimate principal geodesic by applying FPCA on each set of quantile functions ----
  Lres_pg <- lapply( Lqf4omega0, function (mqf) {
    n_qf <- nrow(mqf)
    Lt <- replicate(n_qf, qSup, simplify = FALSE)
    Ly <- split( mqf, seq_len(n_qf) )
    res <- fdapace::FPCA(
      Ly = Ly, Lt = Lt,
      optns = list(
        dataType = "Dense",
        methodSelectK = "FVE",
        FVEthreshold = 0.95
      )
    )
    workGrid <- res$workGrid
    mu <- as.vector(res$mu)
    phi1 <- res$phi[,1]
    lambda1 <- res$lambda[1]
    
    ## when qSup and output grid from FPCA are not the same, obtain functions on qSup ----
    if ( !isTRUE( all.equal(qSup, workGrid, tolerance = 1e-6) ) ) {
      mu <- approx( x = workGrid, y = mu, xout = qSup, rule = 2)$y
      phi1 <- approx( x = workGrid, y = phi1, xout = qSup, rule = 2)$y
    }
    list(
      mu = mu,
      phi1 = phi1,
      lambda1 = lambda1
    )
  })
  # a list of FPCA results for each set of distributions,
  # of which each field is a list of estimates of mean function, 1st eigenfunction and 1st eigenvalue.
  
  ## find the range of support point s of principal geodesic s.t. the quantile function corresponding to g0(s) is non-decreasing ----
  find_s_range <- function(mu, phi, lambda) {
    if (length(mu) != length(phi)) stop("mu and phi must have the same length.")
    if (length(mu) < 2) stop("mu and phi must have length >= 2.")
    if (lambda < 0) stop("lambda must be nonnegative.")
    
    dmu <- diff(mu)
    dphi <- diff(phi)
    sqrt_lambda <- sqrt(lambda)
    
    # Handle zero lambda
    if (lambda == 0) {
      return(if (all(dmu >= 0)) c(-Inf, Inf) else numeric(0))
    }
    
    # Compute constraints
    lower_bounds <- -dmu[dphi > 0] / (sqrt_lambda * dphi[dphi > 0])
    upper_bounds <- -dmu[dphi < 0] / (sqrt_lambda * dphi[dphi < 0])
    
    s_lower <- if (length(lower_bounds)) max(lower_bounds) else -Inf
    s_upper <- if (length(upper_bounds)) min(upper_bounds) else Inf
    
    if (any(dphi == 0 & dmu < 0)) return(numeric(0))  # impossible case
    
    if (s_lower <= s_upper) c(s_lower, s_upper) else numeric(0)
  }
  
  ## find the ranges of support grids of principal geodesics ----
  LgSupRange <- lapply( Lres_pg, function (res) {
    with(res, find_s_range( mu = mu, phi = phi1, lambda = lambda1 ) )
  })
  if ( any( sapply(LgSupRange,length) != 2 ) ) {
    stop ("Principal geodesic given by Log-PCA does not have well-defined quantile functions.")
  }
  
  ## function to get points on a principal geodesic ---
  get_PG <- function (gSup, mu, lambda, phi) {
    mu + sqrt(lambda) * phi %*% t(gSup) 
  }
  
  ## prepare a list of lists of objects along each principal geodesic ----
  Lpg <- lapply( seq_along(LgSupRange), function(k) {
    gSup <- seq( LgSupRange[[k]][1], LgSupRange[[k]][2], length.out = 100 )
    mpg <- get_PG(
      gSup = gSup,
      mu = Lres_pg[[k]]$mu, 
      lambda = Lres_pg[[k]]$lambda1,
      phi = Lres_pg[[k]]$phi1
    )
    ## rows -> qSup
    ## columns -> gSup
    split(t(mpg), seq_along(gSup)) # list of quantile functions
  })
  
  ## function to find (by brute force) the combination of points on each principal geodesic such that the sum of squared distances is minimized ----
  find_min_combo <- function(lists, sqDist = function (x,y) { pracma::trapz( qSup, (x-y)^2 ) }) {
    p <- length(lists)
    idx_grid <- expand.grid(lapply(lists, seq_along))
    
    best_val <- Inf
    best_combo <- NULL
    
    for (i in seq_len(nrow(idx_grid))) {
      combo <- Map(`[[`, lists, idx_grid[i, ] )
      dists <- combn(p, 2, function(idx) sqDist(combo[[idx[1]]], combo[[idx[2]]]) )
      val <- sum(dists)
      
      if (val < best_val) {
        best_val <- val
        best_combo <- combo
      }
    }
    
    list(min_value = best_val, best_combo = best_combo)
  }
  ## function to approximate (by local search) the combination of points on each principal geodesic such that the sum of squared distances is minimized ----
  find_min_combo_approx <- function (
    lists,
    sqDist = function (x,y) { pracma::trapz( qSup, (x-y)^2 ) },
    max_iter = 500
  ) {
    p <- length(lists)
    combo <- lapply(lists, function(L) L[[1]])  # init arbitrarily

    for (t in seq_len(max_iter)) {
      changed <- FALSE
      for (i in seq_len(p)) {
        others <- combo[-i]
        scores <- sapply(lists[[i]], function(x)
          sum(sapply(others, function(y) sqDist(x, y))))
        best_idx <- which.min(scores)
        if (!identical(combo[[i]], lists[[i]][[best_idx]])) {
          combo[[i]] <- lists[[i]][[best_idx]]
          changed <- TRUE
        }
      }
      if (!changed) break
    }

    total_val <- sum(combn(p, 2, function(idx) sqDist(combo[[idx[1]]], combo[[idx[2]]])))
    list(min_value = total_val, best_combo = combo)
  }
  
  ## find the best combination ----
  res_combo <- find_min_combo( lists = Lpg )
  # res_combo <- find_min_combo_approx( lists = Lpg )
  # res_combo$best_combo: a list of vectors holding the quantile functions from each principal geodesic such that the sum of squared distances is minimzed
  
  ## return omega0 as the Fréchet mean of the best combination ----
  optns$qSup <- qSup
  res_omega0 <- frechet::DenFMean(
    qin = do.call(rbind, res_combo$best_combo),
    optns = optns
  )
  
  return(list(
    qout = res_omega0$qout, # vector holding the quantile function of \eqn{\omega_0} evaluated on \code{qSup}
    qSup = qSup, # vector holding the support grid of quantile function of \eqn{\omega_0}
    dout = res_omega0$dout, # vector holding the density function of \eqn{\omega_0} evaluated on \code{dSup}
    dSup = res_omega0$dSup # vector holding the support grid of density function of \eqn{\omega_0}
  ))
}

# get amplitude factors ----
get_amp <- function ( 
    Lqf, # list of quantile functions per component per subject, 
         # e.g., \code{Lqf[[j]][[i]]} is a matrix with rows corresponding to time points and
         # columns corresponding to quantile values on \code{qSup} for subject \eqn{i} and component \eqn{j}.
    qf_omega0, # quantile function of \eqn{\omega_0} evaluated on \code{qSup}
    qSup # support grid of quantile functions in \code{Lqf} and \code{qf_omega0}
) {
  sup_dist_to_omega0 <- sapply( Lqf, function ( res_comp ) {
    sapply( res_comp, function ( res_comp_subj ) {
      max( sqrt( apply(res_comp_subj, 1, get_sqL2, y2 = qf_omega0, x = qSup ) ) )
    })
  })
  c0 <- min(sup_dist_to_omega0)
  
  # output a matrix with columns corresponding to components and rows to subjects
  sup_dist_to_omega0 / c0
}

# get normalized processes ----
normalize <- function (
    Lqf, # list of quantile functions per component per subject, 
    # e.g., \code{Lqf[[j]][[i]]} is a matrix with rows corresponding to time points and
    # columns corresponding to quantile values on a grid in \eqn{[0,1]} for subject \eqn{i} and component \eqn{j}.
    qf_omega0, # quantile function of \eqn{\omega_0} evaluated on the same grid in \eqn{[0,1]} as those in \code{Lqf}.
    A # A matrix holding amplitude factors with columns corresponding to components and rows to subjects
) {
  res <- lapply( seq_along(Lqf), function (j) {
    res_comp <- lapply( seq_along(Lqf[[j]]), function (i) {
      AijInv <- 1/A[i,j]
      t( qf_omega0 * (1-AijInv) + t(Lqf[[j]][[i]]) * AijInv )
    })
    names(res_comp) <- names(Lqf[[j]])
    res_comp
  })
  names(res) <- names(Lqf)
  
  res # output a list of quantile functions per component per subject for the normalized processes, with the same structure as \code{Lqf}.
}

# get subject-level deformation functions ----

get_subjwf <- function (
    obsGrid, # vector of time points at which normalized processes are input.
    Lqf_normd, # List of quantile functions per component per subject for the normalized processes, 
               # e.g., \code{Lqf_normd[[j]][[i]]} is a matrix with rows corresponding to time points in \code{obsGrid} and
               # columns corresponding to quantile values on \code{qSup} for subject \eqn{i} and component \eqn{j}.
    qSup, # support grid on \eqn{[0,1]} of quantile functions in \code{Lqf_normd}.
    nknots = 4, # number of knots; default: 4.
    lambda # vector of regularization parameters for each component; default: integrated Fréchet variance times \eqn{10^{-4}}.
) {
  n_comp <- length( Lqf_normd )
  n_subj <- sapply( Lqf_normd, length )
  if ( any( abs( diff(n_subj) ) > 0 ) ) {
    stop ( "Numbers of subjects are not the same across components." )
  }
  n_subj <- n_subj[1]
  
  missing_lambda <- missing( lambda )
  if ( !missing_lambda ) {
    len_lambda <- length( lambda )
    if ( len_lambda == 1 ) {
      lambda = rep( lambda, n_comp )
      message ( "Only one value of the regularization paramter is input---applied to all components." )
    } else if ( len_lambda != n_comp ) {
      stop ( "Length of input lambda does not match the number of components." )
    }
  }
  
  Lres <- lapply( seq_len( n_comp ), function (j) {
    dat_j <- Lqf_normd[[j]]
    if ( missing_lambda ) {
      
      # compute default choice of regularization paramter
      lambda_j <- 1e-4 * pracma::trapz(
        obsGrid, 
        apply(
          simplify2array( dat_j ), 1, 
          function (qmat) {
            frechet::DenFVar(
              qin = t(qmat),
              supin = qSup
            )$DenFVar
          }
        )
      )
      
    } else {
      lambda_j <- lambda[j]
    }
    WassPWdense(
      tVec = obsGrid, qf = dat_j, qSup = qSup,
      optns = list(
        lambda = lambda_j,
        nknots = nknots
      )
    )
  })
  names(Lres) <- names(Lqf_normd)
  
  # take the average across components to get the final estimate of subject-level deformation functions
  H <- Reduce( '+', lapply( Lres, function (res) res$hInv ) ) / n_comp
  Hinv <- t( apply( H, 1, function (h) {
    approx( h, obsGrid, xout = obsGrid, rule = 2 )$y
  }) )
  
  list(
    H = H, # matrix holding the estimates of subject-level deformation functions, with rows corresponding to subjects and columns to time points in \code{obsGrid}
    Hinv = Hinv, # matrix holding the inverses of estimated subject-level deformation functions, with rows corresponding to subjects and columns to time points in \code{obsGrid}
    obsGrid = obsGrid, # time grid on which subject-level deformation functions and their inverses are evaluated
    Loptns = lapply( Lres, `[[`, "optns" ), # list of control options used for each component
    LtimingWarp = lapply( Lres, `[[`, "timingWarp" ) # list of time cost for each component
  )
}

# get indices of time points in a given grid closest to a desired time grid;
# returns an integer vector of the same length as tout. ----
match_times <- function (
    tin, # vector of input time points
    tout # vector of desired output time points 
) {
  idx <- findInterval(tout, tin)
  # Compare which side is closer
  idx <- ifelse(
    idx == 0, 1,
    ifelse(idx == length(tin), length(tin),
           ifelse(abs(tout - tin[idx]) < abs(tout - tin[idx + 1]), idx, idx + 1))
  )
  idx
}

# align per component using subject-level warping functions to obtain (primitive) component tempo trajectories ----

align_comp <- function (
    obsGrid, # vector of time points at which normalized processes are input.
    Lqf_normd,  # List of quantile functions per component per subject for the normalized processes, 
    # e.g., \code{Lqf_normd[[j]][[i]]} is a matrix with rows corresponding to time points in \code{obsGrid} and
    # columns corresponding to quantile values on \code{qSup} in \eqn{[0,1]} for subject \eqn{i} and component \eqn{j}.
    qSup, # vector of support points of quantile functions in \code{Lqf_normd}.
    Hinv, # matrix holding the inverses of estimated subject-level warping functions, with rows corresponding to subjects and columns to time points in \code{obsGrid}.
    optns = list() # list of optional control parameters: \code{bwDen}, \code{ndSup}, \code{dSup}, \code{delta}, \code{kernelDen}, \code{infSupport}, and \code{denLowerThreshold}. See \code{\link[frechet]{DenFMean}()} for details.
) {
  optns$qSup <- qSup
  
  res <- lapply( Lqf_normd, function (dat_comp) {
    n_subj <- nrow(Hinv)
    res_comp <- apply(
      simplify2array( 
        lapply( seq_len(n_subj), function (i) {
          idx <- match_times( tin = obsGrid, tout = Hinv[i,] )
          dat_comp[[i]][ idx, , drop = FALSE ]
        })
      ),
      1, function (qmat) {
        frechet::DenFMean(
          qin = t(qmat),
          optns = optns
        )
      }
    )
    list(
      qout = do.call(rbind, lapply( res_comp, `[[`, "qout" ) ),
      qSup = res_comp[[1]]$qSup,
      dout = do.call(rbind, lapply( res_comp, `[[`, "dout" ) ),
      dSup = do.call(rbind, lapply( res_comp, `[[`, "dSup" ) )
    )
  })
  
  ## re-organize by types of results
  list(
    qout = lapply( res, `[[`, "qout" ), # a list of matrices for each component, each matrix having rows corresponding to time points in \code{obsGrid} and columns corresponding to quantile values on \code{qSup}
    qSup = res[[1]]$qSup, # vector holding the support grid of quantile functions
    dout = lapply( res, `[[`, "dout" ), # a list of matrices for each component, each matrix having rows corresponding to time points in \code{obsGrid} and columns corresponding to density functions evaluated on corresponding grids in \code{dSup}
    dSup = lapply( res, `[[`, "dSup" ) # a list of matrices for each component, each matrix having rows corresponding to time points in \code{obsGrid} and columns giving support grids of density functions in \code{dout}
  )
}

# get primitive estimates of global deformation functions ----

align_glob <- function (
    obsGrid, # vector of time points at which normalized processes are input.
    Lqf_normd, # List of quantile functions per component per subject for the normalized processes, 
    # e.g., \code{Lqf_normd[[j]][[i]]} is a matrix with rows corresponding to time points in \code{obsGrid} and
    # columns corresponding to quantile values on \code{qSup} for subject \eqn{i} and component \eqn{j}.
    qSup, # support grid on \eqn{[0,1]} of quantile functions in \code{Lqf_normd}.
    nknots = 4, # number of knots; default: 4.
    lambda # vector of regularization parameters for each component; default: integrated Fréchet variance times \eqn{10^{-4}}.
) {
  n_comp <- length( Lqf_normd )
  n_subj <- sapply( Lqf_normd, length )
  if ( any( abs( diff(n_subj) ) > 0 ) ) {
    stop ( "Numbers of subjects are not the same across components." )
  }
  n_subj <- n_subj[1]
  
  missing_lambda <- missing( lambda )
  if ( !missing_lambda ) {
    len_lambda <- length( lambda )
    if ( len_lambda > 1 ) {
      lambda = lambda[1]
      warning ( "The input lambda has more than one elements. Only the first value of lambda is used." )
    } else if ( len_lambda < 1 ) {
      warning ( "Length of input lambda is less than 1---reset as default." )
      missing_lambda <- TRUE
    }
  }
  
  # pool all (i,j) pairs and obtain a list of (number of subjects)(number of components) fields each holding the quantile functions of the normalized process \eqn{X^*_{i,j}}
  dat <- do.call( c, Lqf_normd )
  
  # obtain choice of regularization parameter
  if ( missing_lambda ) {
    lambda <- 1e-4 * pracma::trapz(
      obsGrid, 
      apply(
        simplify2array( dat ), 1, 
        function (qmat) {
          frechet::DenFVar(
            qin = t(qmat),
            supin = qSup
          )$DenFVar
        }
      )
    )
    
  } 
  
  # apply pairwise warping
  res <- WassPWdense(
    tVec = obsGrid, qf = dat, qSup = qSup,
    optns = list(
      lambda = lambda,
      nknots = nknots
    )
  )
  
  
  # structure the results according to components
  H <- lapply( seq_len(n_comp), function (j) {
    res$hInv[(j-1)+seq_len(n_subj),,drop=FALSE]
  })
  Hinv <- lapply( seq_len(n_comp), function (j) {
    res$h[(j-1)+seq_len(n_subj),,drop=FALSE]
  })
  Lqf_normd_aligned <- lapply( seq_len(n_comp), function (j) {
    res_comp <- res$qfAligned[(j-1)+seq_len(n_subj)]
    names(res_comp) <- names(Lqf_normd[[j]])
    res_comp
  })
  names(Lqf_normd_aligned) <- names(H) <- names(Hinv) <- names(Lqf_normd)
  
  list(
    H = H, # list of matrices for each component holding the estimates of global deformation functions, with rows corresponding to subjects and columns to time points in \code{obsGrid}, e.g., \code{H[[j]][i,]} holds the values of estimates of \eqn{G_{i,j}(t)} for \eqn{t} in \code{obsGrid}.
    Hinv = Hinv, # list of matrices for each component holding the inverses of estimated global deformation functions, with rows corresponding to subjects and columns to time points in \code{obsGrid}.
    Lqf_normd_aligned = Lqf_normd_aligned, # list of quantile functions per component per subject for the globally aligned normalized processes, 
    # e.g., \code{Lqf_normd_aligned[[j]][[i]]} is a matrix with rows corresponding to time points in \code{obsGrid} and
    # columns corresponding to quantile values on \code{qSup} for subject \eqn{i} and component \eqn{j}.
    obsGrid = obsGrid, # Vector holding the time grid on which global deformation functions and their inverses are evaluated
    optns = res$optns, # list of control options used for time warping
    timingWarp = res$timingWarp # time cost of time warping
  )
}

# get the latent template ----

get_template <- function (
    obsGrid, # vector of time points at which normalized processes are input.
    Lqf_normd_aligned,  # List of quantile functions per component per subject for the globally aligned normalized processes, 
    # e.g., \code{Lqf_normd_aligned[[j]][[i]]} is a matrix with rows corresponding to time points in \code{obsGrid} and
    # columns corresponding to quantile values on \code{qSup} in \eqn{[0,1]} for subject \eqn{i} and component \eqn{j}.
    qSup, # vector of support points of quantile functions in \code{Lqf_normd}.
    optns = list() # list of optional control parameters: \code{bwDen}, \code{ndSup}, \code{dSup}, \code{delta}, \code{kernelDen}, \code{infSupport}, and \code{denLowerThreshold}. See \code{\link[frechet]{DenFMean}()} for details.
) {
  optns$qSup <- qSup
  
  # estimate the template by the Fréchet mean of aligned normalized processes
  res <- apply( simplify2array( do.call( c, Lqf_normd_aligned ) ), 1, function (qmat) {
    frechet::DenFMean(
      qin = t(qmat),
      optns = optns
    )
  })
  
  list(
    qout = do.call(rbind, lapply( res, `[[`, "qout" ) ), # a matrix holding the values of quantile functions of the template with rows corresponding to time points in \code{obsGrid} and columns corresponding to quantile values on \code{qSup}
    qSup = optns$qSup, # a vector holding the support grid of quantile functions
    dout = do.call(rbind, lapply( res, `[[`, "dout" ) ), # a matrix holding the values of density functions of the template with rows corresponding to time points in \code{obsGrid} and columns corresponding to density functions evaluated on corresponding grids in \code{dSup}
    dSup = do.call(rbind, lapply( res, `[[`, "dSup" ) ) # a matrix holding the support grids of density functions of the template with rows corresponding to time points in \code{obsGrid} and columns giving support grids of density functions in \code{dout}
  )
}

# compute composition of \eqn{f_1\circ f_2} ----
composite <- function(
    x, # vector of support grid points of \code{f1}
    f1, # values of function \eqn{f_1} on \code{x}
    f2 # values of function \eqn{f_2} on certain grid
) {
  approx( x = x, y = f1, xout = f2, rule = 2 )$y
  # output a vector holding the values of \eqn{f_1\circ f_2} on the support grid of \code{f2}
}

# get component-level warping functions by pairwise warping ----
get_compwf <- function (
    obsGrid, # vector of time points at which normalized processes are input.
    qf_tau, # matrix holding the quantile functions of the latent template 
    # with rows corresponding to time points in \code{obsGrid} and 
    # columns corresponding to quantile values on \code{qSup}.
    Lqf_eta,  # list of matrices for each component, 
    # each holding the quantile functions of one component tempo
    # with rows corresponding to time points in \code{obsGrid} and 
    # columns corresponding to quantile values on \code{qSup}.
    qSup, # support grid on \eqn{[0,1]} of quantile functions in \code{Lqf_normd}.
    nknots = 4, # number of knots; default: 4.
    lambda # regularization parameter; default: integrated Fréchet variance times \eqn{10^{-4}}.
) {
  missing_lambda <- missing(lambda)
  if ( !missing_lambda ) {
    len_lambda <- length( lambda )
    if ( len_lambda > 1 ) {
      lambda = lambda[1]
      warning ( "The input lambda has more than one elements. Only the first value of lambda is used." )
    } else if ( len_lambda < 1 ) {
      warning ( "Length of input lambda is less than 1---reset as default." )
      missing_lambda <- TRUE
    }
  }
  n_comp <- length(Lqf_eta)
  if ( missing_lambda ) {
    lambda <- 1e-4 * pracma::trapz(
      obsGrid, 
      apply(
        Reduce(
          '+',
          lapply( Lqf_eta, function ( qf_etaj ) {
            ( qf_etaj - qf_tau )^2
          })
        ) / n_comp,
        1, pracma::trapz, x = qSup
      )
    )
  }
  
  Lres <- lapply( Lqf_eta, function ( qf_etaj ) {
    res_comp <- WassPWdense(
      tVec = obsGrid, qf = list( qf_tau, qf_etaj ), qSup = qSup,
      optns = list(
        lambda = lambda,
        nknots = nknots
      )
    )
    
    list(
      h = res_comp$hInv[2,],
      hInv = res_comp$hInv[1,],
      optns = res_comp$optns,
      timingWarp = res_comp$timingWarp
    )
  })
  
  list(
    H = lapply( Lres, `[[`, "h"), # list of vectors for each component holding the estimates of component-level warping functions on time points in \code{obsGrid}
    Hinv = lapply( Lres, `[[`, "hInv"), # list of vectors for each component holding the inverses of estimated component-level warping functions on time points in \code{obsGrid}
    obsGrid = obsGrid, # time grid on which component-level warping functions and their inverses are evaluated
    Loptns = lapply( Lres, `[[`, "optns" ), # list of control options used for each component
    LtimingWarp = lapply( Lres, `[[`, "timingWarp" ) # list of time cost for each component
  )
}

# get component tempo trajectories ----
# returns a list with:
#   qout: list of matrices per component, rows = time points, cols = quantile values on qSup
#   qSup: support grid of quantile functions
#   dout: list of matrices per component, rows = time points, cols = density values
#   dSup: list of matrices per component, rows = time points, cols = density support grids
get_comptempo <- function (
    obsGrid, # vector of time points at which normalized processes are input.
    res_tau, # output from \code{\link{get_template}()}.
    Psi # list of vectors for each component holding the estimates of component-level deformation functions on time points in \code{obsGrid}
) {
  res <- lapply( Psi, function ( Psi_j ) {
    idx <- match_times( obsGrid, Psi_j )
    list(
      qout = res_tau$qout[ idx, , drop = FALSE ],
      dout = res_tau$dout[ idx, , drop = FALSE ],
      dSup = res_tau$dSup[ idx, , drop = FALSE ]
    )
  })
  
  list(
    qout = lapply( res, `[[`, "qout" ),
    qSup = res_tau$qSup,
    dout = lapply( res, `[[`, "dout" ),
    dSup = lapply( res, `[[`, "dSup" )
  )
}


# get (final estimates of) global deformation functions ----
get_globwf <- function (
    res_H, # output from \code{\link{get_subjwf}()}.
    res_Psi, # output from \code{\link{get_compwf}()}
    obsGrid # vector of support grid points for output global deformation functions
) {
  obsGrid_Psi <- res_Psi$obsGrid
  obsGrid_H <- res_H$obsGrid
  if ( isTRUE( all.equal( obsGrid, obsGrid_Psi ) ) ) {
    PsiInv <- res_Psi$Hinv
  } else {
    PsiInv <- lapply( res_Psi$Hinv, function (Psi_j_inv) {
      approx( x = obsGrid_Psi, y = Psi_j_inv, xout = obsGrid, rule = 2 )$y
    })
  }
  if ( isTRUE( all.equal( obsGrid, obsGrid_H ) ) ) {
    H <- res_H$H
  } else {
    H <- t( apply( res_H$H, 1, function (H_i) {
      approx( x = obsGrid_H, y = H_i, xout = obsGrid, rule = 2 )$y
    }) )
  }
  Hinv <- res_H$Hinv
  Psi <- res_Psi$H
  
  H <- lapply( Psi, function (Psi_j) {
    t(apply( H, 1, function (H_i) {
      composite( x = obsGrid_Psi, f1 = Psi_j, f2 = H_i )
    }))
  })
  Hinv <- lapply( PsiInv, function (Psi_j_inv) {
    t(apply( Hinv, 1, function (H_i_inv) {
      composite( x = obsGrid_H, f1 = H_i_inv, f2 = Psi_j_inv )
    }))
  })
  
  list(
    H = H, # list of matrices for each component holding the estimates of global deformation functions, with rows corresponding to subjects and columns to time points in \code{obsGrid}, e.g., \code{H[[j]][i,]} holds the values of estimates of \eqn{G_{i,j}(t)} for \eqn{t} in \code{obsGrid}.
    Hinv = Hinv, # list of matrices for each component holding the inverses of estimated global deformation functions, with rows corresponding to subjects and columns to time points in \code{obsGrid}.
    obsGrid = obsGrid # copy of the input \code{obsGrid}---a vector holding the time grid on which global deformation functions and their inverses are evaluated.
  )
}

# fit a LDM for distributional data ----
# returns a list with all intermediate and final estimation results:
#   res_omega0:  reference point omega0 (output of get_omega0)
#   A:           n x p amplitude factor matrix
#   Lqf_normd:   normalized processes (same structure as Lqf)
#   res_H:       subject-level deformation functions (output of get_subjwf)
#   res_eta_pre: primitive component tempo trajectories (output of align_comp)
#   res_G_pre:   primitive global deformation functions (output of align_glob)
#   res_tau:     latent template (output of get_template)
#   res_Psi:     component-level deformation functions (output of get_compwf)
#   res_eta:     refined component tempo trajectories (output of get_comptempo)
#   res_G:       refined global deformation functions (output of get_globwf)
#   obsGrid, qSup: copies of the inputs
#   timing:      elapsed computation time
ldm_distn <- function (
  obsGrid, # vector of time points at which processes in \code{Lqf} are evaluated.
  Lqf, # list of quantile functions per component per subject, 
       # e.g., \code{Lqf[[j]][[i]]} is a matrix with rows corresponding to time points and
       # columns corresponding to quantile values on \code{qSup} for subject \eqn{i} and component \eqn{j}.
  qSup, # vector of support points on \eqn{[0,1]} of quantile functions in \code{Lqf}.
  ref_source = "pooled", # character specifying the principal geodesics of which sets of objects (here, distributions) are used to determine the reference point: \code{"start"} (using objects at the first time point from each component, i.e., \eqn{p} sets of objects), \code{"end"} (using objects at the last time point from each component, i.e., \eqn{p} sets of objects), \code{"both"} (using objects at the first and last time points from each component, i.e., \eqn{2p} sets of objects), \code{"pooled"} (default, using objects at the first and last time points pooling all components together, i.e., 2 sets of objects
  nknots, # number of knots used in pairwise synchronization
  lambdaH, # vector of the same length as \code{Lqf} holding the regularization parameters used in the estimation of subject-level deformation functions \eqn{H_i} from each component.
  lambdaG, # scalar holding the regularization parameter used in the estimation of \eqn{W_i}.
  lambdaPsi, # scalar holding the regularization parameter used in the estimation of \eqn{\Psi_j}.
  optns # list of optional control parameters for computing density functions: \code{bwDen}, \code{ndSup}, \code{dSup}, \code{delta}, \code{kernelDen}, \code{infSupport}, and \code{denLowerThreshold}. See \code{\link[frechet]{DenFMean}()} for details.
) {
  start_time <- Sys.time()
  
  ## obtain reference point omega_0 ----
  res_omega0 <- get_omega0( Lqf = Lqf, qSup = qSup, ref_source = ref_source )
  qf_omega0 <- res_omega0$qout
  
  ## compute amplitude factors A_{i,j} ----
  A <- get_amp( Lqf = Lqf, qf_omega0 = qf_omega0, qSup = qSup )
  
  ## compute normalized processes ----
  Lqf_normd <- normalize( Lqf = Lqf, qf_omega0 = qf_omega0, A = A )
  
  ## compute subject-level deformation functions ----
  res_H <- get_subjwf(
    obsGrid = obsGrid,
    Lqf_normd = Lqf_normd,
    qSup = qSup,
    nknots = nknots,
    lambda = lambdaH
  )
  
  ## compute (primitive estimates of) component tempo trajectories ----
  res_eta_pre <- align_comp(
    obsGrid = obsGrid,
    Lqf_normd = Lqf_normd,
    qSup = qSup,
    Hinv = res_H$Hinv,
    optns = optns
  )
  
  ## compute (primitive esimates of) global deformation functions ----
  res_G_pre <- align_glob(
    obsGrid = obsGrid,
    Lqf_normd = Lqf_normd,
    qSup = qSup,
    nknots = nknots,
    lambda = lambdaG
  )
  
  ## compute latent template ----
  res_tau <- get_template(
    obsGrid = obsGrid,
    Lqf_normd_aligned = res_G_pre$Lqf_normd_aligned,
    qSup = qSup,
    optns = optns
  )
  
  ## compute component-level deformation functions ----
  res_Psi <- get_compwf(
    obsGrid = obsGrid,
    qf_tau = res_tau$qout,
    Lqf_eta = res_eta_pre$qout,
    qSup = qSup,
    nknots = nknots
  )
  
  ## compute (resolved estimates of) component tempo trajectories ----
  res_eta <- get_comptempo(
    obsGrid = obsGrid,
    res_tau = res_tau,
    Psi = res_Psi$H
  )
  
  ## compute final estimates of global deformation functions ----
  res_G <- get_globwf(
    res_H = res_H, 
    res_Psi = res_Psi, 
    obsGrid = obsGrid
  )
  
  timing <- Sys.time() - start_time
  
  ## return ----
  list(
    res_omega0 = res_omega0, 
    A = A,
    Lqf_normd = Lqf_normd,
    res_H = res_H,
    res_eta_pre = res_eta_pre,
    res_G_pre = res_G_pre,
    res_tau = res_tau,
    res_Psi = res_Psi,
    res_eta = res_eta,
    res_G = res_G,
    obsGrid = obsGrid,
    qSup = qSup,
    timing = timing
  )
}

# get population-level cross-component alignment maps ----
get_pxcalign <- function (
    j1,j2, # indices or labels corresponding to which pairs of components the population-level cross-component alignment maps are to be computed.
    obsGrid, # vector of time points at which deformation functions are evaluated.
    Psi, # list of vectors for each component holding the estimates of component-level deformation functions on time points in \code{obsGrid}
    PsiInv # list of vectors for each component holding the inverses of estimated component-level deformation functions on time points in \code{obsGrid}
) {
  composite( x = obsGrid, f1 = PsiInv[[j1]], f2 = Psi[[j2]] )
  # output a vector of the population-level cross-component alignment maps aligning component \code{j1} to component \code{j2} evaluated on \code{obsGrid}.
}

# get subject-level cross-component alignment maps ----
get_sxcalign <- function (
    j1,j2, # indices or labels corresponding to which pairs of components the population-level cross-component alignment maps are to be computed.
    obsGrid, # vector of time points at which deformation functions are evaluated.
    G, # list of matrices for each component holding the estimates of global deformation functions, with rows corresponding to subjects and columns to time points in \code{obsGrid}, e.g., \code{G[[j]][i,]} holds the values of estimates of \eqn{G_{i,j}(t)} for \eqn{t} in \code{obsGrid}.
    Ginv # list of matrices for each component holding the estimated inverses of global deformation functions, with rows corresponding to subjects and columns to time points in \code{obsGrid}, e.g., \code{G[[j]][i,]} holds the values of estimates of \eqn{G^{-1}_{i,j}(t)} for \eqn{t} in \code{obsGrid}.
) {
  n_subj <- nrow(G[[1]])
  t(sapply( seq_len( n_subj ), function (i) {
    composite( x = obsGrid, f1 = Ginv[[j1]][i,], f2 = G[[j2]][i,] )
  }))
  # output a matrix of the subject-level cross-component alignment maps aligning component \code{j1} to component \code{j2}
  # with rows corresponding to subjects and columns corresponding to time points in \code{obsGrid}.
}
