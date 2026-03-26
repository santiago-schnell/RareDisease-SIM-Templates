# MAMS "drop-a-loser" template (stylized)
#
# Structure:
#   - K experimental arms + shared control
#   - Interim: drop the worst-performing experimental arm (by observed mean diff)
#   - Final: compare best remaining arm to control (no multiplicity adjustment here)
#
# Arguments:
#   K          : number of experimental arms
#   n_interim  : TOTAL sample size at interim (across all K experimental + 1 control
#                arms combined; must be divisible by K + 1)
#   n_max      : MAXIMUM total sample size if NO arms were dropped (must be > n_interim
#                and divisible by K + 1). After dropping one arm the actual
#                final N will be n_interim + K * ((n_max - n_interim) / (K + 1)),
#                which is LESS than n_max -- that is the efficiency gain.
#   deltas     : true treatment effects for each experimental arm (length K)
#   sigma      : common SD
#
# Watch-for: multiplicity + time trends in shared controls + premature dropping.
# Payoff: fewer wasted patients and faster down-selection of candidates.
#         Actual final N = n_interim + K * n_add_per_arm < n_max (by design).

simulate_mams_drop_a_loser <- function(K = 3,
                                       n_interim = 30,
                                       n_max = 60,
                                       deltas = rep(0, K),
                                       sigma = 1,
                                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  stopifnot(length(deltas) == K)
  # n_interim and n_max are TOTAL sample sizes across all K + 1 arms.
  # n_max is the ceiling if no arm had been dropped; after dropping one arm,
  # the realised final N is smaller (that is the efficiency gain).
  stopifnot(n_interim %% (K + 1) == 0)
  stopifnot(n_max > n_interim)
  stopifnot((n_max - n_interim) %% (K + 1) == 0)

  n_per_arm_i   <- n_interim / (K + 1)
  n_add_per_arm <- (n_max - n_interim) / (K + 1)
  n_final_actual <- n_interim + K * n_add_per_arm  # one arm dropped

  # Generate interim data
  y_c_i <- rnorm(n_per_arm_i, 0, sigma)
  y_e_i <- lapply(seq_len(K), function(k) rnorm(n_per_arm_i, deltas[k], sigma))

  diffs_i <- sapply(y_e_i, function(y) mean(y) - mean(y_c_i))
  drop_k  <- which.min(diffs_i)
  keep    <- setdiff(seq_len(K), drop_k)

  # Generate additional data for kept arms and control only
  y_c_f_add <- rnorm(n_add_per_arm, 0, sigma)
  y_e_f_add <- lapply(seq_len(K), function(k) {
    if (k %in% keep) rnorm(n_add_per_arm, deltas[k], sigma) else NULL
  })

  y_c <- c(y_c_i, y_c_f_add)
  y_e <- lapply(seq_len(K), function(k) {
    if (k %in% keep) c(y_e_i[[k]], y_e_f_add[[k]]) else NULL
  })

  diffs_f <- rep(NA_real_, K)
  for (k in keep) diffs_f[k] <- mean(y_e[[k]]) - mean(y_c)

  best_k <- keep[which.max(diffs_f[keep])]

  list(
    dropped_arm    = drop_k,
    kept_arms      = keep,
    interim_diffs  = diffs_i,
    final_diffs    = diffs_f,
    selected_arm   = best_k,
    n_final_actual = n_final_actual   # actual N used (< n_max due to drop)
  )
}
