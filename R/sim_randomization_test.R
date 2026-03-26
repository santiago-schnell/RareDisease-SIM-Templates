# sim_randomization_test.R
# Randomization (permutation) test template for a two-arm experiment
#
# Purpose:
#   Compute a two-sided randomization-based p-value for the observed difference
#   in means. Unlike asymptotic z- or t-tests, this approach derives its null
#   distribution directly from the randomization scheme: it permutes treatment
#   labels the same way the original allocation did, making it valid even at
#   very small N without relying on large-sample approximations.
#
# Statistical basis:
#   The null hypothesis is that treatment assignment has no effect on outcomes
#   (the sharp null: every individual's outcome is the same regardless of arm).
#   Under this null, any relabelling of {0,1} assignments to the same N units
#   is equally likely. The p-value is the proportion of permutations whose
#   absolute mean difference is at least as extreme as the observed value.
#
# Arguments:
#   y      : numeric vector of outcomes (length N).
#   trt    : integer/numeric vector of treatment indicators, same length as y;
#            must contain only 0 (control) and 1 (treatment).
#   n_perm : number of random permutations used to approximate the null
#            distribution. Default 10000 gives stable p-values to ~2 decimal
#            places. For very small N (< ~15 total), consider complete
#            enumeration instead (choose(N, N/2) permutations).
#   seed   : optional integer seed for reproducibility within a single call.
#
# Returns a named list:
#   estimate : observed difference in means (treatment mean - control mean).
#              This is an effect SIZE estimate, not a test statistic -- report
#              it alongside the p-value; the p-value alone conveys no information
#              about magnitude or clinical relevance.
#   p_value  : two-sided randomization p-value. Interpretation: the probability,
#              under the sharp null, of observing a mean difference at least as
#              large in absolute value as the one actually observed.
#              NOTE: this p-value is NOT a substitute for an effect estimate
#              with confidence interval -- always report `estimate` as well.
#
# Watch-for: the randomization-based p-value is valid only under the same
#   allocation mechanism used in the actual trial. Stratified or covariate-
#   adaptive randomization requires a matching stratified permutation scheme
#   (not implemented here). Mismatched schemes can inflate type I error.
# Payoff: exact type I error control at small N with no normality assumption;
#   robust to non-standard outcome distributions.

perm_test_diff_means <- function(y, trt, n_perm = 10000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(y) == length(trt))
  stopifnot(all(trt %in% c(0, 1)))

  # Observed test statistic: difference in group means
  obs <- mean(y[trt == 1]) - mean(y[trt == 0])

  # Build the permutation null distribution by randomly reassigning treatment
  # labels while keeping group sizes fixed (sampling WITHOUT replacement
  # preserves the same number of treated and control units in each permutation)
  diffs <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    trt_b    <- sample(trt, replace = FALSE)   # permute labels, keep sizes
    diffs[b] <- mean(y[trt_b == 1]) - mean(y[trt_b == 0])
  }

  # Two-sided p-value: proportion of permutations as or more extreme than observed
  p_two_sided <- mean(abs(diffs) >= abs(obs))

  list(estimate = obs, p_value = p_two_sided)
}
