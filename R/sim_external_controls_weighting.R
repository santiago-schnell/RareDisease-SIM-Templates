# sim_external_controls_weighting.R
# Synthetic external-controls simulation with measured confounding and PS weighting
#
# Purpose:
#   Simulate a single-arm trial (treated arm only; no concurrent control) and
#   a large external control cohort subject to confounding by indication.
#   Demonstrate propensity-score (PS) weighting as a bias-reduction strategy,
#   and return the Box 4 minimum diagnostic set (SMD, ESS, weight quantiles)
#   alongside naive and weighted treatment effect estimates.
#
#   This is a planning / stress-test tool. In a real program the confounder
#   structure, outcome model, and external data source would all be pre-specified
#   in the statistical analysis plan.
#
# Data-generating model:
#
#   Single measured confounder X:
#     X_trt  ~ Normal(0,   1)   [trial-treated arm: centred distribution]
#     X_ext  ~ Normal(0.5, 1)   [external controls: shifted by 0.5 SD on average]
#     (0.5 SD shift = "moderate covariate imbalance" between sources)
#
#   Outcome model (linear, additive):
#     Y_trt = true_delta + beta * X_trt + eps,   eps ~ Normal(0, sigma^2)
#     Y_ext = 0          + beta * X_ext + eps,   eps ~ Normal(0, sigma^2)
#     beta = 0.6  (moderate confounding strength: a 1-SD difference in X
#                  translates to a 0.6-unit difference in Y, comparable to
#                  the outcome's SD; this is a realistic but detectable bias)
#
#   True estimand: E[Y_trt] - E[Y_{control if untreated}] = true_delta
#   Naive estimand: E[Y_trt] - E[Y_ext] ≠ true_delta because X_ext ≠ X_trt
#
# Bias decomposition (approximate):
#   Naive bias ≈ beta * (mean(X_ext) - mean(X_trt)) = 0.6 * 0.5 = -0.3
#   (external controls appear ~0.3 units WORSE on Y, making naive treatment
#    effect look smaller -- the classic confounding-by-indication pattern)
#
# Arguments:
#   n_trt      : sample size of the treated (trial) arm.
#   n_ext      : sample size of the external control cohort. Larger n_ext
#                reduces variance of weighted estimates but does not reduce bias
#                from unmeasured confounding or non-overlap.
#   true_delta : true treatment effect (mean difference). Bias = estimate - true_delta.
#   seed       : optional integer seed for reproducibility within a single call.
#
# Returns a named list (all Box 4 Item 6 diagnostics included):
#   naive_est           : unadjusted mean difference (biased by confounding).
#   weighted_est        : PS-weighted mean difference (unbiased for measured X).
#   true_delta          : true effect (for bias calculation in simulations).
#   x_trt_mean          : mean of X in treated arm.
#   x_ext_mean          : mean of X in external controls (raw, before weighting).
#   x_ext_weighted_mean : mean of X in external controls after PS weighting
#                         (should be close to x_trt_mean if weighting succeeds).
#   ess_ext             : effective sample size of the weighted external cohort
#                         (Kish formula: (sum w)^2 / sum(w^2)). Reduction from
#                         n_ext reflects weight variability / non-overlap.
#   ess_ext_fraction    : ess_ext / n_ext. Values below ~0.3 suggest severe
#                         positivity problems; consider trimming or restricting
#                         to the region of common support.
#   smd_unweighted      : standardized mean difference of X before weighting
#                         (|SMD| > 0.1 is often considered imbalanced).
#   smd_weighted        : standardized mean difference of X after weighting
#                         (target: |SMD| < 0.1, ideally < 0.05).
#   weight_quantiles    : quantiles of the PS weight distribution for external
#                         controls. Large maximum weights (> ~5) signal non-overlap
#                         and can inflate variance; trimming may be warranted.
#
# Watch-for: this template adjusts only for the ONE measured confounder X.
#   Unmeasured confounders produce residual bias that NO weighting can correct.
#   The E-value and quantitative bias analysis (cited in the paper) translate
#   residual bias into a "how strong would an unmeasured confounder need to be?"
#   statement. Non-overlap (positivity violations) is signalled by very low
#   ess_ext_fraction and extreme weight_quantiles -- in that regime the weighted
#   estimator is unstable regardless of sample size.
# Payoff: when X fully captures confounding and overlap is adequate, the
#   weighted estimator recovers the true delta; compare bias_naive vs bias_wtd
#   in the Monte Carlo OC output (run_demo.R).

simulate_external_control_weighting <- function(n_trt      = 30,
                                                n_ext      = 300,
                                                true_delta = 0.5,
                                                seed       = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # ── Simulate covariate X ────────────────────────────────────────────────────
  # Moderate covariate shift: external controls have X shifted +0.5 SD on average
  x_trt <- rnorm(n_trt, mean = 0,   sd = 1)
  x_ext <- rnorm(n_ext, mean = 0.5, sd = 1)

  # ── Simulate outcomes ────────────────────────────────────────────────────────
  # Confounding strength beta = 0.6: a 1-SD difference in X causes 0.6-unit
  # difference in Y -- comparable to sigma, so confounding is clearly detectable
  beta  <- 0.6
  sigma <- 1

  y_trt <- true_delta + beta * x_trt + rnorm(n_trt, 0, sigma)

  # External controls receive no treatment; their mean Y differs from trial
  # controls because X_ext ≠ X_trt (confounding by indication)
  y_ext <- 0 + beta * x_ext + rnorm(n_ext, 0, sigma)

  # ── Naive estimator (unadjusted; biased by ~beta * 0.5 = -0.3) ─────────────
  naive <- mean(y_trt) - mean(y_ext)

  # ── Propensity score estimation ─────────────────────────────────────────────
  # PS = P(being in treated arm | X): logistic model on the pooled sample
  # "1" = trial-treated; "0" = external control
  df     <- data.frame(group = c(rep(1, n_trt), rep(0, n_ext)),
                       x     = c(x_trt, x_ext))
  ps_fit <- glm(group ~ x, data = df, family = binomial())
  ps     <- predict(ps_fit, type = "response")

  # ── Inverse probability weights for external controls ───────────────────────
  # Odds weight: w = ps/(1-ps) reweights external controls to have the same
  # covariate distribution as the treated arm (ATT estimand: average treatment
  # effect in the treated population)
  ps_ext  <- ps[(n_trt + 1):(n_trt + n_ext)]
  w_ext   <- ps_ext / (1 - ps_ext)

  # Normalize weights to mean 1 for numerical stability (does not change the
  # weighted mean estimator but keeps weights on an interpretable scale)
  w_ext <- w_ext / mean(w_ext)

  # Weighted external control mean
  y_ext_w     <- weighted.mean(y_ext, w_ext)
  weighted_est <- mean(y_trt) - y_ext_w

  # ── Diagnostics (Box 4 Item 6 minimum set) ──────────────────────────────────

  # Kish effective sample size: (sum w)^2 / sum(w^2)
  # Measures equivalent unweighted sample size given weight variability;
  # substantial reduction from n_ext indicates non-overlap or extreme weights
  ess_ext  <- (sum(w_ext)^2) / sum(w_ext^2)

  # Weight distribution quantiles: flag extreme weights (>5 are often trimmed)
  w_quant  <- quantile(w_ext, probs = c(0, 0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1))

  # Standardized mean difference (SMD) of X, before and after weighting
  # Target after weighting: |SMD| < 0.1 (Austin 2009, Stat Med)
  smd_unweighted  <- (mean(x_trt) - mean(x_ext)) /
                     sqrt((var(x_trt) + var(x_ext)) / 2)

  # Weighted variance of X in external controls (for weighted SMD denominator)
  var_x_ext_w     <- sum(w_ext * (x_ext - weighted.mean(x_ext, w_ext))^2) /
                     sum(w_ext)
  smd_weighted    <- (mean(x_trt) - weighted.mean(x_ext, w_ext)) /
                     sqrt((var(x_trt) + var_x_ext_w) / 2)

  list(
    naive_est           = naive,
    weighted_est        = weighted_est,
    true_delta          = true_delta,
    x_trt_mean          = mean(x_trt),
    x_ext_mean          = mean(x_ext),
    x_ext_weighted_mean = weighted.mean(x_ext, w_ext),
    ess_ext             = ess_ext,
    ess_ext_fraction    = ess_ext / n_ext,
    smd_unweighted      = smd_unweighted,
    smd_weighted        = smd_weighted,
    weight_quantiles    = w_quant
  )
}
