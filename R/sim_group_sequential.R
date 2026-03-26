# sim_group_sequential.R
# Two-look group sequential template (stylized)
#
# Purpose:
#   Simulate a single replicate of a two-look group sequential trial and return
#   the stopping decision and test statistics. Use inside a Monte Carlo loop
#   (see scripts/run_demo.R) to estimate operating characteristics: P(early
#   efficacy stop), P(futility stop), and expected sample size E[N].
#
# Design assumptions (intentional simplifications for planning use):
#   - Continuous endpoint; two-arm trial with equal allocation at each look.
#   - Normal outcome model with KNOWN sigma (avoids t-distribution complexity).
#   - n1 participants enrolled before Look 1; n2 ADDITIONAL participants enrolled
#     before Look 2 (final). Total maximum N = n1 + n2.
#   - One-sided test for superiority (z statistics compared to positive thresholds).
#
# Arguments:
#   n1          : total participants at Look 1 (split equally: n1/2 per arm).
#   n2          : ADDITIONAL participants enrolled between Look 1 and Look 2
#                 (n2/2 per arm). Maximum total N = n1 + n2.
#   delta       : true treatment effect (mean difference, treatment - control).
#                 Set delta = 0 to estimate type I error; delta > 0 for power.
#   sigma       : common within-arm SD (assumed known; default 1 for standardised
#                 effect sizes).
#   efficacy_z  : z-statistic thresholds for declaring efficacy at [Look 1, Look 2].
#                 Defaults c(2.5, 2.0) are an approximate O'Brien-Fleming-style
#                 spending pattern (more conservative at early looks). These are
#                 NOT calibrated boundaries -- replace with gsDesign/rpact output
#                 for production use.
#   futility_z  : z-statistic thresholds for stopping for futility (non-binding)
#                 at [Look 1, Look 2]. Defaults c(-0.5, -0.2) correspond to
#                 mild futility rules (stop only if evidence strongly against
#                 treatment). Negative values mean "stop if z falls below this".
#   seed        : optional integer seed for reproducibility within a single call.
#                 For Monte Carlo loops, set the seed once before replicate().
#
# Returns a named list:
#   stop       : always TRUE (function always terminates; kept for consistency).
#   decision   : character string, one of:
#                  "efficacy@look1"      -- stopped early for efficacy
#                  "futility@look1"      -- stopped early for futility
#                  "efficacy@final"      -- efficacy declared at final look
#                  "futility@final"      -- futility declared at final look
#                  "continue_no_crossing"-- completed without crossing any boundary
#   z1         : z-statistic at Look 1.
#   z2         : z-statistic at Look 2 (NA if stopped at Look 1).
#
# Watch-for: Boundary thresholds must be calibrated by simulation under the
#   actual data-generating model (variance, missingness, non-normality). Nominal
#   alpha is not preserved if variance differs from the planning assumption.
# Payoff: E[N] < n1 + n2 when early stopping is frequent -- the efficiency gain
#   is the difference between maximum N and expected N under the true effect.

simulate_group_sequential <- function(n1,
                                      n2,
                                      delta,
                                      sigma = 1,
                                      efficacy_z = c(2.5, 2.0),
                                      futility_z = c(-0.5, -0.2),
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Helper: generate one arm's data (equal allocation assumed)
  gen_arm <- function(n, mean) rnorm(n, mean = mean, sd = sigma)

  # ── Look 1 ──────────────────────────────────────────────────────────────────
  y_t1 <- gen_arm(n1 / 2, delta)   # treatment arm at Look 1
  y_c1 <- gen_arm(n1 / 2, 0)       # control arm at Look 1

  # Two-sample z-statistic: (xbar_T - xbar_C) / SE
  # SE = sigma * sqrt(1/(n1/2) + 1/(n1/2)) = sigma * sqrt(4/n1)
  z1 <- (mean(y_t1) - mean(y_c1)) / (sigma * sqrt(4 / n1))

  if (z1 >= efficacy_z[1])
    return(list(stop = TRUE, decision = "efficacy@look1",  z1 = z1, z2 = NA))
  if (z1 <= futility_z[1])
    return(list(stop = TRUE, decision = "futility@look1",  z1 = z1, z2 = NA))

  # ── Look 2 (final) ──────────────────────────────────────────────────────────
  # Cumulative data: Look-1 participants PLUS n2/2 new participants per arm
  y_t2 <- c(y_t1, gen_arm(n2 / 2, delta))
  y_c2 <- c(y_c1, gen_arm(n2 / 2, 0))

  # SE at final look uses total cumulative N = n1 + n2
  z2 <- (mean(y_t2) - mean(y_c2)) / (sigma * sqrt(4 / (n1 + n2)))

  if (z2 >= efficacy_z[2])
    return(list(stop = TRUE, decision = "efficacy@final",  z1 = z1, z2 = z2))
  if (z2 <= futility_z[2])
    return(list(stop = TRUE, decision = "futility@final",  z1 = z1, z2 = z2))

  # Neither boundary crossed at either look
  list(stop = TRUE, decision = "continue_no_crossing", z1 = z1, z2 = z2)
}
