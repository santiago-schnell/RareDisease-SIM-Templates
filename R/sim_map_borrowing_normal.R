# sim_map_borrowing_normal.R
# Stylized MAP / robust-MAP borrowing template for Normal outcomes
#
# Purpose:
#   Compute the posterior mean and variance of the current-trial control mean
#   (mu_c) under a two-component robust mixture prior: an informative MAP
#   component derived from historical data, plus a vague (near-flat) component
#   that protects against prior-data conflict (non-exchangeability).
#   This is a planning / stress-test tool, not a production analysis engine.
#
# Statistical model (simplified conjugate Normal-Normal):
#
#   Likelihood:    y_c | mu_c  ~  Normal(mu_c,  sigma^2 / n_c)
#                  [observed current control data as a sufficient summary]
#
#   MAP prior:     mu_c        ~  Normal(mu_h,  sigma^2/n_h + tau^2)
#                  [historical mean y_h shrunk by between-study heterogeneity tau]
#
#   Vague prior:   mu_c        ~  Normal(0,  1e6)
#                  [nearly flat; captures prior-data conflict and acts as a
#                   "safety valve" against over-borrowing]
#
#   Robust mixture (prior weight w on MAP, 1-w on vague):
#     pi(mu_c) = w * Normal(mu_h, sigma^2/n_h + tau^2)
#              + (1-w) * Normal(0, 1e6)
#
#   The posterior mixture weight w_post is updated via Bayes' theorem using
#   the marginal likelihood of y_c under each component (Normal-Normal
#   conjugacy: marginal is Normal(m0, v0 + v_like)).
#
# Arguments:
#   y_h   : observed historical control mean (sufficient statistic).
#   n_h   : historical control sample size (governs precision of MAP component).
#   y_c   : observed current control mean (sufficient statistic).
#   n_c   : current control sample size.
#   sigma : common within-arm SD, assumed known and shared across sources.
#           (Planning assumption; sensitivity to this is worth exploring.)
#   tau   : between-study heterogeneity SD. Governs how much the MAP prior
#           discounts the historical mean relative to the current data.
#           tau = 0   => no discounting (treat sources as identical).
#           tau = 0.5 => moderate heterogeneity (typical planning assumption).
#           tau >> sigma => MAP prior becomes very diffuse; borrowing is minimal.
#   w     : prior mixture weight on the MAP (informative) component, 0 < w < 1.
#           w = 0.8 (default) means 80% prior probability that sources are
#           exchangeable. Reduced to w_post after seeing y_c if conflict is
#           detected.
#
# Returns a named list:
#   mean    : posterior mean of mu_c under the robust mixture posterior.
#   var     : posterior variance of mu_c (law of total variance across components:
#             Var = E[Var|component] + Var[E|component]).
#   w_post  : posterior weight on the MAP (informative) component.
#             w_post near 1 => data are compatible with the historical source
#               (borrowing is active; posterior mean pulled toward y_h).
#             w_post near 0 => severe prior-data conflict detected; posterior
#               relies almost entirely on current data via the vague component.
#             Monitoring w_post across a scenario grid is the primary diagnostic
#             for whether borrowing is warranted in a given analysis.
#
# Watch-for: prior-data conflict (y_c very different from y_h) can still inflate
#   false positives if w is set too high and tau is too small. Design-stage
#   simulations should include explicit conflict scenarios (see run_demo.R).
# Payoff: when sources are exchangeable (y_c close to y_h), the MAP component
#   dominates, reducing posterior variance relative to using only current data
#   -- quantified by comparing E[sd_post] in compatible vs conflict scenarios.

map_posterior_mu_c <- function(y_h, n_h,
                               y_c, n_c,
                               sigma = 1,
                               tau   = 0.5,
                               w     = 0.8) {

  # ── Prior parameters ────────────────────────────────────────────────────────

  # MAP component: Normal(y_h, sigma^2/n_h + tau^2)
  # sigma^2/n_h = sampling variance of historical mean; tau^2 = between-study var
  v_map <- sigma^2 / n_h + tau^2
  m_map <- y_h

  # Vague component: nearly flat Normal centred at 0 with very large variance
  v_vague <- 1e6
  m_vague <- 0

  # ── Likelihood variance ─────────────────────────────────────────────────────
  v_like <- sigma^2 / n_c   # sampling variance of current control mean

  # ── Posterior under each component (Normal-Normal conjugacy) ────────────────
  # For prior Normal(m0, v0) and likelihood Normal(mu_c, v_like):
  #   posterior precision  = 1/v0 + 1/v_like
  #   posterior mean       = v_post * (m0/v0 + y_c/v_like)
  post_norm <- function(m0, v0) {
    v_post <- 1 / (1 / v0 + 1 / v_like)
    m_post <- v_post * (m0 / v0 + y_c / v_like)
    list(mean = m_post, var = v_post)
  }

  post_map   <- post_norm(m_map,   v_map)
  post_vague <- post_norm(m_vague, v_vague)

  # ── Posterior mixture weights (Bayes factor update) ──────────────────────────
  # Marginal likelihood of y_c under each prior component:
  #   p(y_c | component) = Normal(y_c; m0, v0 + v_like)   [predictive distribution]
  mlik_map   <- dnorm(y_c, mean = m_map,   sd = sqrt(v_map   + v_like))
  mlik_vague <- dnorm(y_c, mean = m_vague, sd = sqrt(v_vague + v_like))

  # Updated weight on MAP component proportional to prior weight * marginal lik
  w_post <- (w * mlik_map) / (w * mlik_map + (1 - w) * mlik_vague)

  # ── Mixture posterior moments ────────────────────────────────────────────────
  # Posterior mean (weighted average of component means)
  mean_mix <- w_post * post_map$mean + (1 - w_post) * post_vague$mean

  # Posterior variance via law of total variance:
  #   Var[mu_c] = E[Var[mu_c | component]] + Var[E[mu_c | component]]
  #             = sum_k w_k * v_k  +  sum_k w_k * (m_k - mean_mix)^2
  # Equivalently: E[mu_c^2] - (E[mu_c])^2, computed as:
  var_mix <- w_post       * (post_map$var   + post_map$mean^2)   +
             (1 - w_post) * (post_vague$var + post_vague$mean^2) -
             mean_mix^2

  list(mean   = mean_mix,
       var    = var_mix,
       w_post = w_post)
}
