# run_demo.R — Demo runner for rare-disease scarcity simulation templates
#
# Run from repository root: source("scripts/run_demo.R")
#
# Structure:
#   PART 1: Single-replicate sanity checks for each template
#   PART 2: Monte Carlo operating-characteristics (OC) loop over a scenario grid
#
# All code uses base R; no packages required.

source("R/sim_group_sequential.R")
source("R/sim_mams_drop_a_loser.R")
source("R/sim_randomization_test.R")
source("R/sim_map_borrowing_normal.R")
source("R/sim_external_controls_weighting.R")

# ─── PART 1: Single-replicate sanity checks ──────────────────────────────────

cat("\n=================================================================\n")
cat("PART 1: Single-replicate sanity checks\n")
cat("=================================================================\n")

set.seed(2026)

cat("\n--- Group sequential (single replicate) ---\n")
gs <- simulate_group_sequential(n1 = 20, n2 = 20, delta = 0.5, sigma = 1)
print(gs)

cat("\n--- MAMS drop-a-loser (single replicate) ---\n")
# Note: n_max is the ceiling if no arm is dropped; actual final N is smaller
mams <- simulate_mams_drop_a_loser(
  K = 3, n_interim = 40, n_max = 80, deltas = c(0, 0.2, 0.6)
)
print(mams)
cat("Efficiency gain: actual final N =", mams$n_final_actual,
    "vs n_max =", 80, "\n")

cat("\n--- Randomization test (single replicate) ---\n")
N   <- 20
trt <- rep(c(0, 1), each = N / 2)
y   <- rnorm(N, mean = 0.2 * trt, sd = 1)
rt  <- perm_test_diff_means(y = y, trt = trt, n_perm = 2000)
print(rt)

cat("\n--- MAP / robust-MAP (single replicate, prior-data conflict scenario) ---\n")
# y_c = 0.4 vs y_h = 0.1: detectable conflict -> posterior weight shifts toward vague
post_conflict <- map_posterior_mu_c(
  y_h = 0.1, n_h = 60, y_c = 0.4, n_c = 10, sigma = 1, tau = 0.4, w = 0.8
)
cat("Prior-data conflict scenario: w_post (borrowing weight) =",
    round(post_conflict$w_post, 3), "\n")
print(post_conflict)

cat("\n--- MAP / robust-MAP (single replicate, exchangeable scenario) ---\n")
post_compat <- map_posterior_mu_c(
  y_h = 0.1, n_h = 60, y_c = 0.12, n_c = 10, sigma = 1, tau = 0.4, w = 0.8
)
cat("Compatible scenario: w_post (borrowing weight) =",
    round(post_compat$w_post, 3), "\n")
print(post_compat)

cat("\n--- External controls weighting (single replicate) ---\n")
cat("  Setup: single-arm trial (no concurrent control) + external cohort\n")
cat("  with confounding by indication (X_ext shifted +0.5 SD vs X_trt)\n")
ec <- simulate_external_control_weighting(n_trt = 30, n_ext = 300, true_delta = 0.5)
cat("Naive estimate:", round(ec$naive_est, 3),
    "| Weighted:", round(ec$weighted_est, 3),
    "| True delta:", ec$true_delta, "\n")
cat("SMD unweighted:", round(ec$smd_unweighted, 3),
    "| SMD after weighting:", round(ec$smd_weighted, 3), "\n")
cat("ESS (fraction of n_ext):", round(ec$ess_ext_fraction, 3), "\n")

# ─── PART 2: Monte Carlo operating characteristics ───────────────────────────

cat("\n=================================================================\n")
cat("PART 2: Monte Carlo operating characteristics\n")
cat("(N_SIM replicates per scenario; may take a few seconds)\n")
cat("=================================================================\n")

N_SIM <- 500   # Number of Monte Carlo replicates per scenario.
               # 500  => fast preview (~seconds); SE of a proportion ~ 0.02.
               # 2000 => planning document quality; SE ~ 0.01.
               # 5000 => publication-quality; SE ~ 0.007 for adaptive designs.
               # For adaptive designs and Bayesian borrowing, use >= 2000.

# Load scenario grid
grid <- read.csv("inst/scenario_grid_example.csv", stringsAsFactors = FALSE)
cat("\nScenario grid:\n")
print(grid)

# ── 2a: Group sequential OC ──────────────────────────────────────────────────
cat("\n--- Group sequential OC (n1=20, n2=20, sigma=1) ---\n")
cat(sprintf("%-25s %8s %8s %8s\n", "scenario", "P(eff)", "P(fut)", "E[N]"))
for (i in seq_len(nrow(grid))) {
  sc <- grid[i, ]
  results <- replicate(N_SIM, {
    r <- simulate_group_sequential(
      n1 = 20, n2 = 20, delta = sc$delta, sigma = sc$sigma
    )
    # Return: stopped for efficacy (1/0), stopped for futility (1/0), N used
    eff  <- as.integer(grepl("efficacy", r$decision))
    fut  <- as.integer(grepl("futility", r$decision))
    n_used <- if (grepl("look1", r$decision)) 20 else 40
    c(eff, fut, n_used)
  })
  p_eff  <- mean(results[1, ])
  p_fut  <- mean(results[2, ])
  e_n    <- mean(results[3, ])
  cat(sprintf("%-25s %8.3f %8.3f %8.1f\n", sc$scenario, p_eff, p_fut, e_n))
}

# ── 2b: MAMS drop-a-loser OC ─────────────────────────────────────────────────
cat("\n--- MAMS drop-a-loser OC (K=3, n_interim=40, n_max=80, sigma=1) ---\n")
cat("  Effect profile: arm1=null(0), arm2=small(0.2), arm3=moderate(delta)\n")
cat(sprintf("%-25s %8s %8s %8s\n", "scenario", "P(sel3)", "E[N_act]", "N_max"))
for (i in seq_len(nrow(grid))) {
  sc <- grid[i, ]
  results <- replicate(N_SIM, {
    r <- simulate_mams_drop_a_loser(
      K = 3, n_interim = 40, n_max = 80,
      deltas = c(0, 0.2, sc$delta), sigma = sc$sigma
    )
    c(as.integer(r$selected_arm == 3), r$n_final_actual)
  })
  p_sel3  <- mean(results[1, ])
  e_n_act <- mean(results[2, ])
  cat(sprintf("%-25s %8.3f %8.1f %8d\n", sc$scenario, p_sel3, e_n_act, 80L))
}

# ── 2c: Randomization test OC ────────────────────────────────────────────────
cat("\n--- Randomization test OC (N=20, n_perm=999) ---\n")
cat(sprintf("%-25s %8s\n", "scenario", "P(p<0.05)"))
for (i in seq_len(nrow(grid))) {
  sc <- grid[i, ]
  p_vals <- replicate(N_SIM, {
    trt_i <- rep(c(0, 1), each = 10)
    y_i   <- rnorm(20, mean = sc$delta * trt_i, sd = sc$sigma)
    perm_test_diff_means(y = y_i, trt = trt_i, n_perm = 999)$p_value
  })
  cat(sprintf("%-25s %8.3f\n", sc$scenario, mean(p_vals < 0.05)))
}

# ── 2d: MAP borrowing OC: bias and ESS across scenarios ──────────────────────
cat("\n--- MAP / robust-MAP OC (n_h=60, n_c=10, sigma=1, tau=0.4, w=0.8) ---\n")
cat("  Measures shrinkage of posterior mean toward historical vs current data\n")
cat(sprintf("%-25s %8s %8s %8s\n", "scenario", "E[mu_post]", "E[w_post]", "E[sd_post]"))
y_h_ref <- 0.1
for (i in seq_len(nrow(grid))) {
  sc <- grid[i, ]
  results <- replicate(N_SIM, {
    y_c_i <- rnorm(1, sc$delta, sc$sigma / sqrt(10))
    post  <- map_posterior_mu_c(
      y_h = y_h_ref, n_h = 60, y_c = y_c_i, n_c = 10,
      sigma = sc$sigma, tau = 0.4, w = 0.8
    )
    c(post$mean, post$w_post, sqrt(post$var))
  })
  cat(sprintf("%-25s %8.3f %8.3f %8.3f\n",
              sc$scenario,
              mean(results[1, ]),
              mean(results[2, ]),
              mean(results[3, ])))
}

# ── 2e: External controls weighting OC ───────────────────────────────────────
cat("\n--- External controls weighting OC (n_trt=30, n_ext=300) ---\n")
cat(sprintf("%-25s %8s %8s %8s %8s\n",
            "scenario", "bias_naive", "bias_wtd", "E[SMD_w]", "E[ESS_frac]"))
for (i in seq_len(nrow(grid))) {
  sc <- grid[i, ]
  results <- replicate(N_SIM, {
    r <- simulate_external_control_weighting(
      n_trt = 30, n_ext = 300, true_delta = sc$delta
    )
    c(r$naive_est - r$true_delta,
      r$weighted_est - r$true_delta,
      r$smd_weighted,
      r$ess_ext_fraction)
  })
  cat(sprintf("%-25s %8.3f %8.3f %8.3f %8.3f\n",
              sc$scenario,
              mean(results[1, ]),
              mean(results[2, ]),
              mean(results[3, ]),
              mean(results[4, ])))
}

cat("\n=================================================================\n")
cat("Done. Increase N_SIM at the top of this script for smoother estimates.\n")
cat("=================================================================\n")
