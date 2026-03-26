# How the code maps to the paper

**Manuscript:** *Decision-grade evidence under scarcity in rare diseases: a decision-centric, validity-first framework for design, analysis, and borrowing*  
Santiago Schnell. *Orphanet Journal of Rare Diseases* (under review).

The paper's central chain is:

**decision -> estimand / target trial -> data sources -> design -> analysis -> borrowing -> decision-grade conclusion**

The repository does not implement the full manuscript. Instead, it provides representative templates that let readers simulate operating characteristics and stress-test assumptions for several high-value method families.

---

## How to read this mapping

For each template below, the mapping identifies:
1. the main manuscript section it supports,
2. the specific box / figure / practical theme it operationalizes,
3. the operating-characteristic outputs that connect the code to the paper's narrative.

---

## Block I: increase information per participant (design efficiency)

### `sim_group_sequential.R`

| Paper element | Connection |
|---|---|
| Figure 1, Lever I | Interim monitoring as a way to gain information efficiency at fixed recruitment ceilings |
| Design section | Two-look group sequential trial as a simple prototype for adaptive monitoring |
| Box 2 | Makes the decision target and operating characteristics explicit |
| Box 6 | Supports the minimum simulation set: null calibration, target effect, expected N, and threshold sensitivity |
| Main watch-for | Boundary choices can misbehave if variance, event rate, or dropout differ from planning assumptions |
| Main payoff | Expected sample size can fall below maximum N when early stopping is plausible |

### `sim_mams_drop_a_loser.R`

| Paper element | Connection |
|---|---|
| Figure 1, Lever I | Multi-arm designs can screen multiple candidates without separate independent trials |
| Design section | MAMS logic for down-selection under scarcity |
| Box 5 exemplar logic | Reflects the same operational idea as rare-disease multi-arm screening programs |
| Box 6 | Allows users to simulate selection probability, expected N, and early dropping under favorable and adverse scenarios |
| Main watch-for | Premature dropping of a truly effective arm; complexity in shared-control settings |
| Main payoff | `n_final_actual` quantifies how much sample-size saving is achieved relative to `n_max` |

---

## Block II: improve analytic efficiency

### `sim_randomization_test.R`

| Paper element | Connection |
|---|---|
| Figure 1, Lever II | Design-valid inference when asymptotic approximations are unstable at very small N |
| Analysis section | Randomization-based inference as a small-sample option that remains tied to the allocation mechanism |
| Box 2 | Separates effect estimation from thresholding; the output includes an effect estimate plus a p-value |
| Box 6 | Supports null calibration and small-sample stress testing |
| Main watch-for | The permutation scheme must match the actual randomization mechanism |
| Main payoff | Empirical type I error and power can be checked without relying on large-sample theory |

---

## Block III: increase information beyond the trial

### `sim_map_borrowing_normal.R`

| Paper element | Connection |
|---|---|
| Figure 1, Lever III | Bayesian borrowing via MAP / robust-MAP priors |
| Borrowing section | Stylized borrowing template for a current and historical control mean |
| Box 2 | Makes prior assumptions, conflict scenarios, and claim calibration explicit |
| Box 6 | Supports heterogeneity and prior-data-conflict stress tests |
| Main watch-for | Overconfident borrowing when `tau` is too small or conflict is not examined |
| Main payoff | When sources are compatible, posterior uncertainty decreases; when conflict appears, robustification should reduce borrowing |

### `sim_external_controls_weighting.R`

| Paper element | Connection |
|---|---|
| Figure 3 | External controls are treated as a transportability workflow, not a default efficiency shortcut |
| Borrowing / external-controls section | Propensity-score weighting to align an external cohort with the treated group on measured covariates |
| Box 4 | Returns the minimum diagnostics emphasized in the manuscript: balance, weight distribution, and ESS |
| Box 6 | Supports non-overlap, residual-confounding, and ESS stress tests |
| Main watch-for | Good-looking weighted estimates are not enough if overlap collapses or key confounders are unavailable |
| Main payoff | Weighting can reduce measured-confounder bias when transportability is plausible, while ESS diagnostics show the cost of that adjustment |

---

## Scenario grid (`inst/scenario_grid_example.csv`)

The example scenario grid is intentionally simple. It is there to demonstrate the manuscript's logic that simulations should not be run only under ideal assumptions.

A minimum manuscript-consistent grid should include:
- a **null** scenario for calibration,
- a **target-effect** scenario matched to the intended decision,
- at least one **validity-stress** scenario such as heterogeneity, prior-data conflict, or non-overlap.

Users should adapt the grid to their endpoint, recruitment ceiling, and external-data risks.

---

## What is intentionally not implemented

To keep the repository lightweight and dependency-free, several methods discussed in the paper are not directly coded here, including:
- crossover / aggregated N-of-1 designs,
- longitudinal mixed-effects models,
- commensurate and power priors beyond the stylized MAP template,
- network meta-analysis,
- bias-reduced logistic regression,
- basket-trial borrowing structures,
- response-adaptive randomization.

Those omissions are intentional. The repository is a support package for the paper's framework, not an attempt to implement every method reviewed in the manuscript.
