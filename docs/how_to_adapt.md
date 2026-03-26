# How to adapt the templates to your own disease and design

This guide is for readers who want to move from the paper's conceptual framework to a protocol-facing simulation plan. The organizing logic is the same as in the manuscript:

**decision -> estimand / target trial -> data sources -> design -> analysis -> borrowing**

The templates are intentionally simplified. Your task is not to make the code look sophisticated; it is to make the assumptions explicit, stress-test them, and report outputs that are fit for the intended decision.

---

## Before you change any code

Write down the answers to these questions first.

### 1. What decision is being calibrated?
Examples:
- exploratory signal detection,
- dose or regimen down-selection,
- supportive contextualization,
- confirmatory evidence generation,
- post-marketing follow-up or safety augmentation.

The intended decision determines which outputs matter most. Confirmatory uses place more weight on false-positive control; supportive uses may care more about bias, overlap, and uncertainty calibration.

### 2. What is the estimand?
At minimum, specify the five ICH E9(R1)-style attributes:
- population,
- variable / endpoint,
- handling of intercurrent events,
- population-level summary measure,
- sensitivity strategy linked to the estimand.

Do not skip this step. In rare-disease work, changing the estimand often changes which design, analysis, or borrowing strategy is defensible.

### 3. What data sources are actually available?
List the concurrent trial data, historical controls, registry or natural-history data, and any mechanistic or expert inputs that affect planning assumptions.

### 4. What sample size is realistically achievable?
Rare-disease recruitment ceilings are real. Simulate the program you can actually run, not the program you wish were feasible.

---

## Minimum simulation reporting standard

A protocol-facing run should document, at minimum:

1. **Decision target and estimand**  
   What decision is being supported, what estimand is targeted, and what evidence role the external data will play, if any.

2. **Base-case data-generating mechanism**  
   Recruitment, allocation, outcome model, follow-up, intercurrent events, missingness, and analysis model.

3. **Reference scenarios**  
   Include at least:
   - a null / no-benefit scenario,
   - a target-effect scenario,
   - a near-threshold or MCID scenario when relevant.

4. **Validity-stress scenarios**  
   Add method-relevant failures such as heterogeneity, model misspecification, informative missingness, non-overlap, calendar drift, residual confounding, or prior-data conflict.

5. **Minimum outputs**  
   Report the outputs that matter for your method family: type I error, power / decision success, expected sample size, bias, precision / RMSE, interval calibration, and borrowing / ESS diagnostics when relevant.

6. **Reproducibility details**  
   Record the scenario grid, number of Monte Carlo replicates, random-seed strategy, and session information.

This repository is designed so those items can be recovered from a simulation package built around the included templates.

---

## Template 1: Group sequential design (`sim_group_sequential.R`)

### Parameters to replace

| Parameter | Synthetic default | Replace with |
|---|---|---|
| `delta` | 0.5 | Your clinically meaningful treatment effect in endpoint units (or in standardized units if you keep `sigma = 1`) |
| `sigma` | 1 | Your endpoint variability from natural history or prior trials |
| `n1` | 20 | Participants available at the interim look |
| `n2` | 20 | Additional accrual planned after interim |
| `efficacy_z` | `c(2.5, 2.0)` | Boundaries calibrated in appropriate design software; these defaults are illustrative only |
| `futility_z` | `c(-0.5, -0.2)` | Non-binding or binding futility boundaries justified in the protocol |

### What to check
- false-positive control under the null,
- power or decision success under the target-effect scenario,
- expected sample size relative to the maximum N,
- sensitivity to variance inflation or event-rate misspecification.

### Common mistake
Do not lift the default boundaries into a real protocol. They are demonstrative, not protocol-ready.

---

## Template 2: MAMS drop-a-loser (`sim_mams_drop_a_loser.R`)

### Parameters to replace

| Parameter | Synthetic default | Replace with |
|---|---|---|
| `K` | 3 | Number of experimental arms |
| `n_interim` | 40 | Total accrued participants before down-selection |
| `n_max` | 80 | Maximum total N if no arm were dropped |
| `deltas` | `c(0, 0.2, 0.5)` | True effects for each arm across plausible scenarios |
| `sigma` | 1 | Outcome variability |

### What to check
- correct-arm selection probability,
- sample-size savings (`n_max - n_final_actual`),
- robustness when one promising arm is only modestly better than competitors,
- behavior under time trends or shared-control instability, if those are plausible.

### Scope note
This template implements a simple drop-a-loser structure. More complex multi-stage or platform settings will need additional bookkeeping and, often, a richer scenario grid.

---

## Template 3: Randomization test (`sim_randomization_test.R`)

### Parameters to replace

| Parameter | Synthetic default | Replace with |
|---|---|---|
| `y` | synthetic Normal outcome | Simulated outcomes consistent with your endpoint and estimand |
| `trt` | balanced allocation | Your actual allocation pattern |
| `n_perm` | 10000 | A value high enough to stabilize the randomization p-value for the use case |

### What to check
- empirical false-positive rate under the null,
- decision stability at the achievable N,
- agreement or disagreement with asymptotic alternatives,
- whether the permutation scheme matches the actual randomization scheme.

### Important constraint
If your study uses stratified or covariate-adaptive randomization, the permutation mechanism must respect that design. A simple unrestricted permutation is not automatically valid.

---

## Template 4: MAP / robust-MAP borrowing (`sim_map_borrowing_normal.R`)

### Parameters to replace

| Parameter | Synthetic default | Replace with |
|---|---|---|
| `y_h` | 0.1 | Historical control summary |
| `n_h` | 60 | Historical sample size |
| `y_c` | scenario-specific | Current-control summary under the simulated scenario |
| `n_c` | 10 | Concurrent control size, if any |
| `sigma` | 1 | Outcome SD |
| `tau` | 0.4 to 0.5 | Heterogeneity / exchangeability assumption |
| `w` | 0.8 | Prior probability that borrowing is appropriate |

### What to check
- posterior borrowing weight under compatible and conflicting scenarios,
- sensitivity to low / medium / high `tau`,
- whether borrowing remains conservative when prior-data conflict is introduced,
- how claim strength changes if borrowing is downgraded.

### Key advice
Never present a single `tau` run as if it were enough. The defensible question is how conclusions behave across a justified range of heterogeneity assumptions.

---

## Template 5: External controls with PS weighting (`sim_external_controls_weighting.R`)

### Parameters to replace

| Parameter | Synthetic default | Replace with |
|---|---|---|
| `n_trt` | 30 | Treated-arm size in the current study |
| `n_ext` | 300 | External-control cohort size |
| `true_delta` | 0.5 | Treatment effect used for simulation benchmarking |
| covariate shift | 0.5 SD | Plausible baseline imbalance between current and external data |
| confounding strength | 0.6 | Relationship between the measured confounder and the outcome |
| `sigma` | 1 | Outcome residual SD |

### What to check
- unweighted versus weighted bias,
- balance before and after weighting,
- effective sample size after weighting,
- whether non-overlap or extreme weights would force a downgrade from augmentation to contextualization.

### Extension note
The toy template uses a single confounder. Real applications usually need multiple covariates and a more explicit target-trial alignment exercise.

---

## Building a scenario grid that matches the paper

A minimum rare-disease scenario grid should include:

| Scenario purpose | What it helps calibrate |
|---|---|
| Null / no-benefit | false-positive control |
| Near-threshold / MCID | decision sensitivity at the clinical margin |
| Target effect | primary operating characteristics |
| Heterogeneity / model misspecification | robustness of the planned analysis |
| Prior-data conflict or non-overlap | robustness of borrowing and external-control claims |

You can add optimistic or pessimistic scenarios, but do not omit the adverse ones.

---

## Checklist before using results in a protocol or SAP

- [ ] Decision and estimand are written down before simulation starts.
- [ ] Scenario grid includes at least one validity-stress scenario.
- [ ] `N_SIM` is large enough for the intended use.
- [ ] The outputs reported match the evidentiary role being claimed.
- [ ] Session information has been captured via `scripts/session_info.R`.
- [ ] Any external-data use includes overlap / ESS diagnostics and a downgrade rule if assumptions fail.
- [ ] Protocol or SAP text states which scenarios were examined and how conclusions were calibrated from them.

---

## Getting support

This draft package may be shared before a public repository URL exists. Use the included documentation locally, and add the issue tracker or public contact route when the repository is released.
