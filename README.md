# Decision-grade evidence in rare diseases: R simulation templates

This repository accompanies the review manuscript:

> **Decision-grade evidence in rare diseases: A validity-first framework for design, small-sample analysis, and evidence borrowing**  
> Santiago Schnell  
> *Orphanet Journal of Rare Diseases* (under review)

The manuscript argues that rare-disease methods should be chosen along an auditable alignment chain:

**decision -> estimand / target trial -> data sources -> design -> analysis -> borrowing -> decision-grade conclusion**

The templates in this repository support that framework. They are intentionally lightweight, synthetic-data examples for calibrating operating characteristics and stress-testing assumptions before a protocol or statistical analysis plan (SAP) is finalized.

---

## What this code is - and is not

**This repository is:**
- A small set of executable **R templates** covering representative method families from the paper's three efficiency levers:
  1. increase information per participant,
  2. improve analytic efficiency,
  3. increase information beyond the trial.
- A **planning aid** for simulation-based operating-characteristic work.
- A **reproducibility support package** for the manuscript: all simulations use synthetic inputs, and the documentation explains the assumptions, expected inputs, and key outputs.

**This repository is not:**
- A validated clinical-trial software package.
- Production software or regulatory-grade code.
- A substitute for a protocol-specific simulation plan designed with an experienced statistician.

---

## Minimal reproducibility standard

For this repository to be useful in a protocol, SAP, or methods appendix, a reader should be able to identify:

1. **Which script was run** and for what decision problem.
2. **What synthetic inputs were assumed** (sample sizes, effect sizes, variability, heterogeneity, conflict, overlap, stopping rules, etc.).
3. **Which scenario grid was examined**, including at minimum a null scenario, a target-effect scenario, and one validity-stress scenario.
4. **What outputs were generated**, such as type I error, power, expected sample size, bias, precision / RMSE, interval calibration, and borrowing / ESS diagnostics when relevant.
5. **Which software environment was used**, captured via `scripts/session_info.R`.

These templates are designed to meet that standard when paired with a scenario grid, a short methods note, and archived session information.

---

## Repository structure

| Path | Purpose |
|---|---|
| `R/` | Core simulation templates for representative design, analysis, and borrowing problems |
| `scripts/run_demo.R` | Demonstration driver that loops over a scenario grid and prints operating-characteristic summaries |
| `scripts/session_info.R` | Records the R session for reproducibility |
| `inst/scenario_grid_example.csv` | Example scenario grid that can be replaced with disease-specific assumptions |
| `docs/mapping_to_paper.md` | Maps each template to the manuscript's sections, boxes, and figures |
| `docs/how_to_adapt.md` | Protocol-facing guide for adapting the templates to a disease-specific program |
| `CITATION.cff` | Citation metadata for the repository |

---

## Quick start

**Requirements:** R >= 4.2. The templates use base R only.

```r
# From the repository root:
source("scripts/run_demo.R")

# Capture the software environment used for the run:
source("scripts/session_info.R")
```

`run_demo.R` has two roles:
- **Sanity check:** execute each template once so the user can inspect the returned object structure.
- **Operating-characteristic loop:** evaluate the templates across a scenario grid and summarize key planning outputs.

Before using any numbers in a protocol, increase `N_SIM` beyond the default demo value and verify Monte Carlo precision is adequate for the decision at hand.

---

## Templates at a glance

| Template | Lever | Typical role in the paper | Key outputs |
|---|---|---|---|
| `sim_group_sequential.R` | I - Design efficiency | Interim monitoring / sample-size efficiency | probability of efficacy stop, probability of futility stop, expected N |
| `sim_mams_drop_a_loser.R` | I - Design efficiency | Multi-arm selection under scarcity | selected arm, dropped arm, actual versus maximum N |
| `sim_randomization_test.R` | II - Analytic efficiency | Design-valid small-sample inference | effect estimate, randomization p-value |
| `sim_map_borrowing_normal.R` | III - Borrowing | MAP / robust-MAP borrowing with prior-data conflict checks | posterior mean, posterior variance, posterior borrowing weight |
| `sim_external_controls_weighting.R` | III - Borrowing | External-control transportability with measured-confounder adjustment | naive and weighted estimates, SMDs, ESS fraction, weight diagnostics |

---

## How these templates should be used

Use the templates to answer questions like:

- Does a group sequential design reduce expected sample size enough to justify its complexity?
- Is a randomization-based analysis more stable than an asymptotic alternative at the achievable N?
- How sensitive is borrowing to heterogeneity or prior-data conflict?
- When external controls are considered, how much bias remains, and what does weighting do to overlap and ESS?

The correct output is not a single favorable number. The goal is a **stress-tested scenario set** that helps calibrate what evidentiary claim is defensible.

---

## Expected inputs and outputs

### Inputs
Each template assumes synthetic values for some combination of:
- sample sizes and allocation ratios,
- effect sizes and endpoint variability,
- heterogeneity or exchangeability parameters,
- alignment / overlap assumptions for external controls,
- stopping rules or adaptation rules where relevant.

### Outputs
Across templates, the main outputs are:
- type I error or false-positive rate,
- power or decision success rate,
- expected sample size,
- bias / precision or RMSE,
- posterior borrowing weight or posterior uncertainty,
- balance / overlap diagnostics and effective sample size for external controls.

These outputs are exactly the quantities the manuscript recommends reporting when methods are used to support high-stakes rare-disease decisions.

---

## Linking the code to the manuscript

The manuscript is a conceptual framework paper rather than a software paper. The code is therefore intentionally supportive and proportionate:
- the paper explains **when** a method is worth considering,
- the repository helps users examine **how it behaves** under realistic and adverse assumptions.

See:
- `docs/mapping_to_paper.md` for the manuscript mapping,
- `docs/how_to_adapt.md` for protocol-facing adaptation guidance.

---

## Release and citation notes

Use the metadata in `CITATION.cff` when citing the repository.

---

## License

MIT
