# Dynamic Risk Surface Modelling of Cytopenia Recovery and Infection Vulnerability During Bispecific Antibody Therapy in High-Risk Lymphoma

## Overview

This repository contains the computational model, input assumptions, model-generated results, figures, and manuscript associated with the study:

**“Dynamic Risk Surface Modelling of Cytopenia Recovery and Infection Vulnerability During Bispecific Antibody Therapy in High-Risk Lymphoma.”**

The study is a **literature-parameterized computational proof-of-concept** designed to examine the dynamic relationship between cytopenia recovery and broader infection vulnerability during bispecific-antibody therapy.

No patient-level clinical data were used for model fitting, parameter estimation, or validation. All numerical results in this repository are outputs of the specified mathematical model and computational experiments and should not be interpreted as observed clinical incidences, validated risk estimates, treatment effects, or clinical decision thresholds.

---

## Study Objective

The study examines whether a dynamic reserve–stress–hazard system can represent and interrogate time-dependent separation between cytopenia depth and broader infection vulnerability, while identifying structural assumptions and parameters that materially influence this behaviour.

The framework is intended for hypothesis generation, mathematical interrogation, reproducibility, and future observed-data calibration rather than bedside prediction or treatment recommendation.

---

## Mathematical Framework

The model contains seven continuous state variables representing:

- `R(t)` — marrow reserve
- `N(t)` — neutrophil reserve
- `P(t)` — platelet reserve
- `G(t)` — immune competence
- `C(t)` — inflammatory burden
- `B(t)` — microbial pressure
- `H(t)` — cumulative infection hazard

These states are coupled through a nonlinear ordinary differential equation (ODE) system.

Cumulative model-generated infection risk is obtained from the cumulative hazard as:

`P_inf(t) = 1 - exp[-H(t)]`

A computational stability coordinate combines reserve preservation with inflammatory and microbial burden. A provisional three-state latent-regime representation summarizes the trajectory as **controlled**, **strained**, or **overloaded**.

The latent-state probabilities are numerical model summaries and are not empirically estimated clinical-state prevalences.

---

## Treatment-Stress Forcing

Treatment-associated effects are represented using a normalized, dimensionless structural forcing function `E(t)`.

`E(t)` is **not** intended to represent:

- serum drug concentration,
- receptor occupancy,
- a validated pharmacokinetic profile,
- or a drug-specific clinical dosing regimen.

Its amplitudes, timing, and decay coefficient are structural modelling assumptions.

Alternative weekly-after-step-up and q21-after-step-up schedules are included only as structural sensitivity experiments to evaluate dependence of model behaviour on the assumed forcing schedule.

They should not be interpreted as comparative clinical-regimen analyses.

---

## Computational Analyses

The executable workflow includes:

1. 84-day reference simulation
2. Prespecified scenario perturbations
3. Frailty–dose-intensity dynamic risk surface
4. Global sensitivity analysis using 600 Latin-hypercube parameter sets and PRCC
5. Local elasticity analysis
6. 500-run stochastic perturbation analysis
7. Latent-regime trajectory analysis
8. Synthetic AIC/BIC mechanism-recovery experiment
9. Numerical robustness assessment
10. Structural exposure-schedule sensitivity analysis

The synthetic mechanism-recovery experiment uses pseudo-observations generated from the mechanistic reference trajectory and therefore represents an **internal structural recovery check**, not external or clinical validation.

---

## Selected Model-Generated Results

Under the declared reference parameterization:

- Neutrophil-reserve nadir: **0.724** at day **74.5**
- Marrow-reserve nadir: **0.773**
- Platelet-reserve nadir: **0.769**
- Immune-competence nadir: **0.761**
- Minimum computational stability: **0.708** at day **74.6**
- Time below the provisional `S = 0.72` boundary: **10.2 model-days**
- Model-generated 84-day cumulative infection-risk output: **0.342 (34.2%)**

At day 84, the model-generated latent-state probabilities were:

- Controlled: **0.345**
- Strained: **0.502**
- Overloaded: **0.154**

These values are mathematical outputs under the declared model assumptions. They are **not observed patient outcomes or validated clinical probabilities**.

---

## Sensitivity and Uncertainty

Global sensitivity was evaluated using **600 Latin-hypercube samples**.

Leading PRCCs for the model-generated 84-day infection-risk output included:

- Frailty: **0.963**
- Baseline hazard (`h0`): **0.963**
- Microbial input (`lambda_B`): **0.877**
- Neutrophil-related hazard coefficient (`beta_N`): **0.726**
- Relative dose intensity: **0.691**
- Microbial clearance (`clear_B`): **-0.763**

A separate **500-run stochastic perturbation analysis** produced a median model-generated risk of **34.9%**, with a 2.5th–97.5th percentile range of **27.3%–42.2%**.

This range represents simulation uncertainty under the declared perturbation assumptions and is **not a clinical confidence or prediction interval**.

---

## Structural Schedule Sensitivity

Model behaviour was materially dependent on the assumed timing of the normalized treatment-stress forcing.

| Structural forcing schedule | N nadir | Minimum S | Model-generated risk |
|---|---:|---:|---:|
| Reference q14 | 0.724 | 0.708 | 0.342 |
| Weekly after step-up | 0.555 | 0.552 | 0.461 |
| q21 after step-up | 0.790 | 0.765 | 0.306 |

These are **structural sensitivity experiments only**. They do not establish that one clinical dosing regimen is safer or more hazardous than another.

---

## Repository Structure

```text
bispecific-dynamic-risk-surface-model/
├── code/
│   ├── dynamic_risk_surface_model.py
│   ├── run_full_pipeline.py
│   └── make_figures.py
│
├── data/
│   ├── Published_Clinical_Anchors.csv
│   └── Model_Parameters.csv
│
├── figures/
│   └── manuscript figures
│
├── manuscript/
│   └── Dynamic_Risk_Surface_Full_Manuscript.pdf
│
├── result/
│   ├── Reference_Trajectory_Revised.csv
│   ├── Dynamic_Risk_Surface_Revised.csv
│   ├── LHS_600_Results_Revised.csv
│   ├── Global_Sensitivity_PRCC_Revised.csv
│   ├── Local_Elasticity_Revised.csv
│   ├── Scenario_Results_Revised.csv
│   ├── Stochastic_Perturbation_500_Revised.csv
│   ├── Synthetic_Model_Comparison_AIC_BIC_Revised.csv
│   ├── Numerical_Robustness_Revised.csv
│   ├── Exposure_Schedule_Sensitivity.csv
│   └── Analysis_Summary_Revised.json
│
├── .gitignore
├── CITATION.cff
├── LICENSE
└── README.md
