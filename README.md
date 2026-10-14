# Mathematical modelling of infectious diseases in R

**From pen-and-paper SEIR models to calibration, contact matrices, variants of concern and stochastic outbreaks: a reproducible R workflow.**

[![R tests](https://github.com/epic-demiologist/infectious-disease-modelling-r/actions/workflows/tests.yml/badge.svg)](https://github.com/epic-demiologist/infectious-disease-modelling-r/actions/workflows/tests.yml)

This repository works through the course *Mathematical modelling of infectious diseases* (Dr. Veronika Jäger, University of Münster; taught in the MSc Epidemiology programme at LMU Munich, 2026). In the course, models were built in Berkeley Madonna. Here, every lecture and practical is rebuilt from scratch in **R**, with:

- unit tests that check the code against the worked examples on the lecture slides,
- a separate **Python** implementation that cross-checks the R results,
- explicitly labelled assumptions in a single parameter file.

The running example is the one used in the course: *"We are back in March 2020. What should we expect from SARS-CoV-2, and which control measures will work?"*

> **Disclaimer.** This is a learning and portfolio project. Parameters flagged `ASSUMPTION` in [`R/parameters.R`](R/parameters.R) are illustrative. They are not estimates for any real population or policy decision.

---

## Contents

| # | Topic (course session) | Script | Methods |
|---|---|---|---|
| 01 | Difference vs. differential equations (L2, L3, P1) | [`01_difference_vs_differential.R`](analysis/01_difference_vs_differential.R) | Discrete-time SEIR; Euler vs. Runge-Kutta 4 vs. adaptive LSODA |
| 02 | Basic SEIR model for SARS-CoV-2 (P2) | [`02_seir_sarscov2.R`](analysis/02_seir_sarscov2.R) | ODE SEIR, final-size equation, herd immunity threshold, sensitivity analysis |
| 03 | ICU capacity (P3) | [`03_icu_capacity.R`](analysis/03_icu_capacity.R) | SEIR extended with a pre-ICU and ICU pathway |
| 04 | Interventions and vaccination (L4, P4) | [`04_interventions.R`](analysis/04_interventions.R) | Time-varying contact rates (NPIs), R(t); endemic SIR with births, deaths and vaccination; critical vaccination coverage |
| 05 | Non-random mixing (L5, P5) | [`05_contact_matrices.R`](analysis/05_contact_matrices.R) | Age-structured SEIR, WAIFW/contact matrix, next-generation matrix, school closure vs. shielding |
| 06 | Variants of concern (P6) | [`06_variants_of_concern.R`](analysis/06_variants_of_concern.R) | Two-strain SEIR; estimating growth advantage from simulated sequencing data with logistic regression |
| 07 | Calibration and validation (L7) | [`07_calibration_validation.R`](analysis/07_calibration_validation.R) | Least-squares (RSS) fitting, RSS profile, parametric bootstrap, out-of-sample validation |
| 08 | Stochastic vs. deterministic models (L8) | [`08_stochastic_models.R`](analysis/08_stochastic_models.R) | Reed-Frost (individual-based, chain-binomial, exact), Gillespie SIR, probability of a major outbreak |

The model functions are in [`R/models.R`](R/models.R). Epidemiological quantities (HIT, final size, next-generation matrix, RSS) are in [`R/utils.R`](R/utils.R). All scenario parameters are in [`R/parameters.R`](R/parameters.R).

---

## Quick start

```bash
git clone https://github.com/epic-demiologist/infectious-disease-modelling-r.git
cd infectious-disease-modelling-r

# R >= 4.1 with: deSolve, ggplot2, dplyr, tidyr, scales, ggrepel, testthat
Rscript -e 'install.packages(c("deSolve","ggplot2","dplyr","tidyr","scales","ggrepel","testthat"))'

Rscript run_all.R                                   # all results/ and figures/ (~2 min)
Rscript -e 'testthat::test_dir("tests/testthat")'   # 10 test blocks, 28 expectations
python python/seir_crosscheck.py                    # independent Python check (numpy, scipy)
```

Optional: if [`socialmixr`](https://cran.r-project.org/package=socialmixr) is installed, script 05 automatically uses the empirical **POLYMOD contact matrix for Germany** (Mossong et al. 2008) in place of the built-in illustrative matrix.

---

## Model structure

<p align="center"><img src="figures/model_diagram.svg" width="720" alt="SEIR model with ICU pathway"></p>

$$
\frac{dS}{dt} = -\beta\, c(t)\, I S,\qquad
\frac{dE}{dt} = \beta\, c(t)\, I S - fE,\qquad
\frac{dI}{dt} = fE - rI,\qquad
\frac{dR}{dt} = (1-p_{ICU})\, rI + U/d_{ICU}
$$

$\beta = R_0/(D N)$, $f = 1/\text{pre-infectious period}$, $r = 1/\text{infectious period}$, and $c(t)$ is the share of pre-pandemic contacts that remain under NPIs. The ICU pathway (dashed box in the diagram) uses $dP/dt = p_{ICU}\,rI - P/d_{pre}$ and $dU/dt = P/d_{pre} - U/d_{ICU}$.

---

## Results

### 01 · Step size matters

The pen-and-paper example from Practical 1 is reproduced exactly (S₁ = 9,900, E₁ = 100). With that example's parameters (implied R₀ = 1,000), the daily difference equations break down: S becomes negative on day 27. For the SARS-CoV-2 model, a 1-day difference equation overestimates peak prevalence by **4.1%** and delays the peak by 6 days. RK4 with the same 1-day step is within 0.001% of the adaptive solver.

<img src="figures/01_integrators.png" width="720">

### 02 · What to expect in March 2020

Assumptions: R₀ = 2.5, pre-infectious period 5.5 d, infectious period 3 d, N = 83.2 million. The unmitigated epidemic peaks on **day 98** with 8.1% of the population infectious at the same time. Final size is **89.3%**, which matches the analytical final-size equation *z = 1 − e^(−R₀z)* to one decimal place. This is **29 percentage points beyond** the herd immunity threshold of 60%, because the epidemic overshoots.

<img src="figures/02_seir_baseline.png" width="720">

The peak timing is very sensitive to the two parameters that were least certain early in the pandemic:

<img src="figures/02_sensitivity_peak.png" width="620">

### 03 · ICU capacity

With illustrative ICU assumptions (0.5% of infections need ICU care, 7-day delay, 10-day stay, 10,000 beds), demand peaks at **8.5× capacity** and stays above capacity for **75 days**.

<img src="figures/03_icu_unmitigated.png" width="720">

### 04 · Interventions

| Scenario | Peak ICU beds | Days over capacity | Attack rate |
|---|---:|---:|---:|
| No intervention | 84,731 | 75 | 89.3% |
| 40% contact reduction from day 30 on | 28,510 | 95 | 58.3% |
| 60% reduction day 30–150, then lifted | 84,447 | 75 | 89.2% |
| 60% day 30–150, then 40% sustained | 28,288 | 96 | 58.2% |

Strict NPIs hold R(t) at 1. Once they are lifted, the population is almost as susceptible as before, so the wave is only postponed. NPIs buy time (for example, for vaccines). They do not lower the final burden unless they are sustained.

<img src="figures/04_npi_icu.png" width="720">
<img src="figures/04_npi_rt.png" width="720">

**Vaccination in an endemic setting.** The model is a measles-like infection with births and deaths: R₀ = 15, life expectancy 70 years. The simulated equilibria match the analytical results from the lecture, *S\* = 1/R₀* and *I\* = μ(1 − 1/R₀ − p)/(r + μ)*, to six significant figures. Coverage below *p_crit = 1 − 1/R₀ = 93.3%* leads to a "honeymoon period" with damped epidemic cycles and then a lower endemic level. Coverage above p_crit leads to elimination.

<img src="figures/04_vaccination_endemic.png" width="720">
<img src="figures/04_critical_coverage.png" width="620">

### 05 · Who acquires infection from whom

The model has four age groups with a reciprocal contact matrix, and q is calibrated so that the dominant eigenvalue of the next-generation matrix equals 2.5.

| Scenario | R₀ | Attack 0–19 | 20–39 | 40–59 | 60+ |
|---|---:|---:|---:|---:|---:|
| Baseline | 2.50 | 93.8% | 89.0% | 85.7% | 61.8% |
| School closure (−70% contacts among 0–19) | 2.18 | 81.3% | 88.1% | 84.6% | 60.5% |
| Shielding 60+ (−60% of their contacts) | 2.42 | 92.5% | 86.7% | 82.2% | **25.7%** |

Shielding barely changes R₀, but it reduces infections in the oldest group by more than half. A single summary number such as R₀ hides who is protected.

<img src="figures/05_contact_matrix.png" width="480"> <img src="figures/05_age_structured_scenarios.png" width="480">

### 06 · Variants of concern

NPIs are held constant at a level that keeps the wild type at R = 0.95. A variant with 50% higher transmissibility (R = 1.43) grows underneath the declining epidemic. Total cases fall for three months and then turn upward on **day 97**, although nothing about the measures has changed. A logistic regression on simulated weekly sequencing data (400 cases/week) recovers the growth advantage: **0.052/day (95% CI 0.046–0.059)** against a model truth of 0.047/day.

<img src="figures/06_voc_incidence.png" width="720">
<img src="figures/06_voc_share_logistic.png" width="720">

### 07 · Calibration ≠ validation

The model is fitted by least squares to days 0–40 of synthetic reported-case data (true R₀ = 2.7). The estimate is R₀ = **2.69** (bootstrap 95% CI 2.61–2.78), and the initial number infectious is recovered as 51 (true 50). Forecasts for the held-out days 41–60 have a mean absolute percentage error of 1.5%. Only 80% of held-out points fall inside the 95% Poisson prediction interval, because that interval ignores parameter uncertainty. This shows why validation has to use data that were not used for fitting.

<img src="figures/07_calibration_validation.png" width="720">
<img src="figures/07_rss_profile.png" width="620">

### 08 · Stochastic vs. deterministic

**Reed-Frost (Lecture 8 example).** One measles case is introduced among 10 susceptibles, with p = 0.15. The individual-based and chain-binomial simulations (10,000 runs each) agree with the exact final-size distribution, which is computed by enumerating the Markov chain. The outcome is bimodal: in 19.7% of introductions nobody else is infected, while 9 of 10 is the most likely outbreak size.

<img src="figures/08_reed_frost.png" width="720">

**Gillespie SIR.** With R₀ > 1 a large epidemic is likely but not certain. The simulated probability of a major outbreak matches branching-process theory, 1 − 1/R₀: for example, 0.494 vs. 0.500 at R₀ = 2. A deterministic model cannot show early extinction.

<img src="figures/08_major_outbreak.png" width="620">
<img src="figures/08_stochastic_vs_deterministic.png" width="720">

---

## Verification

| Check | How |
|---|---|
| Lecture worked examples (pen-and-paper SEIR, HIT for influenza and varicella, mean age at infection, Reed-Frost step probabilities) | `tests/testthat/test-models.R` |
| Mass conservation; ODE final size = analytical final size | unit tests |
| RK4 converges to the adaptive solver; Euler error decreases with step size | unit tests |
| Endemic equilibrium = analytical S\*, I\* | unit tests + `results/04_endemic_equilibrium_check.csv` |
| Next-generation matrix reduces to q·c·D under proportionate mixing | unit tests |
| Exact Reed-Frost distribution sums to 1, P(no spread) = (1 − p)^S₀ | unit tests |
| R and Python give the same SEIR peak and final size | `python/seir_crosscheck.py` |

Every push runs the tests and the Python cross-check on GitHub Actions ([`.github/workflows/tests.yml`](.github/workflows/tests.yml)).

## Repository structure

```
├── R/
│   ├── models.R          # SEIR (difference, ODE, ICU), endemic SIR+vaccination,
│   │                     # age-structured SEIR, two-strain SEIR, Reed-Frost, Gillespie
│   ├── utils.R           # R0/HIT/final size, NGM, RSS, plotting theme
│   └── parameters.R      # all scenario parameters (assumptions flagged)
├── analysis/             # 01-08, one script per course session
├── tests/testthat/       # 10 test blocks, 28 expectations
├── python/               # independent cross-check of the R implementation
├── results/              # CSV outputs (regenerated by run_all.R)
├── figures/              # PNG figures (regenerated by run_all.R)
└── run_all.R
```

## References

- Vynnycky E, White RG. *An Introduction to Infectious Disease Modelling.* Oxford University Press, 2010.
- Anderson RM, May RM. *Infectious Diseases of Humans: Dynamics and Control.* Oxford University Press, 1991.
- Mossong J, et al. Social contacts and mixing patterns relevant to the spread of infectious diseases. *PLoS Med* 2008;5(3):e74.
- Ferguson NM, et al. Report 9: Impact of non-pharmaceutical interventions (NPIs) to reduce COVID-19 mortality and healthcare demand. Imperial College London, 16 March 2020.
- Halloran ME, Struchiner CJ, Longini IM. Study designs for evaluating different efficacy and effectiveness aspects of vaccines. *Am J Epidemiol* 1997;146(10):789–803.
- Habbema JDF, et al. Model building on the basis of empirical data. 1996 (model development steps as presented in the course).
- Soetaert K, Petzoldt T, Setzer RW. Solving differential equations in R: package deSolve. *J Stat Softw* 2010;33(9):1–25.

Course material: *Mathematical modelling of infectious diseases*, Dr. Veronika Jäger (University of Münster), Munich 2026. Lecture slides are not redistributed here.

## Author

**Bikram Kumar Singh**, MSc Epidemiology (LMU Munich) · [GitHub @epic-demiologist](https://github.com/epic-demiologist)

Released under the MIT License.
