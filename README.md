# LOGOS Theory B — Mathematical Optimization in Clinical Pharmacy

**From Mathematical Theory to Model-Based Regimen Design**

**Author: Vangelis D. Karalis.**

**Release v0.5.1 contains synthetic data only.** The research description below
introduces the manuscript's scope. This public package implements the core
model comparisons using newly simulated observations; it does not reproduce
all empirical analyses, tables or figures. See the [scientific alignment
review](docs/SCIENTIFIC_ALIGNMENT.md) for the comparison with the manuscript.

**LOGOS (from Ancient Greek λόγος: reason, rational account; Longitudinal Optimization of Global Outcomes and Strategies)** is a mathematical and computational framework for investigating whole-regimen pharmacokinetic/pharmacodynamic (PK/PD) treatment strategies.

This repository presents the applications of LOGOS Theory B to finite-course regimen design, integrating mathematical optimization, dynamical PK/PD modeling, individual pharmacological observations, and numerical verification.

The research evaluates the incremental value of additional treatment-design freedoms through nested mathematical strategy classes under matched drug amounts and longitudinal exposure constraints.

**Core research areas**

- **Constrained mathematical optimization:** Quantification of objective increments across nested admissible regimen classes, supported by analytical and numerical bounds.
- **Individualized PK/PD modeling:** Model fitting and simulation using remifentanil and warfarin pharmacological observations.
- **Temporal dose allocation:** Optimization of infusion profiles, administration timing, and dose redistribution over finite treatment horizons.
- **Conditional treatment cessation:** Isolation of stopping-specific value through mathematically matched continuation-constrained comparisons.
- **Longitudinal system dynamics:** Evaluation of residual drug exposure, delayed pharmacodynamic responses, and post-administration effects.
- **Computational verification:** Sensitivity analyses, bootstrap uncertainty estimation, independent numerical integration, and reproducibility checks.

The accompanying Python implementations provide a reproducible computational workflow for model fitting, constrained optimization, counterfactual regimen simulations, and mathematical verification.

**Scientific objective:** To determine how additional temporal treatment-design freedoms change attainable PK/PD model objectives under explicitly matched resource, exposure, and feasibility constraints.

The findings establish conditional mathematical and fitted-model optimization opportunities, rather than clinically validated improvements in therapeutic efficacy or safety.

*Research software for mathematical and computational pharmacology; not a validated clinical dosing system.*

## Scope of the synthetic release

In this release, pharmacological observations and fitted profiles are synthetic.
The computational-verification research area above describes the broader study;
the complete sensitivity campaigns and warfarin residual/parameter bootstrap
from the manuscript are not bundled. The included MTX bootstrap summarizes
synthetic available pairs. [Validation scope](docs/VALIDATION.md) lists exclusions.

Version **0.5.1**. Reproducible mathematical examples of temporal redistribution,
policy-conditional stopping, dose saving and target feasibility.

All bundled observations and profile parameters are **synthetic**, generated
from invented parameter ranges with seed `20261012`. No patient records or
parameters fitted to patient records are distributed. The synthetic cohort
contains four profiles per PK/PD application and twelve course-ordered MTX
biomarker profiles. These examples are model computations, not evidence of
clinical efficacy or safety and not dosing recommendations.

## Install and reproduce

Supported interpreter: **CPython 3.12 only** (`>=3.12,<3.13`). Python 3.12.14 is
the release verification interpreter. Other versions are not declared supported.

```bash
python3.12 -m venv .venv
# macOS / Linux:
source .venv/bin/activate
# Windows PowerShell instead: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.lock.txt
python -m pip install --no-build-isolation --no-deps .
python -m unittest discover -s tests -v
logos-b --output outputs/run1
logos-b --output outputs/run2
python scripts/compare_runs.py outputs/run1 outputs/run2
```

Each output directory must be empty. The complete workflow generates new data,
fits PK then PD, optimizes nested classes, evaluates continuation policies,
computes the theoretical dose frontier, generates tables/figures and runs
independent numerical checks. It does not consume frozen output tables.
`--profiles N` changes synthetic PK/PD cohort size; `--seed N` changes invented
profiles and observation noise. Release expectations apply to the default seed
and cohort size. Parameters and all optimized schedules are exported.

## Outputs

- `data/`: regenerated synthetic profiles, observations and MTX course records.
- `fitted_parameters.json`, `tables/fit_diagnostics.csv`: sequential fits;
  oral absorption and lag are known simulation constants.
- `tables/nested_classes.csv`: within-profile successive increments and
  supporting-plane lower bounds, in minutes or hours as explicitly labeled.
- `tables/continuation.csv`, `controls.json`: matched final-quarter policies
  and full schedules, including candidate gaps and independent ODE errors.
- `canonical.json`: analytical loading, maintenance and coast benchmark.
- `B_dose_frontier*`: fixed-target dose comparisons under an injury cap,
  including target residuals, dose lower bounds and complete-tail checks.
- `mtx_summary.json`: available-pair descriptive synthetic biomarker summary.
- `figures/`: four PNG figures of the synthetic and theoretical examples.
- `verification.json`, `run_manifest.json`: numerical checks, versions, seed
  and output hashes. Failure produces a nonzero exit status.

A complete fresh-run example is bundled in `examples/reference_run/`.
The local release checks are recorded in `docs/REPRODUCIBILITY.json`.
Verify the extracted release with `python scripts/verify_manifest.py`.

See [model definitions](docs/MODELS.md), [data dictionary](data/README.md),
[validation scope](docs/VALIDATION.md) and [release instructions](docs/RELEASE.md).

## License and repository topics

The MIT license applies to all code, documentation and generated synthetic data
in this release. Retain the copyright and permission notice when reusing them.

Repository topics (also provided as `github-topics.json`):
`logos`, `pharmacokinetics`, `pharmacodynamics`, `pkpd`, `dose-optimization`,
`dosing-regimens`, `mathematical-modeling`, `computational-pharmacology`,
`model-informed-precision-dosing`, `treatment-optimization`, `drug-safety`.
