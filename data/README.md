# Synthetic data dictionary

All values are newly simulated; IDs begin with SYN_. License: MIT.
No clinical records, source identifiers, copied observation arrays or fitted
patient parameters are included. Invented ranges are computational examples;
they are not population parameter estimates.

`profiles.json`: generator seed, synthetic flag and two parameter dictionaries.
Remifentanil keys: ID, T_min, rate_ug_min, k1_min, k2_min, w1, w2, ke0_min,
K_ng_ml, Hill, E0 and response_fraction. Warfarin keys: ID, dose_mg, k_h,
scale_L_inv, ka_h, lag_h, kout_h and IC50_mg_L. See docs/MODELS.md for units.

`observations.csv`: ID, time, kind, value, unit. Remifentanil time is minutes;
warfarin time is hours. Cp has drug-specific units; EEG is a synthetic readout;
PCA is percent. Cp noise is multiplicative lognormal. EEG and PCA noise is
additive Gaussian. Baseline PCA is fixed at 100. Every profile is generated
independently of any real participant.

`mtx_courses.csv`: synthetic ID, course 1–4, sampling_day, MTX_PG2_5_nmol_L.
Blank values are intentionally missing. Generator uses invented baselines,
logarithmic accumulation and observation noise; it is not a calibrated MTX
PK model. Cohorts are not linked.

The bundled files are the default generator output. Every pipeline execution
regenerates them from the source code and seed.
