# Scientific alignment with the manuscript

Final review: 14 October 2026. Release: v0.5.0.
Primary manuscript: *Applications of LOGOS Theory B to Finite-Course Regimen
Design* (`4_Theory_B_Apps.docx`). Its supplied supplement and `repository.zip`
were used as supporting sources. Exact source hashes and benchmark comparisons
are recorded in `SCIENTIFIC_ALIGNMENT.json`; source documents and real records
are excluded from this public package.

## Conclusion

The release implements the core mathematical definitions and a reproducible
synthetic workflow. The theoretical canonical benchmark agrees with the
source archive within absolute tolerance 1e-9. The public PK/PD observations,
parameters, cohort counts and biomarker values are invented. Consequently,
this review confirms alignment of the declared equations and comparison
contracts, **not full empirical reproduction of the manuscript**. Complete
synthetic workflow reproduction and manuscript reproduction are distinct.

## Equations and constraints

| Manuscript definition | Release implementation | Review |
| --- | --- | --- |
| Main Eq. (1), ΔJ = J larger − J smaller | `pipeline.py`, `tables/nested_classes.csv` | Differences calculated within each profile before aggregation |
| Main Eq. (2), positive biexponential PK convolution | `models.py`, `RemModel.step/matrix` | Equivalent two-component linear ODE/convolution |
| Main Eq. (3), effect-site delay and Hill EEG | `models.py`, `fitting.py`, `synthetic.py` | Same equations; parameters fitted to synthetic observations |
| Main Eq. (4), integral of Ce/(K+Ce) through W | `RemModel.objective` | Same exponent-one proxy and minutes; fitted Hill exponent is not substituted |
| Effect-site tail AUC from W to infinity | `RemModel.tail` | Analytical sum c_i(W)/k_i + Ce(W)/ke0 |
| Main Eq. (5), oral absorption/lag and PCA turnover | `WarModel.cp/response/objective` | Same model, mg/L and hours; objective horizon 120 h |
| Main Eq. (6), concave supporting-plane upper bound | `models.py`, `continuation.py` | Tangent maximization over sampled outer relaxation; numerical assumptions retained |
| Supplement S1, quadratic response and measure input | `benchmark.py` | Same canonical parameters and analytic reference |
| Supplement S2, turnover convolution identity | `WarModel.objective` | Same weighted coverage identity; also handles nonunit R0 |
| Fixed amount and nested families | `models.py`, `pipeline.py` | Infusion caps 1/1.25/1.5/2, 20 bins; oral 24/12/6 h grids through 72 h |
| Final-quarter continuation floors | `continuation.py` | Same ε=0.05/0.25, five infusion bins and four of thirteen oral times; other permissions held common |
| Injury and fixed-target frontier | `benchmark.py` | Same Z′=C−0.2Z, cap 2.5, initial atom, 80 bins, rate cap 20 and targets 10–13 |

Exposure constraints are reference guardrails, not calibrated clinical injury
constraints. Supporting gaps describe finite numerical families, not rigorous
interval certificates or unrestricted continuous-control optima. A negative
increment lower bound remains unresolved.

## Parameters and denominators

| Quantity | Manuscript | Public release |
| --- | --- | --- |
| Remifentanil / warfarin profiles | 61 / 32 fitted source profiles | 4 / 4 independently generated profiles by default |
| MTX longitudinal / rescue cohorts | 43 children / 67 separate records | 12 simulated course series; no rescue cohort |
| Warfarin pooled ka / lag | 2.203766 h⁻¹ / 0.876524 h | Known generator constants 1.2 h⁻¹ / 0.4 h; no pooled estimation |
| Infusion starting states | Zero in nominal campaign | Zero in default release |
| PCA starting fraction | Nominal 1; separate observed-baseline sensitivity | Generated baseline 1; no missing baseline simulation |
| PK/PD parameters and noise | Estimated from source observations | Invented ranges and seeded simulated noise, followed by sequential fits |
| Seed | Separate original campaign seeds | 20261012 for generator; independent seed offsets for synthetic summaries |

Generator ranges are explicit in `synthetic.py`: infusion duration 8–16 min,
reference rate 80–140 µg/min, k1 0.025–0.06 min⁻¹, k2 0.4–0.9 min⁻¹,
w1 0.001–0.002, w2 0.006–0.01, ke0 0.15–0.4 min⁻¹, K 5–12 ng/mL,
Hill exponent 1.2–2.2, E0=22 and response fraction 0.75. Oral dose is 40–70 mg,
k 0.012–0.025 h⁻¹, scale 0.09–0.16 L⁻¹, kout 0.015–0.035 h⁻¹ and
IC50 2–5 mg/L. These ranges are not inferred population estimates.
Units of w_i belong to the stated input-to-concentration convolution.

## All manuscript tables and figures

| Item | Final alignment and scope |
| --- | --- |
| Main Table 1, evidence/claim map | Mathematical, synthetic and clinical evidence distinctions documented; no empirical claims reproduced |
| Main Table 2, sources/denominators | Replaced by synthetic data dictionary; source counts are not the release cohort |
| Main Table 3, horizons | W, oral administration through 72 h, outcome through 120 h and post-input exposure checks align; original observation support is not reproduced |
| Main Table 4, total gains | Synthetic gains exported in `nested_classes.csv`; original empirical values/intervals not reproduced |
| Main Table 5, successive increments | Same paired estimand, new synthetic values; original distributions and summary table not reproduced |
| Main Table 6, robustness/post-input contribution | Full original sensitivity, pooled estimation, leave-one-profile-out and 1,000-replicate residual/parameter bootstrap campaigns not included |
| Main Table 7, continuation floors | Same policy and lower-bound definitions, new synthetic increments |
| Main Figure 1, paired-increment boxplots | Release uses two synthetic total-relative-gain line figures; different quantity, cohorts and presentation |
| Supplement Table S1, rescue biomarkers | Separate rescue/mucositis analysis not included |
| Supplement Table S2, stopping increments | Same policies; synthetic `continuation.csv` replaces empirical results |
| Supplement Table S3, dose frontier | Same theoretical problem; fresh numerical candidates and constant comparator |
| Supplement Table S4, MTX completeness | Synthetic missingness and available pairs exported; empirical counts not reproduced |
| Supplement Table S5, residuals and dose bounds | Fresh target/cap residuals, lower/upper dose bounds, ODE errors and controls exported |
| Supplement Table S6, grid/tolerance/margin sensitivity | Complete 160-bin and larger-margin campaign not included |
| Supplement Table S7, 240/336 h observation support | Not included; original empirical support cannot be transferred to invented observations |
| Supplement Figure S1, canonical and delayed injury | Same mathematical curves; release uses different colors and extends injury display to time 12 rather than 8; manuscript's cessation marker and early-window shading are absent |
| Supplement Figure S2, MTX trajectories and paired changes | Release shows synthetic course trajectories with missing-value gaps; original median curve and paired-change panel not reproduced |

All three embedded manuscript figures and all four fresh release figures were
visually inspected. Release plots are readable and explicitly label synthetic
or theoretical quantities. Visual inspection does not establish equality of
empirical distributions.

## Theoretical numerical comparison

The source canonical benchmark and the fresh public implementation agree on
B=12.613903616328283, plateau 3, cessation 3.2046345387760944,
J=28.826261662949886, Jref=25.460925472460776, and the analytical stopping
increment 0.13947807099533094. The delayed unit-bolus injury peak is
0.668740304976422 at time 2.0117973905426254. The safe constant-family
ceiling is 11.27093028990327.

Fresh fixed-target candidate doses round to the manuscript's six-decimal
values: 2.928587, 3.266454, 3.614249 and 3.973756 for targets 10, 11, 12
and 13. Targets 10 and 11 permit dose-saving comparison against a feasible
constant family; 12 and 13 demonstrate target restoration. Floating-point
candidate bounds may differ slightly from archived witnesses. The release
exports its own gaps and does not transfer the manuscript's reported maximum
gap as an unconditional guarantee for new runs or platforms.

## Reproducibility and remaining limits

See `REPRODUCIBILITY.json` for the fresh isolated CPython 3.12 verification,
two complete runs, independent checks and comparisons to bundled reference
results. This covers the entire **distributed synthetic workflow**.
Full empirical manuscript reproduction, the omitted campaigns above and
clinical validation remain outside this release. Any future manuscript or
source revision requires a new alignment review.
