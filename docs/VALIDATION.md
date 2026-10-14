# Validation scope

This is an independently runnable **synthetic demonstration** of the declared
models. All packaged observation data and generator parameters are invented.
The release supports CPython 3.12; the local clean-install interpreter was
3.12.14 on macOS arm64. GitHub Actions is configured for Linux CPython 3.12;
execution of that hosted workflow is not claimed before upload.

The release checks deterministic data generation, the removable equal-rate
singularity in effect-site convolution, independently integrated effect-site
tail AUC, oral amount/AUC conservation, the turnover convolution identity,
and analytical quadratic values. The pipeline checks all fit/solver statuses,
matched amounts, exposure constraints, nested candidate values, independent
objectives, and injury-constrained dose targets and dose lower bounds.
Continuation comparisons independently enforce amount, floor and exposure
constraints and verify objectives against augmented ODE integration.

The analytical reference gives:

| Quantity | Value |
| --- | ---: |
| Resource cap B | 12.613903616328283 |
| Plateau C | 3 |
| Cessation time | 3.2046345387760944 |
| Flexible objective | 28.826261662949886 |
| Constant objective | 25.460925472460776 |
| Analytical stopping increment | 0.13947807099533094 |
| Delayed injury peak | 0.668740304976422 |
| Safe constant-family ceiling | 11.27093028990327 |

Numerical optimization may differ slightly across BLAS implementations and
platforms. Compare scientific outputs at rtol=1e-9, atol=1e-10 using the bundled
comparison script. PNG byte hashes can vary with fonts; they are not numerical
reproduction criteria. The release manifest supplies hashes for its artifacts.

## Limits

No original clinical aggregate, participant denominator or observed empirical
figure is reproduced by the synthetic cohort. Its results must not be substituted
for real-data results. In particular, this release does not run a complete PK/PD
uncertainty campaign, pooled absorption estimation, leave-one-profile-out
sensitivities, the 1,000-replicate residual/parameter bootstrap, rescue biomarkers
or Firth mucositis regression. The MTX example is a descriptive simulated course
series. These scopes are separate from full reproduction of this package's
complete synthetic workflow.

The quadratic frontier is an 80-bin-plus-atom numerical problem; it provides
candidate doses and supporting tangent-LP bounds, not unrestricted continuous
control optima. The cap margin near target 13 is numerical. Drug exposure
constraints are reference guardrails, not injury constraints or clinical safety
certificates. No counterfactual clinical regimen validation is performed.

The final manuscript comparison is recorded in [SCIENTIFIC_ALIGNMENT.md](SCIENTIFIC_ALIGNMENT.md),
with source hashes and analytical numerical checks in `SCIENTIFIC_ALIGNMENT.json`.
The clean-environment verification was repeated on 14 October 2026.
