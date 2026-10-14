# Mathematical definitions

Every comparison holds model parameters, initial state, objective and exposure
guardrails fixed within each synthetic profile. The two drug examples have
separate units and are not pooled.

## Infusion example

For i=1,2, `c_i' = -k_i c_i + w_i u`, `Cp = c_1+c_2`,
`Ce' = ke0 (Cp-Ce)`, with zero initial states. Time is minutes; input is µg/min;
concentrations are ng/mL. `w_i` is the conversion weight in this convolution.
The synthetic EEG observation function is
`E0 [1 - f Ce^h/(K^h+Ce^h)]`.

The optimization proxy uses exponent one:
`J = integral_0^W Ce/(K+Ce) dt`, in minutes. EEG fitting estimates h, but h
is not substituted into this concave proxy. Twenty equal input bins have
nonnegative rate factors x. `sum(x)=20` preserves exact amount `rate*W`.
Nested classes use rate caps 1, 1.25, 1.5 and 2 times the reference rate.
Cap 1 with that equality is exactly the constant reference.

Reference full-path Cp/Ce peaks and the effect-site residual exposure integral
`integral_W^infinity Ce dt` are common guardrails.
The latter equals `sum_i c_i(W)/k_i + Ce(W)/ke0` and is computed analytically.
Peak searches include the post-input tail. Numerical constraint exchange adds
violating extrema. This is a floating-point full-path audit, not a rigorous
interval-arithmetic certificate or calibrated injury model.

## Oral turnover example

For dose fraction x_j at t_j, plasma contribution after lag is
`D x_j s ka [exp(-k z)-exp(-ka z)]/(ka-k)`, `z=t-t_j-lag>=0`,
and zero before lag. Concentration is mg/L and time is hours.
`R' = kout [K/(K+Cp)-R]`, `R(0)=1` for generated observations.
`J = integral_0^120 (1-R) dt`, in hours.

Known generator constants ka=1.2 h^-1 and lag=0.4 h are fixed in fitting.
Individual elimination/scale and turnover/K are refitted from noisy generated
observations. A single dose at zero is nested in nonnegative dose fractions on
0–72 h grids of 24, 12 and 6 h. `sum(x)=1` matches exact invented amount D.
All comparisons retain the single-dose reference's full-path Cp peak cap.

The convolution identity used for numerical optimization is
`J = integral_0^120 Cp/(K+Cp) [1-exp(-kout(120-t))] dt`
plus `(1-R0) [1-exp(-120 kout)]/kout` for nonunit initial R.
Independent augmented turnover ODE integration verifies the objective.

## Stopping attribution and numerical gaps

At epsilon=0.05 and 0.25, the smaller infusion family requires factors >=epsilon
in five final-quarter bins (reserved amount epsilon/4).
The smaller oral family requires fractions >=epsilon/13 at 54,60,66,72 h
(reserved amount 4 epsilon/13). The larger family removes only these floors.
These are hypothetical within-drug obligations, not cross-drug equivalents.
An early last input alone is not a stopping-value estimate.

For a concave objective, a supporting-plane maximization over sampled linear
constraints supplies a numerical upper gap g. The sampled polytope is an outer
relaxation of the continuous constraint set. The successive lower increment is
`J_large-J_small-g_small`; continuation bounds additionally subtract independent
ODE discrepancies from both candidates. A negative lower bound is unresolved,
not rounded into a positive claim. Quadrature and floating-point assumptions
remain; there is no unrestricted-control optimality claim.

## Quadratic theoretical benchmark

`dC=-C dt+dmu`, `g(C)=4C-C^2/2`, outcome horizon T=4, zero initial state,
nonnegative input measures and resource cap B=12.613903616328. An initial
loading atom is allowed. Shadow price q=1 gives plateau c=3, coast length
`acosh(4/3)` and cessation `4-acosh(4/3)`. The no-loading constant reference
matches the realized amount. A separate mandatory minimum rate epsilon=0.5
isolates its analytical stopping increment.

A separate injury example uses `Z'=C-0.2Z`, `Z(0)=0`. A unit bolus gives
`Z=(exp(-0.2t)-exp(-t))/0.8`, demonstrating a delayed peak.
The minimum-dose frontier targets J=10,11,12,13 with injury cap 2.5, the same
available cap B, an initial atom and 80 equal infusion bins with upper bounds
20. The code minimizes atom plus integrated input, checks continuous injury
extrema through time 20 and the decreasing two-exponential tail thereafter,
and computes a tangent-LP dose lower bound. Independent ODE integration checks
benefit and injury. An infeasible constant-family target has undefined dose
saving, represented as null, and is labeled target restoration.

## Biomarker context

Synthetic MTX course values describe persistence, available measurements and
paired course-four minus course-one changes, with 3,000 seeded pair bootstrap
resamples. Missing observations remain unavailable. Course is categorical;
connecting lines break at missing values. This component neither estimates a
causal treatment effect nor identifies an optimal cessation time.
