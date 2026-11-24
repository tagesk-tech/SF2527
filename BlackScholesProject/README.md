# SF2527 Black--Scholes project

Reproducible MATLAB implementation and audit of Wang et al. (2025),
"A second-order finite difference method for the Black--Scholes model
without far-field boundary conditions."

The project tests the paper's two-sweep Saul'yev/ADE method on a grid that
loses one rightmost point per time step. It keeps paper reproduction separate
from independent verification: the paper's printed exact European-call value
is inconsistent with the standard Black--Scholes formula, so all verification
uses the corrected analytical oracle.

## Run

From MATLAB, in the repository root:

```matlab
run('tests/run_tests.m')
run('experiments/run_all.m')
```

Generated CSV tables and figures are written to `results/` and intentionally
excluded from version control. The scripts use only base MATLAB functions.

## Repository map

- `src/ade_shrinking_grid.m`: the paper's opposite sweeps and averaging.
- `src/forward_euler_shrinking_grid.m`: independent explicit baseline.
- `src/black_scholes_exact.m`: analytical prices and call/put Greeks.
- `experiments/`: paper reproduction, temporal convergence, spatial-operator
  consistency, the deliberately retained coupled full-solver scan, put
  validation, and a stability stress test.
- `tests/run_tests.m`: exact-oracle, reproduction, parity, refinement, and
  planted-fault checks.
- `docs/VALIDATION_LOG.md`: the recorded falsification-test correction.
- `validation-plan.json`: method-level verification protocol and claim status.

Course notes, presentation scope, results interpretation, and the contribution
template live in the SF2527 section of the Obsidian vault. This repository owns
the runnable MATLAB implementation and its technical validation artifacts.

## Scientific status

`UNTESTED` means a script has not been run. `PRELIMINARY`, `CONDITIONAL`,
`VERIFIED`, and `REFUTED` are used only as defined in `validation-plan.json`.
Agreement with a plot or an old output is not treated as proof of correctness.

## Source boundary

The journal PDF is not redistributed in this repository. The DOI is
<https://doi.org/10.1016/j.jfs.2025.101477>. Equations (8)--(10) were
implemented from the paper; the appendix code was used only as an audit target.
