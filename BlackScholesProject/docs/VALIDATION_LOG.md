# Validation log

## 2026-11-24 — planted reaction-sign fault

The first frozen test required the wrong-sign reaction run at `Nt=640` to
have more than three times the correct solver's target error. It failed:

| Run | Target price | Absolute error |
|---|---:|---:|
| Correct reaction sign | 3.88713341646 | 0.04214262368 |
| Planted wrong sign | 3.91117272103 | 0.01810331911 |

This does not validate the defect. The two errors cancel at this coarse
resolution: the correct discretization underprices the option, while the wrong
reaction sign raises the value toward the analytical target.

The protocol was therefore changed rather than the failure being hidden. The
falsification test now uses the already-declared finest paper resolution,
`Nt=1280`, and requires both:

- the faulty target error to exceed the correct target error;
- the faulty and correct target prices to differ by more than `2e-2`.

At this resolution the preliminary observed values were `0.01041522088` and
`0.01332899458`, with a price separation of `0.02374421546`. The validation
plan records `protocol_frozen_before_outcomes=false`; later outcomes cannot be
described as prospectively frozen.
