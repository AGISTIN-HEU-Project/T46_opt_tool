# Battery Ageing Model Validation — A123 ANR26650M1

Validates a hypertuned, cross-validated per-curve cycling-ageing model against the
A123 ANR26650M1 (LFP) cell's own manufacturer datasheet.

Based on: J. Nájera, J.R. Arribas, R.M. de Castro, C.S. Núñez, "Semi-empirical
ageing model for LFP and NMC Li-ion battery chemistries," *Journal of Energy
Storage*, vol. 72, 108016, 2023.

## Files

| File | What it is |
|---|---|
| `datasheet_validation.py` | The script. Fits the model to the datasheet's three real curves and generates the chart + CSV below. |
| `datasheet_validation_hypertuned.png` | The validation chart (model fit vs. real datasheet points, per curve). |
| `datasheet_validation_hypertuned.csv` | The chart's underlying data (fitted SOH curves, point by point). |

## Regenerating

```bash
python3 datasheet_validation.py
```

Regenerates `datasheet_validation_hypertuned.png` and `.csv` from scratch. No
other inputs are needed — the three real datasheet curves are hardcoded in the
script (`PROTOCOLS` dict).

## The three real datasheet curves

| Curve | Charge / Discharge | Temperature | Real measured points |
|---|---|---|---|
| Orange | +2.5 A / −2.5 A | 23°C | 8 |
| Green | +3 A / −5 A | 45°C | 7 |
| Blue | +3 A / −20 A | 23°C | 7 |

## Key result

Each curve is fit independently (own k_cyc, β, n — no shared constants) and
checked against the datasheet's own trend at 80% SOH (end-of-life):

| Curve | Model EOL | Datasheet's own EOL | Gap at 80% SOH |
|---|---|---|---|
| Blue | 1645 | 1617 | 28 cycles (best) |
| Green | 2278 | 2182 | 96 cycles |
| Orange | 4754 | 3512 | 1242 cycles (worst) |

## Where this model works best

The model is most accurate for **blue-type** operation: slow-to-moderate
charging and fast discharging at room temperature (EOL gap 28 cycles, trend
RMSE 0.17%).

| Condition | Best-fit range (blue-type) | Measured datasheet point |
|---|---|---|
| Cell / chemistry | A123 ANR26650M1, LFP, 2.5 Ah, 3.3 V | A123 ANR26650M1 |
| Charge current | Below 5 A (below 2C) | +3 A (1.2C) |
| Discharge current | Above 5 A (above 2C), up to 20 A (8C) | −20 A (8C) |
| Effective C-rate | Around 2C (one full cycle in about 1 h) | 2.09C (0.96 h per cycle) |
| Temperature | Room temperature, ≈ 20–25°C | 23°C |
| Depth of discharge | Full cycles, 100% DOD | 100% DOD |
| SOH range | 100% down to 80% SOH (end of life) | 100% → 82.6% measured |

Typical match: a BESS charged gently (e.g. from PV) and discharged hard (peak
shaving, pump start-up, high-price market hours).

Between −5 A and −20 A discharge at 23°C the datasheet has no measured curve;
values in that range are the model's interpolation.

**Caveat:** orange's real data barely ages over its measured range (~4.4% loss
over 1150 cycles), so its own points don't tightly constrain the fit's curve
shape — that's why its predicted EOL drifts furthest from the datasheet's own
trend.

**Fix:** a regularized fit for orange is already in the code (commented out). It
pulls orange's bend n toward blue/green's average n (1.43) by minimizing
`RMSE + λ·(n − 1.43)²`, cutting orange's gap from 1242 to 54 cycles (λ = 3.0).
To use it, uncomment `orange_regularized_addendum()` at the bottom of the script.

See `Datasheet_Validation_Hypertuned_Memo.pdf` for the full write-up.
