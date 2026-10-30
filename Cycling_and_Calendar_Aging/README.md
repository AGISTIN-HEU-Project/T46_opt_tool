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

**Caveat:** orange's real data barely ages over its measured range (~4.4% loss
over 1150 cycles), so its own points don't tightly constrain the fit's curve
shape — that's why its predicted EOL drifts furthest from the datasheet's own
trend. A regularized fix (pulling orange's shape toward blue/green's own
curvature) exists in the code, commented out, ready to re-enable if needed.

See `Datasheet_Validation_Hypertuned_Memo.pdf` for the full write-up.
