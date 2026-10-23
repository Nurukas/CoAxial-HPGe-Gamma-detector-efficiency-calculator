# Example output

Produced by `hpgeff check / run / transfer examples/demo_detector2.yaml` (10⁶ histories per
energy). The transferred efficiencies use `reference_curve_SYNTHETIC.csv`, which is made-up
test data — they show the format only, **not** your sample's efficiency. Re-run `transfer`
with the measured reference curve to get real numbers.

| file | what it is |
|---|---|
| `check_report.txt` | pre-flight report (volumes, overlaps, nuclear data, short debug transport) |
| `glass_in_pmma_cup_xz.png` + legend | OpenMC slice through the axis |
| `*_mc_efficiency.csv`, `mc_efficiency.png` | Monte Carlo FEP efficiencies |
| `glass_in_pmma_cup_efficiency_transferred.csv`, `efficiency_transferred.png` | transfer to the (synthetic) reference |
