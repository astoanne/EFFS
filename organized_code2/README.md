# organized_code2

This folder is a paper-faithful rebuild of the reproducible code path for the WCCI paper results.

What is included:

- `run_stock_paper.R`: stock experiment for `AAPL -> NVDA`
- `plot_stock_rule_curves.R`: draws paper-style antecedent/consequent rule curves for the final stock EEFFS rule base
- `run_agriculture_paper.R`: agricultural experiment for `solar_radiation -> bare_soil_temp`
- `plot_agriculture_rule_curves.R`: draws paper-style antecedent/consequent rule curves for the agricultural EEFFS run
- `main_STR_FFS_more_pruning.R`: EEFFS core used by the paper-style runs
- `real_data/`: bundled datasets needed by the paper runners
- `FLM_interaction-master/`: local helper code for the PLS baseline
- `fAPLS-master/`: local helper code for the fAPLS baseline
- `reference_outputs/`: exact historical CSVs and PNGs that match the paper tables/figures
- `validate_paper_results.R`: compares new outputs against the exact paper reference files
- `run_all_paper.R`: convenience runner for both experiments

This folder is intended to be self-contained as a GitHub folder.
You can move or publish `organized_code2` by itself without keeping the larger `EFFS` parent folder.

The paper runners include only the baselines used in the final paper:

- `EFFNS`
- `PLS-int`
- `PLS-main`
- `fAPLS`
- `SLASSO`

`OVGLASSO` is not executed by the `organized_code2` runners.

Important note:

- In the paper tables, the `EEFFS` and `EFFNS` rows come from `train.real.mse.c` and `test.real.mse.c`.
- The batch baselines come from `train.real.mse` and `test.real.mse`.

How to run:

```r
Rscript run_stock_paper.R
Rscript plot_stock_rule_curves.R
Rscript run_agriculture_paper.R
Rscript plot_agriculture_rule_curves.R
Rscript validate_paper_results.R
```

On Windows, if `Rscript` is not on your `PATH`, you can run:

```powershell
& "C:/Program Files/R/R-4.4.1/bin/Rscript.exe" "C:/path/to/organized_code2/run_stock_paper.R"
& "C:/Program Files/R/R-4.4.1/bin/Rscript.exe" "C:/path/to/organized_code2/plot_stock_rule_curves.R"
& "C:/Program Files/R/R-4.4.1/bin/Rscript.exe" "C:/path/to/organized_code2/run_agriculture_paper.R"
& "C:/Program Files/R/R-4.4.1/bin/Rscript.exe" "C:/path/to/organized_code2/plot_agriculture_rule_curves.R"
& "C:/Program Files/R/R-4.4.1/bin/Rscript.exe" "C:/path/to/organized_code2/validate_paper_results.R"
```

Do not run the older `organized_code/run-stock.R` or `organized_code/run-agriculture.R` if you want the paper-faithful results. Those belong to the newer reorganized workflow and can differ from the paper.

Generated files are written to:

- `outputs/stock/`
- `outputs/agriculture/`

The stock plotting script writes:

- `outputs/stock/stock_rule_curves_final.png`
- `outputs/stock/stock_rule_curves_max_rules.png`
- `outputs/stock/stock_rule_curves_max_rules_compact.png`
- `outputs/stock/stock_rule_curve_plots.csv`

The agricultural plotting script writes:

- `outputs/agriculture/agriculture_rule_curves_final.png`
- `outputs/agriculture/agriculture_rule_curves_max_rules.png`
- `outputs/agriculture/agriculture_rule_curves_max_rules_compact.png`
- `outputs/agriculture/agriculture_rule_curve_plots.csv`

Reference files copied from the exact paper-matching runs:

- `reference_outputs/stock_paper_reference.csv`
- `reference_outputs/agriculture_paper_reference.csv`
- corresponding rule-trajectory PNGs

If a fresh rerun differs from the reference CSVs, `validate_paper_results.R` will show the mismatch. This can happen when package behavior changes across environments, even when the script settings are the same.

Only the local file dependencies are bundled here. You still need the required R packages installed in your environment.
