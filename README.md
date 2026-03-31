# EFFS

Code and data for the paper:

> An Efficient Evolving Functional Fuzzy System for Online Function-on-Function Regression with Usage-Based Rule Pruning

## Repository contents

- `run_stock_paper.R`: reproduces the stock experiment
- `run_agriculture_paper.R`: reproduces the agricultural experiment
- `run_all_paper.R`: runs both paper experiments
- `plot_stock_rule_curves.R`, `plot_agriculture_rule_curves.R`: generate the rule-curve figures
- `validate_paper_results.R`: compares generated outputs with the bundled reference files
- `main_STR_FFS_more_pruning.R`: EEFFS core implementation
- `real_data/`: bundled datasets
- `reference_outputs/`: paper-matching CSVs and figures
- `FLM_interaction-master/`, `fAPLS-master/`: helper code for baseline methods

## Quick start

Install the required R packages, then run:

```r
Rscript run_stock_paper.R
Rscript run_agriculture_paper.R
Rscript plot_stock_rule_curves.R
Rscript plot_agriculture_rule_curves.R
Rscript validate_paper_results.R
```

Generated files are written to:

- `outputs/stock/`
- `outputs/agriculture/`

## Notes

- The paper runners include the baselines used in the final paper: `EFFNS`, `PLS-int`, `PLS-main`, `fAPLS`, and `SLASSO`.
- In the result CSVs, `EEFFS` and `EFFNS` use `train.real.mse.c` and `test.real.mse.c`, while the batch baselines use `train.real.mse` and `test.real.mse`.
