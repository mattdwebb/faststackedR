# Earlier development notes

These are the previous landing-page notes, preserved on 19 November 2026 for
provenance. They may describe superseded restrictions. Use the [current README](../README.md),
[inference guide](../INFERENCE.md), and [change history](../CHANGELOG.md) for current behavior.

# faststackedR

Native base-R corrective-weighted Wing/Freedman/Hollingsworth stacked DiD.
MIT, copyright 2026 Matt Webb. Development reference; see VALIDATION.md.
Install with `R CMD INSTALL .`. No non-base runtime dependency.

```r
library(faststackedR)
r <- fit(data, "y", "time", unit="id", gvar="g", cluster="id",
         event_window=c(-2,1), inference="native")
print(c(estimate=r$estimate, native_se=r$native_se))
```

Run `Rscript examples/example.R` after installation, or `Rscript tests/validate.R`
from the repository root. Optional official comparisons require STACKED_ORACLE
pointing at the frozen author checkout, data.table, and fixest.

`fit(data,y,time,unit=NULL,gvar=NULL,treatment=NULL,cluster=NULL,covariates=NULL,
event_window=NULL,data_mode='auto',...)` uses column names. Missing/+Inf cohort
means never; zero is a real date. Panel needs unit; RCS needs gvar and cluster.
Auto chooses panel if unit is supplied, otherwise RCS. Numeric covariates enter
at their contemporaneous values, labelled as an adjustment extension.

Returns list with estimate, target_label, native_se, native_p, native_df (df
alias), events/event_study (data.frame), event_covariance (matrix), stack,
cohort weights, support and covariate metadata. Point-only native fields are
NA. Jackknife adds jk_se, jk_p and deletion details. RI adds method_p and detail
lists. No source data mutation or calls to Python/Stata.

See METHOD_SPECIFICATION.md, INFERENCE.md and THIRD_PARTY_NOTICES.md.
