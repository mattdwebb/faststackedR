library(faststackedR)
base <- if(dir.exists("tests")) "tests" else "."
d <- read.csv(file.path(base,"singleton_cohorts.csv"))
expected <- read.csv(file.path(base,"singleton_reference.csv"))
for(cv in 0:1) {
  args <- list(y="y",time="time",gvar="g",cluster="state",data_mode="rcs",event_window=c(-2,1),covariates=if(cv) "x" else NULL,weights="w")
  z <- do.call(fit,c(list(data=d,inference="jackknife"),args))
  ref <- subset(expected,adjusted==cv)
  physical <- vapply(sort(unique(d$state)),function(i) do.call(fit,c(list(data=d[d$state!=i,],inference="point"),args))$estimate,numeric(1))
  stopifnot(abs(z$estimate-ref$estimate[1])<1e-9,abs(z$jk_se-ref$jk_se[1])<1e-9,
            max(abs(z$jackknife$estimates-ref$estimate[-1]))<1e-9,
            max(abs(physical-z$jackknife$estimates))<1e-11,z$jackknife$df==7)
}
cat("SINGLETON_COHORT_R_PYTHON_PARITY_PASS\n")
