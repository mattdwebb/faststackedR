.check <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

.point <- function(d, lo, hi, mode, control, xs, native=TRUE) {
  events <- seq.int(lo, hi)
  cohorts <- sort(unique(d$g[is.finite(d$g)]))
  stacks <- list()
  dropped <- numeric()
  for (a in cohorts) {
    if (a + lo < min(d$t) || a + hi > max(d$t)) {
      dropped <- c(dropped, a)
      next
    }
    clean <- d$g > a + hi
    if (control == "never_treated") clean <- is.infinite(d$g)
    if (control == "not_yet_treated") clean <- clean & is.finite(d$g)
    s <- d[(d$g == a | clean) & d$t >= a + lo & d$t <= a + hi, , drop=FALSE]
    s$a <- rep(a, nrow(s)); s$e <- s$t - a; s$z <- as.numeric(s$g == a)
    if (length(unique(s$z)) < 2) {
      dropped <- c(dropped, a)
      next
    }
    if (!all(vapply(events, function(e) all(c(0, 1) %in% s$z[s$e == e]), logical(1)))) {
      dropped <- c(dropped, a)
      next
    }
    if (mode == "panel") {
      .check(all(table(s$u) == length(events)), "UNSUPPORTED: unbalanced panel within event window")
      .check(all(vapply(split(s$w, s$u), function(w) length(unique(w)) == 1, logical(1))),
             "UNSUPPORTED: panel weights must be time invariant within stack")
    }
    stacks[[length(stacks) + 1L]] <- s
  }
  .check(length(stacks) > 0, "IDENTIFICATION_FAILURE: no identified cohorts in window")
  s <- do.call(rbind, stacks)
  ta <- tapply(s$w * s$z, s$a, sum)
  shares <- ta / sum(ta)
  s$q <- s$w
  if (mode == "panel") {
    ca <- tapply(s$w * (1 - s$z), s$a, sum)
    ratio <- shares / (ca / sum(ca))
    s$q[s$z == 0] <- s$w[s$z == 0] * ratio[as.character(s$a[s$z == 0])]
  } else {
    cell <- interaction(s$a, s$e, s$z, drop=TRUE)
    total <- interaction(s$e, s$z, drop=TRUE)
    s$q <- s$w * shares[as.character(s$a)] * ave(s$w, total, FUN=sum) / ave(s$w, cell, FUN=sum)
  }
  s$q <- as.numeric(s$q)
  nonref <- events[events != -1]
  E <- outer(s$e, nonref, "==") * 1
  X <- cbind(1, s$z, E, E * s$z, as.matrix(s[, xs, drop=FALSE]))
  storage.mode(X) <- "double"
  A <- X * sqrt(s$q)
  qrA <- qr(A)
  n <- nrow(X); k <- ncol(X)
  .check(qrA$rank == k, "IDENTIFICATION_FAILURE: singular regression design")
  beta <- as.numeric(qr.coef(qrA, s$y * sqrt(s$q)))
  G <- length(unique(s$c))
  V <- matrix(NA_real_, k, k)
  if (native) {
    .check(n > k, "IDENTIFICATION_FAILURE: no residual degrees of freedom")
    .check(G >= 2, "IDENTIFICATION_FAILURE: at least two original clusters required")
    residual <- s$y - as.numeric(X %*% beta)
    scores <- rowsum(X * (s$q * residual), s$c, reorder=FALSE)
    bread <- solve(crossprod(A))
    V <- bread %*% crossprod(scores) %*% bread * G/(G-1) * (n-1)/(n-k)
    V <- (V + t(V))/2
  }
  ix <- 2 + length(nonref) + seq_along(nonref)
  L <- rep(0, k); L[ix[nonref >= 0]] <- 1/(hi+1)
  estimate <- sum(L * beta)
  se <- sqrt(max(0, as.numeric(t(L) %*% V %*% L)))
  p <- if (native && se > 0) 2*stats::pt(-abs(estimate/se), G-1) else NA_real_
  event <- data.frame(event_time=nonref, estimate=beta[ix], native_se=sqrt(pmax(0, diag(V)[ix])))
  event <- rbind(event, data.frame(event_time=-1, estimate=0, native_se=0))
  list(estimate=estimate, target_label=if (length(xs)) "regression-adjusted Q-weighted post-event contrast" else
         "trimmed aggregate ATT: equal mean over post event times", native_se=se, native_p=p,
       df=if(native) G-1 else NA_real_, native_df=if(native) G-1 else NA_real_,
       status="PASS", event_study=event[order(event$event_time), ], event_covariance=V[ix, ix, drop=FALSE],
       event_covariance_times=nonref, n_obs=nrow(d), n_used=length(unique(s$row)), n_stacked=n,
       n_clusters=G, cohorts=as.numeric(names(shares)), cohort_weights=shares, dropped_cohorts=dropped,
       metadata=list(event_window=c(lo, hi), baseline=-1, data_mode=mode, control_group=control,
                     covariates_used=xs, covariate_timing="contemporaneous additive; no transformation",
                     native_vce="CR1: G/(G-1)*(Nstack-1)/(Nstack-K)",
                     native_reference="Student t, G-1: author Stata lincom / underlying fixest; R scalar helper CI uses normal",
                     oracle_commit="3116918a5773c6415acc9bce6699c5c15b3a1792"),
       stack=s[, c("row", "a", "e", "z", "q", "c")])
}

fit <- function(data, y, time, unit=NULL, gvar=NULL, treatment=NULL, cluster=NULL,
                covariates=NULL, event_window=NULL, data_mode="auto", ...) {
  options <- list(...)
  .check(all(names(options) %in% c("weights", "control_group", "inference", "assignments")) &&
           length(names(options)) == length(options), "UNSUPPORTED options")
  control <- if (is.null(options$control_group)) "both" else options$control_group
  .check(control %in% c("both", "never_treated", "not_yet_treated"), "invalid control_group")
  mode <- if (data_mode == "auto") {if (is.null(unit)) "rcs" else "panel"} else data_mode
  .check(mode %in% c("panel", "rcs"), "invalid data_mode")
  .check(mode != "panel" || !is.null(unit), "panel requires unit")
  if (is.null(cluster)) cluster <- unit
  .check(!is.null(cluster), "original cluster required, particularly for RCS")
  .check(length(event_window) == 2 && all(is.finite(event_window)) &&
           all(event_window == floor(event_window)) && event_window[1] <= -1 && event_window[2] >= 0,
         "event_window=(-pre,post) must include -1 and 0")
  lo <- event_window[1]; hi <- event_window[2]
  data <- as.data.frame(data)
  required <- c(y, time, cluster, unit, gvar, treatment, covariates, options$weights)
  .check(all(required %in% names(data)), "missing requested columns")
  .check(nrow(data) > 0, "empty data")
  d <- data.frame(y=data[[y]], t=data[[time]], c=data[[cluster]], row=seq_len(nrow(data))-1)
  d$u <- if (is.null(unit)) d$row else data[[unit]]
  d$w <- if (is.null(options$weights)) 1 else data[[options$weights]]
  xs <- if (length(covariates)) paste0("x", seq_along(covariates)-1) else character()
  for (i in seq_along(xs)) d[[xs[i]]] <- data[[covariates[i]]]
  .check(!anyNA(d), "missing values in outcome/time/IDs/weights/covariates")
  numeric <- d[, c("y", "t", "w", xs), drop=FALSE]
  .check(all(vapply(numeric, is.numeric, logical(1))) && all(is.finite(as.matrix(numeric))), "nonfinite numeric input")
  .check(all(d$w > 0) && all(d$t == floor(d$t)), "positive weights and integer time required")
  if (mode == "panel") {
    .check(!anyDuplicated(d[, c("u", "t")]), "duplicate unit-time rows")
    .check(all(vapply(split(d$c, d$u), function(v) length(unique(v)) == 1, logical(1))), "cluster must be constant within unit")
  }
  .check(!is.null(gvar) || !is.null(treatment), "gvar or treatment required")
  if (!is.null(gvar)) {
    .check(is.numeric(data[[gvar]]), "cohort dates must be numeric")
    d$g <- data[[gvar]]; d$g[is.na(d$g)] <- Inf
    .check(all(d$g == Inf | (is.finite(d$g) & d$g == floor(d$g))), "invalid cohort dates")
  } else {
    .check(mode == "panel", "UNSUPPORTED: RCS treatment-only input; supply cohort dates")
    D <- data[[treatment]]
    .check(all(D %in% c(0,1)), "treatment must be binary")
    first <- tapply(ifelse(D == 1, d$t, Inf), d$u, min)
    d$g <- unname(first[as.character(d$u)])
    .check(!any(d$g == ave(d$t, d$u, FUN=min)), "cannot infer initially treated adoption; supply gvar")
  }
  if (mode == "panel") .check(all(vapply(split(d$g, d$u), function(v) length(unique(v)) == 1, logical(1))), "cohort must be constant within unit")
  if (!is.null(treatment)) .check(all(data[[treatment]] %in% c(0,1)) &&
    all(data[[treatment]] == as.numeric(d$t >= d$g)), "inconsistent or nonabsorbing treatment history")
  methods <- if (is.null(options$inference)) "native" else options$inference
  .check(all(methods %in% c("point", "native", "jackknife", "ri_beta", "ri_t_native")), "UNSUPPORTED inference")
  .check(!"point" %in% methods || identical(methods,"point"), "point cannot be combined with inference methods")
  ri <- intersect(methods, c("ri_beta", "ri_t_native"))
  result <- .point(d, lo, hi, mode, control, xs, native=!identical(methods,"point"))
  result$events <- result$event_study
  result$metadata$covariates_requested <- covariates
  result$metadata$covariates_used <- covariates
  jack <- function(dd, obs) {
    ids <- unique(dd$c)
    values <- vapply(ids, function(id) {
      rr <- .point(dd[dd$c != id, , drop=FALSE], lo, hi, mode, control, xs, native=FALSE)
      rr$estimate
    }, numeric(1))
    se <- sqrt((length(ids)-1)/length(ids)*sum((values-mean(values))^2))
    .check(se > 1e-14, "IDENTIFICATION_FAILURE: zero jackknife variance")
    list(se=se, p=2*stats::pt(-abs(obs$estimate/se), length(ids)-1), df=length(ids)-1,
         estimates=values, clusters=ids)
  }
  if (any(c("jackknife") %in% methods)) {
    result$jackknife <- jack(d, result)
    result$jk_se <- result$jackknife$se; result$jk_p <- result$jackknife$p
  }
  if (length(ri)) {
    aa <- options$assignments
    .check(is.data.frame(aa) || is.matrix(aa), "UNSUPPORTED: RI requires supplied assignments")
    .check(all(vapply(split(d$g, d$c), function(v) length(unique(v)) == 1, logical(1))), "RI requires homogeneous cohort within cluster")
    ids <- unique(d$c)
    .check(length(colnames(aa)) == length(ids) && setequal(colnames(aa), as.character(ids)), "assignment columns must match original clusters")
    aa <- as.matrix(aa[, as.character(ids), drop=FALSE]); aa[is.na(aa)] <- Inf
    observed <- d$g[match(ids, d$c)]
    .check(nrow(aa) > 0 && !anyDuplicated(as.data.frame(aa)), "assignments must be nonempty and unique")
    for (i in seq_len(nrow(aa))) .check(identical(sort(as.numeric(aa[i,])), sort(observed)) &&
      !identical(as.numeric(aa[i,]), observed), "assignments must preserve cohort counts and exclude observed")
    randomized <- matrix(NA_real_, nrow(aa), length(ri), dimnames=list(NULL, ri))
    for (i in seq_len(nrow(aa))) {
      dd <- d; dd$g <- aa[i, match(dd$c, ids)]
      rr <- .point(dd, lo, hi, mode, control, xs)
      for (m in ri) {
        den <- if (m == "ri_t_native") rr$native_se else 1
        .check(den > 1e-14, "IDENTIFICATION_FAILURE: zero randomized SE")
        randomized[i,m] <- rr$estimate / den
      }
    }
    for (m in ri) {
      den <- if (m == "ri_t_native") result$native_se else 1
      .check(den > 1e-14, "IDENTIFICATION_FAILURE: zero observed SE")
      obs <- result$estimate/den
      tolerance <- 64*.Machine$double.eps*pmax(1,abs(randomized[,m]),abs(obs))
      p <- (1+sum(abs(randomized[,m]) >= abs(obs)-tolerance))/(nrow(aa)+1)
      result[[m]] <- list(p=p, statistic=obs, randomized_statistics=randomized[,m], reps=nrow(aa), null=0)
      result[[paste0(m,"_p")]] <- p
    }
  } else .check(is.null(options$assignments), "assignments supplied without RI method")
  result
}
