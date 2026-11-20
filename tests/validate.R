if (file.exists("R/fit.R")) source("R/fit.R") else library(faststackedR)

make_data <- function(rcs=FALSE) {

  d <- expand.grid(time=1:7, id=0:11)

  if (rcs) d <- d[rep(seq_len(nrow(d)), 1+(d$id+d$time)%%3), ]

  d$j <- ave(d$id, interaction(d$id,d$time), FUN=seq_along)-1

  d$g <- ifelse(d$id<4,3,ifelse(d$id<8,5,Inf))

  d$x <- sin(d$id*2.3+d$time*.8+d$j)

  d$x0 <- d$id%%2

  d$w <- 1+d$id/20

  d$y <- d$id*.3+d$time*.2+(d$time>=d$g)*(1+d$id*.2)+.3*d$x+.2*cos(d$id*d$time+d$j)

  d

}

run <- function(d, ...) fit(d, "y", "time", unit="id", gvar="g", event_window=c(-2,1), ...)

d <- make_data()

one <- d; one$cluster <- 1

point <- run(one,cluster="cluster",inference="point")

stopifnot(is.na(point$native_se),is.na(point$native_df),abs(point$estimate-run(d)$estimate)<1e-12)

r <- run(d,inference="jackknife")

vv <- sapply(0:11,function(i) run(d[d$id!=i,])$estimate)

stopifnot(max(abs(vv-r$jackknife$estimates))<1e-10,

          abs(r$jk_se-sqrt(11/12*sum((vv-mean(vv))^2)))<1e-10)

aa <- t(sapply(1:3,function(k) d$g[match((0:11-k)%%12,d$id)]))

colnames(aa) <- as.character(0:11)

ri <- run(d,inference=c("jackknife","ri_beta","ri_t_native"),assignments=aa)

stopifnot(ri$ri_t_native_p>0,ri$ri_t_native_p<=1)

hom <- d; hom$y <- .4*d$id+.3*d$time+2*(d$time>=d$g)

stopifnot(abs(run(hom)$estimate-2)<1e-10)

expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr),error=identity),"error"))

expect_error(run(d,fe="interacted"))

expect_error(run(d[-2,]))

expect_error(run(d,inference="ri_beta",assignments=rbind(aa,aa)))

single <- d[d$id==0|d$id>=4,]

jk <- run(single,inference="jackknife")$jackknife

physical <- vapply(jk$clusters,function(i) run(single[single$id!=i,])$estimate,numeric(1))

stopifnot(max(abs(jk$estimates-physical))<1e-10)

cat("R native checks passed: 8 assertions/groups\n")



oracle <- Sys.getenv("STACKED_ORACLE")

if (nzchar(oracle)) {

  .libPaths(c(".test-library", .libPaths()))

  stopifnot(requireNamespace("fixest",quietly=TRUE),requireNamespace("data.table",quietly=TRUE))

  fixest::setFixest_nthreads(1)

  env <- new.env(parent=globalenv())

  for (f in list.files(file.path(oracle,"R","stacked","R"),"[.]R$",full.names=TRUE)) sys.source(f,env)

  results <- list()

  for (rcs in c(FALSE,TRUE)) for (weighted in c(FALSE,TRUE)) for (adjust in c(FALSE,TRUE)) {

    dd <- make_data(rcs)

    mode <- if(rcs) "rcs" else "panel"

    cv <- if(adjust) c("x","x0") else NULL

    w <- if(weighted) "w" else NULL

    ours <- run(dd,data_mode=mode,weights=w,covariates=cv)

    refdata <- dd; refdata$g[is.infinite(refdata$g)] <- NA_real_

    ss <- env$build_stack(refdata,"time","id","g",kappa_pre=2,kappa_post=1,

                          weight_var=w,data_type=if(rcs) "repeated_cross_section" else "panel")

    rr <- env$stackreg(ss,"y",cluster_var="id",covariates=cv)

    avg <- attr(rr,"avg_post_att")

    native_df <- fixest::degrees_freedom(rr,type="t")

    ct <- fixest::coeftable(rr)

    stopifnot(native_df==ours$native_df,

      max(abs(ct[,"Pr(>|t|)"]-2*stats::pt(-abs(ct[,"t value"]),native_df)))<1e-12,

      abs(avg$ci_upper-avg$estimate-stats::qnorm(.975)*avg$se)<1e-12)

    err <- c(estimate=abs(ours$estimate-avg$estimate),se=abs(ours$native_se-avg$se),

             weights=max(abs(ours$stack$q-ss$q_weight)))

    stopifnot(max(err)<1e-9)

    results[[length(results)+1]] <- data.frame(mode=mode,weighted=weighted,adjusted=adjust,

      estimate=ours$estimate,se=ours$native_se,native_p=ours$native_p,native_df=native_df,

      estimate_error=err[1],se_error=err[2],q_error=err[3])

    name <- paste(mode,weighted,adjust,sep="_")

    write.csv(dd,file.path("tests",paste0(name,".csv")),row.names=FALSE)

  }

  write.csv(do.call(rbind,results),"tests/oracle_results.csv",row.names=FALSE)

  print(do.call(rbind,results))

  # RI outputs provide an independent native cross-language test fixture.

  write.csv(data.frame(method=c("native","jackknife","ri_beta","ri_t_native"),

    value=c(ri$native_se,ri$jk_se,ri$ri_beta_p,ri$ri_t_native_p)),

    "tests/inference_results.csv",row.names=FALSE)

  official <- function(dd) {

    dd$g[is.infinite(dd$g)] <- NA_real_

    st <- env$build_stack(dd,"time","id","g",kappa_pre=2,kappa_post=1)

    attr(env$stackreg(st,"y",cluster_var="id"),"avg_post_att")

  }

  official_jack <- function(dd) {

    ids <- unique(dd$id)

    vals <- vapply(ids,function(id) official(dd[dd$id!=id,])$estimate,numeric(1))

    sqrt((length(ids)-1)/length(ids)*sum((vals-mean(vals))^2))

  }

  oo <- official(d); oj <- official_jack(d)

  stopifnot(abs(oj-ri$jk_se)<1e-10)

  oracle_stats <- t(vapply(seq_len(nrow(aa)),function(i) {

    dd <- d; dd$g <- aa[i, match(dd$id,0:11)]

    oo <- official(dd)

    c(oo$estimate,oo$estimate/oo$se)

  }, numeric(2)))

  differences <- c(jackknife_se=abs(oj-ri$jk_se),

    ri_beta=max(abs(oracle_stats[,1]-ri$ri_beta$randomized_statistics)),

    ri_t_native=max(abs(oracle_stats[,2]-ri$ri_t_native$randomized_statistics)))

  stopifnot(max(differences)<1e-9)

  write.csv(data.frame(check=names(differences),max_error=differences),"tests/oracle_inference_errors.csv",row.names=FALSE)

  print(differences)

  probs <- expand.grid(statistic=c(0,.1,1,2,5,10,100),df=c(1,2,5,10,30,100))

  probs$p <- 2*stats::pt(-abs(probs$statistic),probs$df)

  write.csv(probs,"tests/probability_oracle.csv",row.names=FALSE)

  cat("Official frozen R oracle: 8 scenarios x point/CR1/Q checks passed\n")

}

