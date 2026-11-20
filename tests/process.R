library(faststackedR)
d <- expand.grid(time=1:7,id=0:11)
d$g <- ifelse(d$id<4,3,ifelse(d$id<8,5,Inf))
d$y <- .3*d$id+.2*d$time+2*(d$time>=d$g)
cl <- parallel::makePSOCKcluster(2)
tryCatch({
  invisible(parallel::clusterCall(cl,function(paths) {
    .libPaths(paths)
    library(faststackedR)
    NULL
  }, .libPaths()))
  results <- parallel::parLapply(cl,1:4,function(i,d) {
    fit(d,"y","time",unit="id",gvar="g",event_window=c(-2,1))$estimate
  },d=d)
  stopifnot(all(abs(unlist(results)-2)<1e-10))
  cat("R_PROCESS_PASS: two independent workers, four deterministic fits\n")
},finally=parallel::stopCluster(cl))
