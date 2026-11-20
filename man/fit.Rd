\name{fit}
\alias{fit}
\title{Corrective-Weighted Saturated Stacked Event Study}
\description{Native weighted event study with original-cluster inference.}
\usage{fit(data, y, time, unit = NULL, gvar = NULL, treatment = NULL,
cluster = NULL, covariates = NULL, event_window = NULL,
data_mode = "auto", ...)}
\arguments{
\item{data}{A data.frame.}
\item{y,time}{Outcome and integer calendar time column names.}
\item{unit}{Panel unit column name, optional in repeated cross sections.}
\item{gvar}{First adoption date column; NA or Inf denotes never treated.}
\item{treatment}{Alternative absorbing binary treatment column, panel only.}
\item{cluster}{Original cluster column; defaults to unit.}
\item{covariates}{Numeric additive contemporaneous covariate column names.}
\item{event_window}{Required full event window c(-pre,post), including -1 and 0.}
\item{data_mode}{One of auto, panel, rcs.}
\item{...}{weights (column), control_group, inference, assignments. Unknown
options fail. Inference names: point, native, jackknife, ri_beta, ri_t_native,
, . RI requires count-preserving unique supplied assignments.}}
\value{List with estimate, target_label, native_se, native_p, native_df,
event_study, covariance, support metadata, stack and requested inference.}
\details{The scalar averages all post event coefficients relative to -1.
Panel requires balanced within-stack support. RCS weights correct event-cell
composition. Covariates enter without transformation and require additional
identification assumptions. Native variance is original-cluster CR1. Point
mode skips sandwich computation. See package method and inference documents.}
\examples{
d <- expand.grid(time=1:7,id=0:11)
d$g <- ifelse(d$id<4,3,ifelse(d$id<8,5,Inf))
d$y <- .3*d$id+.2*d$time+2*(d$time>=d$g)
fit(d,"y","time",unit="id",gvar="g",event_window=c(-2,1))$estimate
}
