# R port of core instantaneous-vs-history and IID-vs-AR(1) models from
# Continuous CI models .ipynb. Requires MATLAB/CSV-prepared CI trajectories.
source('_r_helpers.R')
RHO_LOWER <- -.995;RHO_UPPER <- .995
GAMMA_LOWER <- -20;GAMMA_UPPER <- 20
INCLUDE_HISTORY_X_DISPARITY <- TRUE
NUMERIC_VAR_FLOOR <- 1e-12
logistic <- function(x) plogis(x)
logit_clip <- function(p) qlogis(pmin(1-1e-9,pmax(1e-9,p)))
logistic_equilibrium <- function(disparity,alpha,beta) logistic(alpha-beta*disparity)
# Source uses constrained nonlinear least squares [-20,20] for alpha, [0,2] beta.
fit_static_curve <- function(d){
  need_columns(d,c('abs_disparity_deg','mean_p_common'))
  x<-num(d$abs_disparity_deg);y<-num(d$mean_p_common);ok<-is.finite(x)&is.finite(y);x<-x[ok];y<-y[ok]
  if(length(x)<4||length(unique(x))<4)return(NULL)
  a0<-logit_clip(pmin(.99,pmax(.01,y[which.min(x)])))
  fit<-optim(c(a0,.10),function(z)sum((y-logistic_equilibrium(x,z[1],z[2]))^2),method='L-BFGS-B',lower=c(-20,0),upper=c(20,2))
  alpha<-fit$par[1];beta<-fit$par[2];pse<-if(beta>1e-8)alpha/beta else NA_real_;sst<-sum((y-mean(y))^2)
  data.frame(alpha=alpha,beta=beta,static_curve_r2=if(sst>0)1-fit$value/sst else NA_real_,static_curve_pse_deg=if(is.finite(pse)&&pse>=0&&pse<=40)pse else NA_real_,static_curve_pse_extrapolated_deg=pse,pse_in_observed_range=is.finite(pse)&&pse>=0&&pse<=40,static_curve_fit_success=fit$convergence==0)
}
fit_static_by_subject <- function(d,subject='participant_id') {
  need_columns(d,subject)
  do.call(rbind,lapply(split(d,d[[subject]]),function(x){v<-fit_static_curve(x);if(is.null(v))return(NULL);cbind(participant_id=x[[subject]][1],v)}))
}
# Source history-coded mean in logit space: +gamma_dir*h + gamma_int*h*(d-20)/20.
predict_history_mean <- function(disparity,alpha,beta,history_code,gamma_dir=0,gamma_int=0) {
  eta<-alpha-beta*num(disparity);logistic(eta+gamma_dir*history_code+if(INCLUDE_HISTORY_X_DISPARITY)gamma_int*history_code*(num(disparity)-20)/20 else 0)
}
# Profile AR1 innovations likelihood (initial residual omitted, no stationary prior).
ar1_profile_nll <- function(residual_sequences,rho){
  innovations<-unlist(lapply(residual_sequences,function(e){e<-num(e);e<-e[is.finite(e)];if(length(e)<2)return(numeric(0));e[-1]-rho*e[-length(e)]}),use.names=FALSE)
  n<-length(innovations);if(!n)return(Inf)
  var<-max(mean(innovations^2),NUMERIC_VAR_FLOOR)
  .5*n*(log(2*pi*var)+1)
}
fit_ar1_error <- function(residual_sequences){
  opt<-optimize(function(rho)ar1_profile_nll(residual_sequences,rho),c(RHO_LOWER,RHO_UPPER),tol=1e-7)
  n<-sum(vapply(residual_sequences,function(x)max(0,sum(is.finite(x))-1),0L))
  list(rho=opt$minimum,sigma=sqrt(exp(2*opt$objective/n-1)/(2*pi)),loglik=-opt$objective,n_innovations=n,mean_logscore=-opt$objective/n)
}
fit_iid_error <- function(residual_sequences){e<-unlist(residual_sequences,use.names=FALSE);e<-e[is.finite(e)];n<-length(e)
  variance<-max(mean(e^2),NUMERIC_VAR_FLOOR)
  list(sigma=sqrt(variance),loglik=-.5*n*(log(2*pi*variance)+1),n=n)
}
# Trial dataset columns: participant_id, trial_id, analysis_group, time_s,
# disparity_deg, p_common, alpha, beta. The notebook's fixed-offset fits
# supply alpha/beta before fitting sweep history models.
fit_participant_history <- function(d,subject_id=NULL,ar1=TRUE) {
  need_columns(d,c('trial_id','disparity_deg','p_common','alpha','beta','analysis_group'))
  d<-d[d$analysis_group%in%c('sweep_diverge','sweep_converge'),,drop=FALSE]
  parts<-split(d,d$trial_id)
  dlist<-lapply(parts,function(x){x<-x[order(num(x$time_s)),];h<-if(x$analysis_group[1]=='sweep_diverge') .5 else -.5
    list(disparity=num(x$disparity_deg),p=num(x$p_common),alpha=num(x$alpha[1]),beta=num(x$beta[1]),h=h)})
  residuals<-function(g){lapply(dlist,function(x){pred<-predict_history_mean(x$disparity,x$alpha,x$beta,x$h,g[1],g[2]);x$p-pred})}
  objective<-function(g)if(ar1)ar1_profile_nll(residuals(g),fit_ar1_error(residuals(g))$rho) else -fit_iid_error(residuals(g))$loglik
  fit<-optim(c(0,0),objective,method='L-BFGS-B',lower=c(GAMMA_LOWER,GAMMA_LOWER),upper=c(GAMMA_UPPER,GAMMA_UPPER))
  err<-if(ar1)fit_ar1_error(residuals(fit$par)) else fit_iid_error(residuals(fit$par))
  data.frame(participant_id=if(is.null(subject_id))NA_character_ else subject_id,gamma_dir=fit$par[1],gamma_int=fit$par[2],rho=if(ar1)err$rho else NA_real_,sigma=err$sigma,train_loglik=err$loglik,success=fit$convergence==0)
}
run <- function(curve_csv='CI_subject_static_curve.csv',sweeps_csv='CI_prepared_sweeps.csv') {
  curves<-read_input(curve_csv);fits<-fit_static_by_subject(curves)
  save_results(fits,'CI_static_equilibrium_model_parameters_R.csv')
  d<-read_input(sweeps_csv)
  if(!all(c('alpha','beta')%in%names(d))){d<-merge(d,fits[,c('participant_id','alpha','beta')],by='participant_id',all.x=TRUE)}
  out<-do.call(rbind,lapply(split(d,d$participant_id),function(x)fit_participant_history(x,subject_id=x$participant_id[1])))
  save_results(out,'CI_history_AR1_R.csv');invisible(list(static=fits,history=out))
}
if(sys.nframe()==0)run()
