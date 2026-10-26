# R port of Continuous CI v2 trial trajectories: P(common), confidence,
# time courses, fixed-offset PSE, and within-person sweep hysteresis.
# The notebook's other cluster-permutation, breakpoint and sustained-event
# sensitivity analyses require separate validation and are not included here.
source('_r_helpers.R')
PCOMMON_FROM_X_MODE <- 'neg1_to_1'
CONFIDENCE_FROM_Y_MODE <- 'neg1_to_1'
N_PERM <- 10000L;RANDOM_SEED <- 12345L
clip01 <- function(x)pmin(1,pmax(0,x))
convert_x_to_pcommon <- function(x,mode=PCOMMON_FROM_X_MODE){x<-num(x);switch(mode,raw=x,neg1_to_1=clip01((x+1)/2),flip_neg1_to_1=1-clip01((x+1)/2),zero_to_one=clip01(x),flip_zero_to_one=1-clip01(x),stop('Unknown conversion mode'))}
convert_y_to_confidence <- function(y,mode=CONFIDENCE_FROM_Y_MODE){y<-num(y);switch(mode,raw=y,neg1_to_1=clip01((y+1)/2),zero_to_one=clip01(y),stop('Unknown conversion mode'))}
classify_ci_group <- function(s){s<-tolower(as.character(s));ifelse(s=='congruent','congruent',ifelse(grepl('^spatial_10',s),'fixed_10',ifelse(grepl('^spatial_20',s),'fixed_20',ifelse(grepl('^spatial_30',s),'fixed_30',ifelse(grepl('^spatial_40',s),'fixed_40',ifelse(grepl('^sweep_converge',s),'sweep_converge',ifelse(grepl('^sweep_diverge',s),'sweep_diverge',NA_character_)))))))}
read_ci_samples <- function(filename='CI_trial_samples.csv') {
  d<-read_input(filename)
  # Columns correspond to the original task: raw x = same/different continuum,
  # raw y = confidence continuum; no arbitrary reflection of responses.
  need_columns(d,c('participant_id','trial_id','plot_condition','time_s','mouse_x','mouse_y'))
  d$p_common<-convert_x_to_pcommon(d$mouse_x)
  d$confidence<-convert_y_to_confidence(d$mouse_y)
  d$analysis_group<-classify_ci_group(d$plot_condition)
  d
}
subject_first_condition_means <- function(d) {
  k<-c('participant_id','trial_id','analysis_group')
  tr<-aggregate(d[c('p_common','confidence')],d[k],mean_na)
  z<-aggregate(tr[c('p_common','confidence')],tr[c('participant_id','analysis_group')],mean_na)
  list(trial=tr,participant=z)
}
condition_summary <- function(d) {
  x<-split(d,d$analysis_group)
  do.call(rbind,lapply(x,function(z){
    cbind(condition=z$analysis_group[1],outcome=c('p_common','confidence'),do.call(rbind,lapply(z[c('p_common','confidence')],function(x)as.data.frame(as.list(ci_t(x))))))
  }))
}
# Binned subject-first trajectories. Each participant contributes equally.
trajectory_means <- function(d,bin_sec=.5) {
  d$time_bin<-floor(num(d$time_s)/bin_sec)*bin_sec
  tr<-aggregate(d[c('p_common','confidence')],d[c('participant_id','trial_id','analysis_group','time_bin')],mean_na)
  sub<-aggregate(tr[c('p_common','confidence')],tr[c('participant_id','analysis_group','time_bin')],mean_na)
  group<-aggregate(sub[c('p_common','confidence')],sub[c('analysis_group','time_bin')],mean_na)
  list(participants=sub,group=group)
}
# Static fixed-offset mean same-source function: plogis(alpha - beta*disparity).
fit_fixed_offset_pse <- function(subj,offsets=c(10,20,30,40)) {
  x<-subj[subj$analysis_group%in%paste0('fixed_',offsets),,drop=FALSE]
  sp<-split(x,x$participant_id)
  out<-lapply(sp,function(d){
    disp<-num(sub('fixed_','',as.character(d$analysis_group)));prob<-num(d$p_common)
    if(length(unique(disp))<4)return(NULL)
    fun<-function(z)sum((prob-plogis(z[1]-z[2]*disp))^2)
    o<-optim(c(qlogis(pmax(.01,pmin(.99,prob[which.min(disp)]))),.10),fun,method='L-BFGS-B',lower=c(-20,0),upper=c(20,2))
    a<-o$par[1];b<-o$par[2];pse<-if(b>1e-8)a/b else NA_real_
    data.frame(participant_id=d$participant_id[1],alpha=a,beta=b,PSE_deg=if(is.finite(pse)&&pse>=0&&pse<=40)pse else NA_real_,PSE_extrapolated_deg=pse,fit_success=o$convergence==0)
  });do.call(rbind,out)
}
# Uses matched instantaneous disparity bins rather than trial timing to avoid
# labeling simple direction effects as hysteresis.
sweep_hysteresis <- function(d,bin_width_deg=5) {
  need_columns(d,'abs_disparity_deg')
  d<-d[d$analysis_group%in%c('sweep_diverge','sweep_converge') & is.finite(d$abs_disparity_deg),]
  d$disparity_bin<-floor(d$abs_disparity_deg/bin_width_deg)*bin_width_deg
  tr<-aggregate(d[c('p_common','confidence')],d[c('participant_id','trial_id','analysis_group','disparity_bin')],mean_na)
  sub<-aggregate(tr[c('p_common','confidence')],tr[c('participant_id','analysis_group','disparity_bin')],mean_na)
  w<-reshape(sub,idvar=c('participant_id','disparity_bin'),timevar='analysis_group',direction='wide')
  for(o in c('p_common','confidence')){
    need_columns(w,paste0(o,c('.sweep_diverge','.sweep_converge')))
    w[[paste0(o,'_diff')]]<-w[[paste0(o,'.sweep_diverge')]]-w[[paste0(o,'.sweep_converge')]]
  }
  bysub<-aggregate(w[c('p_common_diff','confidence_diff')],list(participant_id=w$participant_id),mean_na)
  summaries<-rbind(data.frame(measure='p_common',as.data.frame(as.list(ci_t(bysub$p_common_diff))),p_perm=signflip_test(bysub$p_common_diff,N_PERM,RANDOM_SEED)),data.frame(measure='confidence',as.data.frame(as.list(ci_t(bysub$confidence_diff))),p_perm=signflip_test(bysub$confidence_diff,N_PERM,RANDOM_SEED)))
  list(bins=w,participants=bysub,summary=summaries)
}
# Replicates 1-constant, 2-linear and 4-parameter continuous piecewise BIC.
fit_ci_trajectory_breakpoint <- function(d,time_col='time_s',metric='p_common',min_segment_s=2) segmented_bic(num(d[[time_col]]),num(d[[metric]]),min_segment_s)
run <- function() {
  d<-read_ci_samples();means<-subject_first_condition_means(d)
  save_results(means$trial,'CIv2_trial_means_R.csv');save_results(means$participant,'CIv2_subject_means_R.csv')
  save_results(condition_summary(means$participant),'CIv2_condition_summary_R.csv')
  curve<-fit_fixed_offset_pse(means$participant)
  if(!is.null(curve))save_results(curve,'CIv2_subject_PSE_R.csv')
  if('abs_disparity_deg'%in%names(d)){
    h<-sweep_hysteresis(d);save_results(h$bins,'CIv2_hysteresis_bins_R.csv')
    save_results(h$participants,'CIv2_hysteresis_subject_R.csv');save_results(h$summary,'CIv2_hysteresis_group_R.csv')
  }
  t<-trajectory_means(d);save_results(t$participants,'CIv2_subject_trajectories_R.csv');save_results(t$group,'CIv2_group_trajectories_R.csv')
  invisible(list(trials=d,means=means,pse=curve))
}
if(sys.nframe()==0)run()
