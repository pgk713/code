# R port of the Continuous Spatial Localization v3 primary lag/bias/capture metrics.
# The 145-cell Python notebook additionally implements many sensitivity,
# breakpoint, reliability, and motor-residualized analyses not yet ported here.
source('_r_helpers.R')
TRIAL_DURATION_SEC <- 18; RESAMPLE_HZ <- 20
TIME_GRID <- seq(0,TRIAL_DURATION_SEC,by=1/RESAMPLE_HZ)
MAX_LAG_SEC <- 4.0; N_PERM <- 10000L; SEED <- 12345L
# Trial CSV: one row per recorded sample with participant/trial/condition/t,
# mouse_x,blob_x,audio_x. Values must match notebook coordinate units.
INPUT_CSV <- 'SL_trial_samples.csv'
trial_metrics <- function(d,participant='participant',trial='trial',condition='condition') {
  needed<-c(participant,trial,condition,'t','mouse_x','blob_x','audio_x');need_columns(d,needed)
  parts<-split(d,interaction(d[[participant]],d[[trial]],drop=TRUE))
  do.call(rbind,lapply(parts,function(x){
    z<-resample_trace(x,TIME_GRID);if(is.null(z))return(NULL)
    lag_a<-lag_scanned_corr(z$mouse_x,z$audio_x,RESAMPLE_HZ,MAX_LAG_SEC)
    lag_v<-lag_scanned_corr(z$mouse_x,z$blob_x,RESAMPLE_HZ,MAX_LAG_SEC)
    mouse_a<-lag_correct(z$mouse_x,z$t,lag_a$lag_sec)
    mouse_v<-lag_correct(z$mouse_x,z$t,lag_v$lag_sec)
    data.frame(participant=x[[participant]][1],trial=x[[trial]][1],condition=x[[condition]][1],
      PeakA=lag_a$peak_corr,LagA_sec=lag_a$lag_sec,
      PeakV=lag_v$peak_corr,LagV_sec=lag_v$lag_sec,
      RMSE_A=rmse(mouse_a,z$audio_x),RMSE_V=rmse(mouse_v,z$blob_x),
      Bias_A=mean_na(mouse_a-z$audio_x),Bias_V=mean_na(mouse_v-z$blob_x),
      VisualCapture=capture_weight(z$mouse_x,z$blob_x,z$audio_x),
      VisualCaptureLagAdjusted=capture_weight(mouse_a,z$blob_x,z$audio_x))
  }))
}
# Matched subject-first differences: avoid pseudoreplication across trial samples.
subject_condition_summary <- function(m) {
  needed<-c('PeakA','LagA_sec','PeakV','LagV_sec','RMSE_A','RMSE_V','Bias_A','Bias_V','VisualCapture','VisualCaptureLagAdjusted')
  groups<-split(m,interaction(m$participant,m$condition,drop=TRUE))
  do.call(rbind,lapply(groups,function(x){
    rec<-data.frame(participant=x$participant[1],condition=x$condition[1],n_trials=nrow(x))
    for(k in needed)rec[[k]]<-mean_na(num(x[[k]]));rec
  }))
}
paired_lag_difference <- function(df){
  d<-df[is.finite(df$LagA_sec)&is.finite(df$LagV_sec),]
  paired_test(d$LagA_sec,d$LagV_sec)
}
# Classification follows fixed spatial offsets and sweep condition naming.
classify_spatial <- function(condition){s<-tolower(as.character(condition));ifelse(grepl('sweep.*diverge|diverge',s),'sweep_diverge',ifelse(grepl('sweep.*converge|converge',s),'sweep_converge',ifelse(grepl('congruent',s),'congruent',ifelse(grepl('10',s),'fixed_10',ifelse(grepl('20',s),'fixed_20',ifelse(grepl('30',s),'fixed_30',ifelse(grepl('40',s),'fixed_40','other')))))))}
# Optional precomputed trial capture trajectory: segmented BIC search.
fit_capture_breakpoint <- function(windows_csv,metric='visual_capture',min_segment_sec=2){
  d<-read_input(windows_csv);need_columns(d,c('participant','trial','window_center_sec',metric))
  x<-split(d,interaction(d$participant,d$trial,drop=TRUE))
  z<-lapply(x,function(sub){fit<-segmented_bic(num(sub$window_center_sec),num(sub[[metric]]),min_segment_sec)
    if(is.null(fit))return(NULL)
    data.frame(participant=sub$participant[1],trial=sub$trial[1],model=fit$model,BIC=fit$bic,breakpoint_s=fit$breakpoint,pre_slope=fit$pre_slope,post_slope=fit$post_slope)
  });do.call(rbind,z)
}
# Signed hysteresis needs *matched disparity bins* within participant.
hysteresis_at_matched_disparity <- function(d,participant='participant',direction='sweep_direction',disparity='disparity_deg',capture='VisualCapture',bin_width=5){
  need_columns(d,c(participant,direction,disparity,capture))
  d$disparity_bin<-floor(num(d[[disparity]])/bin_width)*bin_width
  g<-aggregate(num(d[[capture]]),d[c(participant,direction,'disparity_bin')],mean_na);names(g)[4]<-'capture'
  names(g)[1:2]<-c('participant','direction')
  w<-reshape(g,idvar=c('participant','disparity_bin'),timevar='direction',direction='wide')
  need_columns(w,c('capture.sweep_diverge','capture.sweep_converge'))
  w$diverge_minus_converge<-w$capture.sweep_diverge-w$capture.sweep_converge
  bysub<-aggregate(w$diverge_minus_converge,list(participant=w$participant),mean_na);names(bysub)[2]<-'hysteresis'
  list(bins=w,subject=bysub,group_test=ci_t(bysub$hysteresis),permutation_p=signflip_test(bysub$hysteresis,N_PERM,SEED))
}
run <- function() {
  d<-read_input(INPUT_CSV);metrics<-trial_metrics(d);save_results(metrics,'SL_trial_metrics_R.csv')
  s<-subject_condition_summary(metrics);save_results(s,'SL_subject_condition_R.csv')
  lag<-paired_lag_difference(s);save_results(lag,'SL_auditory_visual_lag_paired_R.csv')
  invisible(list(trials=metrics,summary=s,lag_test=lag))
}
if(sys.nframe()==0)run()
