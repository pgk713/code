# Port of the fNIRS ISPS EI notebook's analytic core AFTER MNE preprocessing.
# Inputs must be pre-extracted MATLAB/CSV HbO ROI time series, not raw SNIRF.
source('_r_helpers.R')
FS <- 1.5; SMOOTH_SEC <- 7.0; ISPS_METRIC <- 'kuramoto'; N_PERM <- 10000L
# Expected long data: participant,trigger,ROI,time,hbo; matched time samples.
prepare_roi_phase <- function(d,subject='participant',trigger='trigger',roi='ROI',time='time',value='hbo') {
  need_columns(d,c(subject,trigger,roi,time,value))
  split_key<-interaction(d[[trigger]],d[[roi]],drop=TRUE)
  lapply(split(d,split_key),function(sub){
    parts<-split(sub,sub[[subject]])
    series<-lapply(parts,function(x){x<-x[order(num(x[[time]])),];num(x[[value]])})
    len<-min(lengths(series));if(len<3)return(NULL)
    phases<-t(vapply(series,function(x)phase_hilbert(x[seq_len(len)]),numeric(len)))
    rownames(phases)<-names(series)
    list(trigger=sub[[trigger]][1],ROI=sub[[roi]][1],phase=phases,time=sort(unique(num(sub[[time]])))[seq_len(len)])
  })
}
# Within-subject circular temporal shifts preserve each individual's structure
# while disrupting cross-subject phase alignment (as in source permutations).
permutation_isps <- function(P,n_perm=N_PERM,seed=12345L) {
  set.seed(seed);P<-as.matrix(P); obs<-isps_kuramoto(P,FS,SMOOTH_SEC,ISPS_METRIC)
  null_mean<-numeric(n_perm)
  for(b in seq_len(n_perm)){
    S<-t(vapply(seq_len(nrow(P)),function(i){k<-sample.int(ncol(P),1)-1L;P[i,((seq_len(ncol(P))+k-1L)%%ncol(P))+1L]},numeric(ncol(P))))
    null_mean[b]<-mean_na(isps_kuramoto(S,FS,SMOOTH_SEC,ISPS_METRIC))
  }
  c(observed_mean=mean_na(obs),p=(1+sum(null_mean>=mean_na(obs)))/(n_perm+1L),null_mean=mean(null_mean))
}
roi_isps_results <- function(d,n_perm=N_PERM) {
  runs<-Filter(Negate(is.null),prepare_roi_phase(d))
  out<-lapply(seq_along(runs),function(i){r<-runs[[i]];p<-r$phase
    trace<-isps_kuramoto(p,FS,SMOOTH_SEC,ISPS_METRIC)
    perm<-permutation_isps(p,n_perm,seed=12345L+i)
    data.frame(trigger=r$trigger,ROI=r$ROI,n_subjects=nrow(p),peak_isps=max(trace,na.rm=TRUE),peak_time=r$time[which.max(trace)],t(perm))
  });z<-do.call(rbind,out);if(is.null(z))return(NULL);z$q_bh<-bh_adjust(z$p);z
}
# Phase-based results assume meaningful analytic phase; compare the source
# notebook's bandpass filtering and artifact exclusions before using in a paper.
run <- function(filename='fnirs_isps_roi_timeseries.csv') {
  d<-read_input(filename);r<-roi_isps_results(d)
  if(!is.null(r))save_results(r,'fnirs_ISPS_ROI_results.csv')
  r
}
if(sys.nframe()==0)run()
