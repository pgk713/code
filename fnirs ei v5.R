# R port of downstream intersubject correlation, IS-RSA, regression and LASSO
# after fNIRS MNE signal extraction. MNE-specific preprocessing intentionally omitted.
source('_r_helpers.R')
# A time-aligned ROI HbO file must contain subject,video,ROI,time,HbO.
# Input series must have the same preprocessing/lag correction as the notebook.
ISC_TIME_SERIES <- 'EI_ROI_HbO_timeseries.csv'
SCORES <- 'EI_subject_survey_scores.csv'
N_PERM <- 10000L
compute_pairwise_isc <- function(d,subject='subject',video='video',roi='ROI',time='time',signal='HbO') {
  need_columns(d,c(subject,video,roi,time,signal))
  keys<-split(d,interaction(d[[video]],d[[roi]],drop=TRUE))
  result<-lapply(keys,function(x){
    subjects<-unique(as.character(x[[subject]]));if(length(subjects)<3)return(NULL)
    traces<-lapply(subjects,function(id){s<-x[as.character(x[[subject]])==id,];s<-s[order(num(s[[time]])),];num(s[[signal]])})
    minlen<-min(lengths(traces));if(minlen<4)return(NULL)
    # Explicitly assume synchronized time bins; DO NOT correlate unequal clocks.
    X<-do.call(cbind,lapply(traces,function(t)t[seq_len(minlen)]));colnames(X)<-subjects
    R<-cor(X,use='pairwise.complete.obs');diag(R)<-1
    list(video=x[[video]][1],ROI=x[[roi]][1],R=R,mean_isc=mean(R[lower.tri(R)],na.rm=TRUE),subjects=subjects)
  });Filter(Negate(is.null),result)
}
pairwise_isc_table <- function(isc_list) {
  do.call(rbind,lapply(isc_list,function(x){ids<-x$subjects;ij<-which(lower.tri(x$R),arr.ind=TRUE)
    data.frame(video=x$video,ROI=x$ROI,subject1=ids[ij[,1]],subject2=ids[ij[,2]],ISC=x$R[ij])
  }))
}
# Brain similarity matrix averaged over ROI/video and score-distance RDM.
score_brain_rsa <- function(isc_edges,scores,score_col,subject_col='subject',n_perm=N_PERM) {
  need_columns(isc_edges,c('subject1','subject2','ISC'))
  need_columns(scores,c(subject_col,score_col))
  x<-aggregate(isc_edges$ISC,isc_edges[c('subject1','subject2')],mean_na)
  ids<-unique(as.character(scores[[subject_col]]));n<-length(ids)
  M<-matrix(NA_real_,n,n);diag(M)<-1;rownames(M)<-colnames(M)<-ids
  for(i in seq_len(nrow(x))){a<-match(x$subject1[i],ids);b<-match(x$subject2[i],ids);if(!is.na(a)&&!is.na(b))M[a,b]<-M[b,a]<-x$x[i]}
  if(anyNA(M[lower.tri(M)]))stop('Missing subject pairs in brain matrix')
  p<-num(scores[[score_col]][match(ids,as.character(scores[[subject_col]]))]);if(anyNA(p))stop('Missing score')
  # Source uses score distances vs (1 - brain similarity), Spearman and 10k permutations.
  res<-permute_rdm(M,p,n_perm=n_perm)
  data.frame(measure=score_col,observed_r=res['observed'],permutation_p=res['p_two_sided'],n_subjects=n)
}
# Participant-level aggregated ROI feature models: surveys only vs surveys+ROI.
fit_hierarchical <- function(df,outcome='Performance',survey=c('AQ','EQ','EmQ'),roi=NULL) {
  need_columns(df,c(outcome,survey))
  if(is.null(roi))roi<-grep('^(Left |Right |L_|R_)',names(df),value=TRUE)
  if(!length(roi))stop('Specify ROI feature columns')
  d<-df[,c(outcome,survey,roi),drop=FALSE];d<-d[complete.cases(d),]
  base<-lm(reformulate(survey,outcome),data=d)
  full<-lm(reformulate(c(survey,roi),outcome),data=d)
  list(baseline=summary(base),full=summary(full),comparison=anova(base,full))
}
lasso_nested_cv <- function(X,y,nfolds=5L,seed=12345L) {
  if(!requireNamespace('glmnet',quietly=TRUE))stop('install.packages("glmnet")')
  X<-as.matrix(X); y<-num(y)
  set.seed(seed); cv<-glmnet::cv.glmnet(X,y,alpha=1,nfolds=min(nfolds,length(y)),standardize=TRUE,intercept=TRUE)
  list(model=cv,selected_lambda=cv$lambda.min,coef=coef(cv,s='lambda.min'))
}
run <- function() {
  d<-read_input(ISC_TIME_SERIES);r<-compute_pairwise_isc(d)
  edges<-pairwise_isc_table(r);save_results(edges,'EI_pairwise_ISC.csv')
  summary<-do.call(rbind,lapply(r,function(z)data.frame(video=z$video,ROI=z$ROI,mean_pairwise_ISC=z$mean_isc,n_subjects=length(z$subjects))))
  save_results(summary,'EI_ISC_ROI_video_summary.csv')
  if(file.exists(file.path(DATA_DIR,SCORES))){
    s<-read_input(SCORES)
    for(metric in intersect(c('AQ','EQ','EmQ','Performance'),names(s))) {
      rr<-score_brain_rsa(edges,s,metric)
      save_results(rr,paste0('EI_IS_RSA_',metric,'.csv'))
    }
  }
  invisible(list(edges=edges,summary=summary))
}
if(sys.nframe()==0)run()
