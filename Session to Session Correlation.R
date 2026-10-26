# R translation of MATLAB-exported data analyses in Session to Session Correlation.ipynb
# MNE/SNIRF loading and preprocessing deliberately omitted.
source('_r_helpers.R')
filter_glm <- function(d, q_threshold=.06, exclude_pairs=TRUE) {
  need_columns(d,c('type','detector','source','q','cond','tstat'))
  d<-d[tolower(trimws(d$type))=='hbo' & !(num(d$detector)%in%16:23) & num(d$q)<=q_threshold,,drop=FALSE]
  if(exclude_pairs) d<-d[!(num(d$source)==14 & num(d$detector)==13) & !(num(d$source)==16 & num(d$detector)==12),,drop=FALSE]
  d
}
reliability_metrics <- function(x,y) {
  v<-cor_test(x,y); sx<-sd_na(x)^2;sy<-sd_na(y)^2
  r<-unname(v['r'])
  c(pearson_r=r,spearman_brown=if(is.finite(r)) if(r<1) (2*r)/(1+r) else 1 else NA_real_,variance_stability=if(is.finite(sx+sy) && max(sx,sy)>0) min(sx,sy)/max(sx,sy) else 0,n_channels=unname(v['n']),p_value=unname(v['p']))
}
split_half_group <- function(file_a='1st_half_GLM_GROUP_ROI_AGG_tvr_clean.csv', file_b='2nd_half_GLM_GROUP_ROI_AGG_tvr_clean.csv') {
  a<-filter_glm(read_input(file_a));b<-filter_glm(read_input(file_b))
  out<-lapply(unique(c(a$cond,b$cond)), function(cond) {
    aa<-a[a$cond==cond,,drop=FALSE];bb<-b[b$cond==cond,,drop=FALSE]
    m<-merge(aa,bb,by=c('source','detector'),suffixes=c('_s1','_s2'))
    if(nrow(m)<3)return(NULL)
    data.frame(condition=cond,t(reliability_metrics(m$tstat_s1,m$tstat_s2)))
  })
  do.call(rbind,out)
}
split_half_subject <- function(file_even='even_subject_level_contrasts.csv',file_odd='odd_subject_level_contrasts.csv') {
  a<-filter_glm(read_input(file_even),exclude_pairs=FALSE);b<-filter_glm(read_input(file_odd),exclude_pairs=FALSE)
  need_columns(a,'subject');need_columns(b,'subject')
  keys<-unique(rbind(a[c('subject','cond')],b[c('subject','cond')]))
  result<-lapply(seq_len(nrow(keys)),function(i){
    subj<-keys$subject[i];cond<-keys$cond[i]
    e<-a[a$subject==subj&a$cond==cond,,drop=FALSE];o<-b[b$subject==subj&b$cond==cond,,drop=FALSE]
    m<-merge(e,o,by=c('source','detector'),suffixes=c('_even','_odd'))
    if(nrow(m)<3)return(NULL)
    data.frame(subject=subj,condition=cond,t(reliability_metrics(m$tstat_even,m$tstat_odd)))
  }); do.call(rbind,result)
}
# Two-way random effects single-measure absolute-agreement ICC(2,1).
# Subjects = targets (rows), two sessions = raters (columns).
icc_2_1 <- function(Y) {
  Y<-as.matrix(Y);Y<-Y[complete.cases(Y),,drop=FALSE]; n<-nrow(Y);k<-ncol(Y)
  if(n<3 || k<2) return(NA_real_)
  grand<-mean(Y); rowm<-rowMeans(Y);colm<-colMeans(Y)
  msr<-k*sum((rowm-grand)^2)/(n-1)
  msc<-n*sum((colm-grand)^2)/(k-1)
  mse<-sum((Y-rowm-colm+grand)^2)/((n-1)*(k-1))
  (msr-mse)/(msr+(k-1)*mse+k*(msc-mse)/n)
}
run <- function() {
  group<-split_half_group();save_results(group,'session_split_half_group.csv')
  subject<-split_half_subject();save_results(subject,'session_split_half_subject.csv')
  if(!is.null(subject) && nrow(subject)) {
    s<-aggregate(subject[c('pearson_r','spearman_brown','variance_stability')],by=list(condition=subject$condition),FUN=mean_na)
    save_results(s,'session_split_half_subject_summary.csv')
  }
  invisible(list(group=group,subject=subject))
}
if(sys.nframe()==0)run()
