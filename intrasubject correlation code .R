# MATLAB CSV analyses corresponding to intrasubject correlation code .ipynb.
# All MNE/SNIRF preprocessing, topomaps and native MNE GLMs are excluded.
source('_r_helpers.R')
# Input: ROI / channel GLM exports with session, condition, subject and beta or tstat.
compare_session_glm <- function(file1='glm_resultsallison1.csv',file2='glm_resultsallison2.csv',statistic='beta') {
  a<-read_input(file1);b<-read_input(file2)
  if('type'%in%names(a))a<-a[tolower(a$type)=='hbo',,drop=FALSE]
  if('type'%in%names(b))b<-b[tolower(b$type)=='hbo',,drop=FALSE]
  keys<-intersect(c('subject','cond','source','detector','ch_name','ChName','ROI'),intersect(names(a),names(b)))
  if(!length(keys))stop('Cannot align sessions: no common GLM identifier columns.')
  m<-merge(a,b,by=keys,suffixes=c('_s1','_s2'))
  v1<-paste0(statistic,'_s1');v2<-paste0(statistic,'_s2')
  need_columns(m,c(v1,v2))
  overall<-data.frame(measure=statistic,t(cor_test(m[[v1]],m[[v2]])))
  save_results(overall,'intrasubject_session_correlation.csv')
  save_results(m,'intrasubject_session_matched.csv')
  if('cond'%in%names(m)){
    z<-do.call(rbind,lapply(split(m,m$cond),function(x)data.frame(condition=x$cond[1],t(cor_test(x[[v1]],x[[v2]])))))
    save_results(z,'intrasubject_by_condition.csv')
  }
  invisible(overall)
}
# Original notebook contains numerous distinct ad-hoc GLM csv paths, so leave
# file choices explicit rather than silently substitute a different analysis.
