# R translation of MATLAB-exported GLM beta analyses from fNIRS Functional Conectivity Analysis.ipynb
# Granger ordering follows the notebook's row order; interpret cautiously: beta rows are NOT time samples.
source('_r_helpers.R')
GLM_CSV <- 'glm_result_group_v3.csv'
# Cross correlation mirrors: cor(series1[:-lag], series2[lag:]) for lag=1.
crosscorr_lag <- function(x,y,lag=1L) {
  n<-min(length(x),length(y));if(lag>=n)return(NA_real_)
  unname(cor_test(x[seq_len(n-lag)],y[(lag+1):n])['r'])
}
regional_correlations <- function(regional_beta,method='pearson',lag=1L) {
  m<-cor_matrix(regional_beta,method)
  cross<-outer(colnames(regional_beta),colnames(regional_beta),Vectorize(function(a,b) crosscorr_lag(regional_beta[[a]],regional_beta[[b]],lag)))
  dimnames(cross)<-list(colnames(regional_beta),colnames(regional_beta))
  list(correlation=m,lagged=cross)
}
partial_correlation_global <- function(beta) {
  x<-as.data.frame(beta);global<-rowMeans(x,na.rm=TRUE)
  out<-matrix(NA_real_,ncol(x),ncol(x),dimnames=list(names(x),names(x)));diag(out)<-1
  for(i in seq_len(ncol(x)))for(j in seq_len(ncol(x)))if(i!=j){
    z<-data.frame(x=num(x[[i]]),y=num(x[[j]]),g=global);z<-z[complete.cases(z),]
    if(nrow(z)>3)out[i,j]<-unname(cor_test(residuals(lm(x~g,data=z)),residuals(lm(y~g,data=z)))['r'])
  };out
}
# Statsmodels grangercausalitytests(..., maxlag=5)['ssr_chi2test'] differs from
# R's conventional F-test; compute equivalent SSR chi-square approximation.
granger_ssr_chi2 <- function(effect,cause,lag=5L) {
  n<-min(length(effect),length(cause));y<-num(effect[seq_len(n)]);x<-num(cause[seq_len(n)])
  if(n<=3*lag+3)return(NA_real_)
  tab<-data.frame(y=y[(lag+1):n]);for(k in seq_len(lag)){
    tab[[paste0('yl',k)]]<-y[(lag+1-k):(n-k)]
    tab[[paste0('xl',k)]]<-x[(lag+1-k):(n-k)]
  };tab<-tab[complete.cases(tab),,drop=FALSE];if(nrow(tab)<3*lag+2)return(NA_real_)
  reduced<-lm(reformulate(paste0('yl',seq_len(lag)),'y'),data=tab)
  full<-lm(reformulate(c(paste0('yl',seq_len(lag)),paste0('xl',seq_len(lag))),'y'),data=tab)
  rss_reduced<-sum(resid(reduced)^2);rss_full<-sum(resid(full)^2)
  chi<-nrow(tab)*(rss_reduced-rss_full)/rss_full
  stats::pchisq(chi,df=lag,lower.tail=FALSE)
}
run <- function() {
  d<-read_input(GLM_CSV)
  save_results(data.frame(column=names(d),type=vapply(d,function(x)class(x)[1],'')),'FC_input_structure.csv')
  # The original notebook creates region-specific beta vectors from its
  # channels / conditions. An explicit `Region`, `beta`, `condition` export
  # enables precisely the same region correlation step without MNE.
  if(!all(c('Region','beta')%in%names(d))) {
    message('For region analyses export a CSV containing Region, beta and optionally condition in row order.')
    return(invisible(d))
  }
  groups<-split(d,if('condition'%in%names(d))d$condition else rep('all',nrow(d)))
  for(cond in names(groups)){
    sub<-groups[[cond]]; l<-split(num(sub$beta),sub$Region);lens<-vapply(l,length,0L)
    if(length(l)<2)next
    mm<-as.data.frame(lapply(l,function(x)c(x,rep(NA_real_,max(lens)-length(x)))),check.names=FALSE)
    out<-regional_correlations(mm)
    save_results(as.data.frame(out$correlation),paste0('FC_beta_Pearson_',cond,'.csv'))
    save_results(as.data.frame(out$lagged),paste0('FC_beta_lagged_',cond,'.csv'))
    save_results(as.data.frame(cor_matrix(mm,'spearman')),paste0('FC_beta_Spearman_',cond,'.csv'))
    save_results(as.data.frame(partial_correlation_global(mm)),paste0('FC_beta_partial_',cond,'.csv'))
    g<-outer(names(mm),names(mm),Vectorize(function(i,j)if(i==j)NA_real_ else granger_ssr_chi2(mm[[i]],mm[[j]],5)))
    dimnames(g)<-list(names(mm),names(mm));save_results(as.data.frame(g),paste0('FC_granger_lag5_p_',cond,'.csv'))
  }; invisible(NULL)
}
if(sys.nframe()==0)run()
