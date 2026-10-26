# R translation of the IAPS resolution study's primary 2-session 2x2 RM design.
source('_r_helpers.R')
# Optional library: readxl for IAPS ratings Excel file.
ENJOY_COL <- 'slider.response'; RT_COL <- 'slider.rt'
VALENCE_NEUTRAL_POINT <- 5.0; BONFERRONI_ALPHA <- .025
classify_resolution <- function(path){p<-tolower(gsub('\\\\','/',as.character(path)));ifelse(grepl('degraded stimuli|degraded',p),'low',ifelse(grepl('/stimuli/|stimuli',p),'high','unknown'))}
classify_valence <- function(x) ifelse(!is.finite(x),'unknown',ifelse(x>5,'pos',ifelse(x<5,'neg','neutral')))
read_valence <- function(path=input_file('VALARratings.xlsx'),sheet='Sheet1') {
  if(!requireNamespace('readxl',quietly=TRUE))stop('install.packages("readxl")')
  d<-as.data.frame(readxl::read_excel(path,sheet=sheet));need_columns(d,c('IAPS','valmn'))
  d$IAPS<-round(num(d$IAPS));d$valmn<-num(d$valmn)
  x<-aggregate(d$valmn,by=list(IAPS=d$IAPS),FUN=mean_na);names(x)[2]<-'iaps_valence';x
}
parse_subject_files <- function(folder=DATA_DIR, valence=read_valence()) {
  files<-sort(list.files(folder,pattern='\\.csv$',full.names=TRUE)); if(!length(files))stop('No participant CSVs')
  result<-lapply(files,function(f){
    pid<-regmatches(tolower(f),regexpr('pid[0-9]{3}',tolower(f)))
    sess<-if(grepl('session1|sess1',tolower(f)))1L else if(grepl('session2|sess2',tolower(f)))2L else NA_integer_
    if(!length(pid)||is.na(sess))return(NULL)
    d<-read.csv(f,check.names=FALSE);stim<-if(sess==1)'cond1' else 'cond2'
    need_columns(d,c('IAPS',stim,ENJOY_COL,RT_COL))
    z<-data.frame(pid=pid,session=sess,IAPS=round(num(d$IAPS)),stim_path=as.character(d[[stim]]),enjoy_rating=num(d[[ENJOY_COL]]),rt=num(d[[RT_COL]]))
    z<-z[!is.na(z$IAPS)&is.finite(z$enjoy_rating)&nzchar(z$stim_path),]
    z<-z[!duplicated(z$IAPS,fromLast=TRUE),] # notebook KEEP_DUPLICATE_POLICY = last
    z$resolution_path<-classify_resolution(z$stim_path);z
  })
  result<-Filter(Negate(is.null),result);if(!length(result))stop('No valid subject / session files')
  out<-do.call(rbind,result)
  out<-merge(out,valence,by='IAPS',all.x=TRUE,sort=FALSE)
  out$valence<-classify_valence(out$iaps_valence);out
}
subject_cell_means <- function(d, outcome='enjoy_rating') {
  d<-d[d$valence%in%c('pos','neg')&d$resolution_path%in%c('high','low'),]
  a<-aggregate(d[[outcome]],d[c('pid','valence','resolution_path')],mean_na)
  names(a)[4]<-'value';names(a)[3]<-'resolution';a
}
rm_anova <- function(cell_means) {
  # statsmodels.AnovaRM uses within-subject effects and classical error terms;
  # R aov(Error(pid/(valence*resolution))) computes corresponding 2x2 tests
  cell_means<-cell_means[complete.cases(cell_means),]
  cell_means$pid<-factor(cell_means$pid);cell_means$valence<-factor(cell_means$valence);cell_means$resolution<-factor(cell_means$resolution)
  complete<-names(which(table(cell_means$pid)==4));cell_means<-cell_means[as.character(cell_means$pid)%in%complete,]
  if(!nrow(cell_means))stop('No complete subjects for 2x2 RM ANOVA')
  summary(aov(value~valence*resolution+Error(pid/(valence*resolution)),data=cell_means))
}
resolution_followups <- function(cells) {
  res<-lapply(c('pos','neg'),function(v){d<-cells[cells$valence==v,]
    w<-reshape(d[,c('pid','resolution','value')],idvar='pid',timevar='resolution',direction='wide')
    cbind(valence=v,paired_test(w[['value.high']],w[['value.low']]))
  });do.call(rbind,res)
}
# Exact standard paired-difference SESOI TOST from the notebook.
tost_dz <- function(d,sesoi_dz=.3,alpha=.05) {
  d<-num(d);d<-d[is.finite(d)];n<-length(d);s<-sd(d);se<-s/sqrt(n);bound<-abs(sesoi_dz)*s
  t_lower<-(mean(d)+bound)/se;t_upper<-(mean(d)-bound)/se
  p_lower<-pt(t_lower,n-1,lower.tail=FALSE);p_upper<-pt(t_upper,n-1)
  crit<-qt(1-alpha,n-1)
  data.frame(n=n,SESOI_dz=sesoi_dz,lower_bound=-bound,upper_bound=bound,t_lower=t_lower,p_lower=p_lower,t_upper=t_upper,p_upper=p_upper,TOST_p=max(p_lower,p_upper),CI_low=mean(d)-crit*se,CI_high=mean(d)+crit*se,equivalent=(max(p_lower,p_upper)<alpha))
}
# JZS Bayes factor numerical integration using the source notebook's integral.
jzs_bf10 <- function(t_value,n_pairs,r=sqrt(2)/2) {
  v<-n_pairs-1;N<-as.numeric(n_pairs)
  integrand<-function(g)(1+N*g*r^2)^(-.5)*(1+t_value^2/(v*(1+N*g*r^2)))^(-(v+1)/2)* (2*pi)^(-.5)*g^(-1.5)*exp(-1/(2*g))
  numerator<-integrate(integrand,0,Inf,subdivisions=500L)$value
  numerator/(1+t_value^2/v)^(-(v+1)/2)
}
run <- function() {
  d<-parse_subject_files();save_results(d,'resolution_trials.csv')
  for(outcome in c('enjoy_rating','rt')){
    cells<-subject_cell_means(d,outcome);save_results(cells,paste0('resolution_',outcome,'_cells.csv'))
    print(rm_anova(cells))
    tests<-resolution_followups(cells);tests$bonferroni_significant<-num(tests$p)<BONFERRONI_ALPHA
    save_results(tests,paste0('resolution_',outcome,'_followup.csv'))
  }
  img<-aggregate(d$enjoy_rating,list(IAPS=d$IAPS,iaps_valence=d$iaps_valence),mean_na)
  save_results(data.frame(t(cor_test(img$iaps_valence,img$x))),'resolution_valence_enjoyment_correlation.csv')
  invisible(d)
}
if(sys.nframe()==0)run()
