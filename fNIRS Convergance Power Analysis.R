# MATLAB/CSV-derived block progression and convergence/power analyses.
source('_r_helpers.R')
paired_block_change <- function(d,subject='Subject',block='block',accuracy='accuracy',rt='rt',start=1L,end=9L) {
  need_columns(d,c(subject,block,accuracy,rt))
  x<-d[d[[block]]%in%c(start,end),,drop=FALSE]
  by<-aggregate(x[c(accuracy,rt)],x[c(subject,block)],mean_na)
  w<-reshape(by,idvar=subject,timevar=block,direction='wide')
  list(accuracy=paired_test(w[[paste0(accuracy,'.',end)]],w[[paste0(accuracy,'.',start)]]),rt=paired_test(w[[paste0(rt,'.',end)]],w[[paste0(rt,'.',start)]]))
}
# Piecewise linear model: early slope and post-knot slope, with participant random effects.
fit_mixed_piecewise <- function(d,knot=5,subject='Subject',block='block',dv='accuracy') {
  need_columns(d,c(subject,block,dv));if(!requireNamespace('lme4',quietly=TRUE))stop('install.packages("lme4")')
  dd<-data.frame(subject=factor(d[[subject]]),block=num(d[[block]]),y=num(d[[dv]]))
  dd<-dd[complete.cases(dd),];dd$hinge<-pmax(0,dd$block-knot)
  lme4::lmer(y~block+hinge+(1+block|subject),data=dd,REML=TRUE,control=lme4::lmerControl(optimizer='bobyqa'))
}
# Exhaustive SSE changepoint selection like the source notebook.
simple_changepoint <- function(block,value,min_segment=2) {
  x<-num(block);y<-num(value);valid<-is.finite(x)&is.finite(y);x<-x[valid];y<-y[valid]
  candidates<-sort(unique(x));candidates<-candidates[candidates>=min(x)+min_segment & candidates<=max(x)-min_segment]
  if(!length(candidates))return(NULL)
  fits<-lapply(candidates,function(k){z<-lm(y~x+I(pmax(0,x-k)));data.frame(knot=k,sse=sum(resid(z)^2),slope_before=coef(z)[2],slope_after=coef(z)[2]+coef(z)[3])})
  z<-do.call(rbind,fits);z[which.min(z$sse),,drop=FALSE]
}
run <- function(filename='power_analysis_combined.csv') {
  d<-read_input(filename)
  save_results(data.frame(field=names(d)),'power_input_columns.csv')
  if(all(c('Subject','block','accuracy','rt')%in%names(d))) {
    comparisons<-paired_block_change(d)
    save_results(do.call(rbind,lapply(names(comparisons),function(x)cbind(outcome=x,comparisons[[x]]))),'power_block1_vs_block9.csv')
    if(requireNamespace('lme4',quietly=TRUE)){
      m<-fit_mixed_piecewise(d)
      utils::capture.output(summary(m),file=file.path(OUTPUT_DIR,'power_mixed_piecewise.txt'))
    }
  }
  invisible(d)
}
if(sys.nframe()==0)run()
