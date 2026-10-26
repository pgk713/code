# MATLAB-exported functional connectivity analysis port.
# Network clustering is based on exported HbO channel-pair FC and ROI labels.
source('_r_helpers.R')
GROUP_CSV <- 'EI_ISC_group_functional_conn.csv'
SUBJECT_CSV <- 'functional_conn_subjectlevel_pruned2hz_cond.csv'
run <- function() {
  d<-read_input(GROUP_CSV);s<-read_input(SUBJECT_CSV)
  d<-prune_hbo_fc(d)
  if('R'%in%names(d)) {
    bycond<-split(d,d$condition)
    for(cond in names(bycond)) {
      dd<-bycond[[cond]]
      # Pair means and ROI means preserve notebook groupby semantics.
      if(all(c('ROI_X','ROI_Y')%in%names(dd))){
        fc<-roi_connectivity(dd,'R');names(fc)[ncol(fc)]<-'R_mean'
        save_results(fc,paste0('matlab_FC_ROI_',cond,'.csv'))
        m<-xtabs(R_mean~ROI_X+ROI_Y,fc);plot_matrix(m,paste('ROI connectivity',cond),paste0('matlab_FC_heatmap_',cond,'.png'))
      }
    }
    if(length(bycond)>=2){
      lev<-names(bycond)[1:2];a<-bycond[[lev[1]]];b<-bycond[[lev[2]]]
      z<-merge(a[,c('pair_id','R')],b[,c('pair_id','R')],by='pair_id')
      if(nrow(z)>2)save_results(data.frame(conditions=paste(lev,collapse=' vs '),t(cor_test(z$R.x,z$R.y))),'matlab_FC_across_conditions.csv')
    }
  }
  if(all(c('TypeOrigin','TypeDest','SourceOrigin','DetectorOrigin','SourceDest','DetectorDest')%in%names(s))){
    s<-prune_hbo_fc(s)
  }
  if(all(c('subject','condition','z_score')%in%names(s))){
    stats<-do.call(rbind,lapply(split(s,list(s$subject,s$condition),drop=TRUE),function(dd){
      data.frame(subject=dd$subject[1],condition=dd$condition[1],mean_z=mean_na(num(dd$z_score)),n=nrow(dd))
    }))
    save_results(stats,'matlab_FC_subject_fisher_z.csv')
    if(length(unique(stats$condition))==2){
      w<-reshape(stats[,c('subject','condition','mean_z')],idvar='subject',timevar='condition',direction='wide')
      save_results(w,'matlab_FC_subject_paired.csv')
      if(ncol(w)==3) save_results(paired_test(w[[2]],w[[3]]),'matlab_FC_condition_paired_ttest.csv')
    }
  }
  invisible(list(group=d,subject=s))
}
if(sys.nframe()==0)run()
