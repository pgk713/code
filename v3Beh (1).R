# Behavioral analyses port for v3Beh (1).ipynb; does not implement MNE.
# Accepts one CSV per participant/session or an already-combined trial table.
source('_r_helpers.R')
load_behavioral_folder <- function(folder) {
  files<-list.files(folder,pattern='\\.csv$',full.names=TRUE,recursive=TRUE)
  if(!length(files))stop('No CSV files in behavioral folder ', folder)
  xs<-lapply(files,function(path){d<-tryCatch(read.csv(path,check.names=FALSE),error=function(e)NULL);if(is.null(d))return(NULL)
    d$source_file<-basename(path);d$source_session<-basename(dirname(path));d
  });xs<-Filter(Negate(is.null),xs)
  # Unequal column layouts are combined by filling unavailable fields with NA.
  allnames<-unique(unlist(lapply(xs,names)))
  do.call(rbind,lapply(xs,function(d){for(n in setdiff(allnames,names(d)))d[[n]]<-NA;d[,allnames,drop=FALSE]}))
}
# Python notebook uses participant-level correctness/RT, block progression,
# paired session accuracy and RT. Specify matching column names if needed.
behavioral_summary <- function(df,subject='participant',condition='condition',correct='correct',rt='rt',block='block') {
  need_columns(df,c(subject,correct,rt))
  k<-c(subject,intersect(c(condition,block,'source_session'),names(df)))
  group<-interaction(df[k],drop=TRUE,lex.order=TRUE)
  chunks<-split(df,group)
  z<-lapply(chunks,function(d){
    rec<-as.list(d[1,k,drop=FALSE]);rec$accuracy<-mean_na(num(d[[correct]]));rec$mean_rt<-mean_na(num(d[[rt]]));rec$n_trials<-nrow(d)
    as.data.frame(rec)
  });out<-do.call(rbind,z);row.names(out)<-NULL;out
}
compare_behavioral_sessions <- function(summary,subject='participant',session='source_session') {
  need_columns(summary,c(subject,session,'accuracy','mean_rt'))
  levels<-unique(summary[[session]]);if(length(levels)!=2L)stop('Exactly two session names required.')
  w<-reshape(summary[,c(subject,session,'accuracy','mean_rt')],idvar=subject,timevar=session,direction='wide')
  save_results(w,'v3_behavioral_sessions_wide.csv')
  out<-rbind(cbind(outcome='accuracy',paired_test(w[[paste0('accuracy.',levels[1])]],w[[paste0('accuracy.',levels[2])]])),cbind(outcome='mean_rt',paired_test(w[[paste0('mean_rt.',levels[1])]],w[[paste0('mean_rt.',levels[2])]])))
  save_results(out,'v3_behavioral_paired_tests.csv');out
}
run <- function(){d<-load_behavioral_folder(input_file('behavioral_data'));s<-behavioral_summary(d);save_results(s,'v3_behavioral_summary.csv');invisible(s)}
# Configure column names before calling run() when the notebook's raw files
# use alternate response/accuracy labels. run() is deliberately opt-in.
