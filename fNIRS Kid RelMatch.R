# Source notebook fNIRS Kid RelMatch.ipynb is entirely an MNE/SNIRF
# preprocessing + GLM workflow: read_raw_snirf, HbO preprocessing, event
# construction, ROI GLM, topographic visualization. Per user request, those
# MNE-specific computations are intentionally NOT replicated in R.
# MATLAB-exported GLM CSVs, if available, can be summarized by this function.
source('_r_helpers.R')
review_matlab_glm <- function(file_name='glm_result_group_kid_relmatch.csv') {
  d<-read_input(file_name)
  if(all(c('type','beta','cond')%in%names(d))){
    d<-d[tolower(d$type)=='hbo',,drop=FALSE]
    by<-intersect(c('cond','ROI','region','subject'),names(d))
    if(length(by)){a<-aggregate(num(d$beta),d[by],mean_na);names(a)[ncol(a)]<-'mean_beta';save_results(a,'kid_relmatch_glm_summary.csv');return(a)}
  }
  d
}
# Call review_matlab_glm() only if an appropriate export actually exists.
