# R translation: MATLAB group functional connectivity analyses (HbO, two conditions)
source('_r_helpers.R')
INPUT_CSV <- 'EI_ISC_group_functional_conn_ROI_tvr.csv'
BAD_DETECTORS <- 16:23
run <- function() {
  d<-read_input(INPUT_CSV)
  need_columns(d,c('condition','R'))
  d<-prune_hbo_fc(d,BAD_DETECTORS)
  conds<-c('positive','negative')
  a<-d[d$condition==conds[1],,drop=FALSE];b<-d[d$condition==conds[2],,drop=FALSE]
  # Source notebook removes mirror duplicates separately and compares rows by position,
  # rather than joining on pair_id; report that exact unaligned comparison as well.
  if(nrow(a)!=nrow(b)) stop('Source notebook would stop here: condition vector lengths do not match.')
  direct<-cor_test(a$R,b$R)
  save_results(data.frame(comparison='positive_vs_negative_original_row_order',t(direct)), 'EI_FC_original_order_correlation.csv')
  # Additional diagnostic: explicitly pair-matched check (not part of original test).
  matched<-merge(a[,c('pair_id','R')],b[,c('pair_id','R')],by='pair_id',suffixes=c('_positive','_negative'))
  if(nrow(matched)>=3)save_results(data.frame(comparison='pair_matched_diagnostic',t(cor_test(matched$R_positive,matched$R_negative))),'EI_FC_pair_matched_diagnostic.csv')
  if(all(c('ROI_X','ROI_Y')%in%names(d))) {
    roi<-roi_connectivity(d,'R');names(roi)[ncol(roi)]<-'R_mean'
    save_results(roi,'EI_FC_ROI_summary.csv')
    for(cc in unique(d$condition)){
      sub<-roi[roi$condition==cc,,drop=FALSE]
      m<-xtabs(R_mean~ROI_X+ROI_Y, data=sub)
      plot_matrix(m,paste('ROI FC:',cc),paste0('EI_FC_ROI_',gsub('[^a-zA-Z0-9]','_',cc),'.png'))
    }
  }
  invisible(direct)
}
if(sys.nframe()==0)run()
