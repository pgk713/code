# The uploaded MNI to EEG Conversion CODE.ipynb contains no MATLAB-data analyses.
# Its only computation is an MNE fsaverage MRI-to-head transform followed by
# nearest-electrode Euclidean distance. No valid MNE-independent transform can
# be inferred from the notebook alone; do not substitute raw MNI millimeters
# for head coordinates in meters. This R function accepts transformed inputs
# exported from MNE (channel positions and MNI->head transformed position).
source('_r_helpers.R')
mni_xyz_mm <- c(45,51,-12)  # Notebook coordinate; NOT a head-space coordinate.
nearest_eeg_electrodes <- function(head_xyz_m, montage_csv, n=3L) {
  d<-utils::read.csv(montage_csv,check.names=FALSE)
  need_columns(d,c('channel','x_m','y_m','z_m'))
  if(length(head_xyz_m)!=3L || any(!is.finite(head_xyz_m)))stop('Provide the true MNE-transformed head coordinates in meters.')
  xyz<-as.matrix(d[,c('x_m','y_m','z_m')]);ds<-sqrt(rowSums((xyz-matrix(head_xyz_m,nrow(xyz),3,byrow=TRUE))^2))
  data.frame(channel=d$channel,distance_m=ds)[order(ds)[seq_len(min(n,nrow(d)))],,drop=FALSE]
}
# No automatic run: no head transform or standard montage is included.
