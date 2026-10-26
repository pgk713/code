# R translation of neural network *downstream statistical analyses* only.
# Source's timm/PyTorch/decord feature generation is Python-specific and
# cannot be recreated by simply loading its weights into this R script.
source('_r_helpers.R')
ALPHA_GRID <- c(.001,.01,.1,1,10,100)
L1_RATIO_GRID <- c(.1,.3,.5,.7,.9)
N_VIDEO_PCS <- 12L;INNER_CV_SPLITS <- 5L;RANDOM_STATE <- 42L
pop_scale <- function(X) {X<-as.matrix(X);mu<-colMeans(X);sigma<-sqrt(colMeans(sweep(X,2,mu,'-')^2));sigma[sigma==0]<-1;sweep(sweep(X,2,mu,'-'),2,sigma,'/')}
# This function expects precomputed pooled/video embeddings in an ordinary
# numeric CSV (one row per video) instead of raw .npy arrays.
across_video_pca <- function(embeddings,index_cols=c('VideoKey','Trigger'),n_pcs=N_VIDEO_PCS) {
  need_columns(embeddings,index_cols)
  X<-as.matrix(embeddings[,setdiff(names(embeddings),index_cols),drop=FALSE]);storage.mode(X)<-'double'
  # Python StandardScaler uses population SD (ddof=0), unlike scale() in R.
  if(anyNA(X))stop('Source PCA does not impute embeddings; input must be complete.')
  pcs<-prcomp(pop_scale(X),center=FALSE,scale.=FALSE,rank.=min(n_pcs,ncol(X),nrow(X)))
  P<-matrix(0,nrow(X),n_pcs);k<-min(ncol(pcs$x),n_pcs);P[,seq_len(k)]<-pcs$x[,seq_len(k),drop=FALSE]
  names_p<-sprintf('VideoPC%02d',seq_len(n_pcs));colnames(P)<-names_p
  var<-pcs$sdev^2;list(features=cbind(embeddings[index_cols],as.data.frame(P)),scree=data.frame(PC=sprintf('PC%02d',seq_along(var)),ExplainedVarianceRatio=var/sum(var),CumulativeVarianceExplained=cumsum(var/sum(var))))
}
# Nested leave-one-video-out elastic net with training-fold-only mean imputation.
# sklearn and glmnet share elastic-net penalty for the same alpha/l1_ratio;
# sklearn's shuffled KFold splits differ from R sample() at the same seed.
nested_loocv_elastic_net <- function(X,y,alpha_grid=ALPHA_GRID,l1_ratio_grid=L1_RATIO_GRID,inner_splits=INNER_CV_SPLITS,seed=RANDOM_STATE) {
  if(!requireNamespace('glmnet',quietly=TRUE))stop('install.packages("glmnet")')
  X<-as.matrix(X);storage.mode(X)<-'double';y<-num(y);n<-length(y);if(nrow(X)!=n)stop('Shape mismatch')
  predict_fold <- function(train,test,a,r) {
    xt<-X[train,,drop=FALSE];xv<-X[test,,drop=FALSE]
    means<-apply(xt,2,mean_na);means[!is.finite(means)]<-0
    for(j in seq_len(ncol(xt))){xt[is.na(xt[,j]),j]<-means[j];xv[is.na(xv[,j]),j]<-means[j]}
    mu<-colMeans(xt);sc<-sqrt(colMeans(sweep(xt,2,mu,'-')^2));sc[sc==0]<-1
    xt<-sweep(sweep(xt,2,mu,'-'),2,sc,'/');xv<-sweep(sweep(xv,2,mu,'-'),2,sc,'/')
    model<-glmnet::glmnet(xt,y[train],alpha=r,lambda=a,intercept=TRUE,standardize=FALSE,thresh=1e-7,maxit=100000)
    as.numeric(predict(model,newx=xv,s=a,exact=FALSE))
  }
  yhat<-rep(NA_real_,n);records<-vector('list',n)
  for(test in seq_len(n)) {
    train<-setdiff(seq_len(n),test);k<-min(inner_splits,length(train));if(k<2)stop('Not enough samples')
    set.seed(seed);order<-sample(seq_along(train));folds<-rep(seq_len(k),length.out=length(train));foldid<-integer(length(train));foldid[order]<-folds
    grid<-expand.grid(alpha=alpha_grid,l1_ratio=l1_ratio_grid);loss<-numeric(nrow(grid))
    for(g in seq_len(nrow(grid))){err<-numeric(k)
      for(f in seq_len(k)) {va<-train[foldid==f];tr<-train[foldid!=f]
        yh<-predict_fold(tr,va,grid$alpha[g],grid$l1_ratio[g]);err[f]<-mean((y[va]-yh)^2)
      };loss[g]<-mean(err)
    }
    best<-which.min(loss);yhat[test]<-predict_fold(train,test,grid$alpha[best],grid$l1_ratio[best])
    records[[test]]<-data.frame(fold=test,test_index_python=test-1,enet_alpha=grid$alpha[best],enet_l1_ratio=grid$l1_ratio[best],best_inner_score_neg_mse=-loss[best])
  }
  list(yhat=yhat,r=unname(cor_test(y,yhat)['r']),r2=1-sum((y-yhat)^2)/sum((y-mean(y))^2),best_params=do.call(rbind,records))
}
run <- function(embeddings_csv='video_embeddings_mean.csv',outcome_csv='video_performance.csv') {
  E<-read_input(embeddings_csv);Y<-read_input(outcome_csv)
  z<-across_video_pca(E);save_results(z$features,'video_pca_features.csv');save_results(z$scree,'video_pca_scree.csv')
  # Outcome matching requires explicit keys rather than aligning unequal row orders.
  need_columns(Y,c('VideoKey','Trigger','performance'))
  d<-merge(z$features,Y[,c('VideoKey','Trigger','performance')],by=c('VideoKey','Trigger'))
  fit<-nested_loocv_elastic_net(d[,grep('^VideoPC',names(d)),drop=FALSE],d$performance)
  save_results(data.frame(VideoKey=d$VideoKey,Trigger=d$Trigger,observed=d$performance,predicted=fit$yhat),'neural_net_downstream_predictions.csv')
  save_results(fit$best_params,'neural_net_nested_cv_params.csv');print(c(r=fit$r,r2=fit$r2));invisible(fit)
}
if(sys.nframe()==0)run()
