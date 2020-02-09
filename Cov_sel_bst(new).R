#############################################################################################
##                      Covariate selection based on xgboost                            ###
#############################################################################################
### Pacakges required
install.packages("xgboost")
library(xgboost)
install.packages("ParamHelpers")
library(ParamHelpers)
install.packages("mlr")
library(mlr)

#### Helpful functions
fill_value<-function(vec_name,gain)
{## This function is used to fill the value for selected paramters
     vec_value<-ifelse(vec_name%in%names(gain),gain[vec_name],0)
     return(vec_value)}

BST_vimp<-function(x,y,objective_1="regr.xgboost",objective_2="reg:linear")
    { ### objective_1= "regr.xgboost"or "classif.xgboost"
      ### objective_2= "reg:linear"or "regr:logistic"
      if(is.null(x)) return(NULL)
      if(dim(x)[2] == 0) return(NULL)
      if(sum(is.na(x), is.na(y)) > 0) stop("No missing values are allowed.")
      x<-as.matrix(x)
      if(objective_2 =="binary:logistic")
      { y<-as.factor(y)
        data<-data.frame(x,y)
        ml_task <-mlr::makeClassifTask(data=data, target = "y")}
      if(objective_2 =="reg:linear")
      { data<-data.frame(x,y)
        ml_task <-mlr::makeRegrTask(data=data, target = "y")}
      ###: hyper pramater turnning for  xg-boost
      cv_folds <-mlr:: makeResampleDesc("CV", iters = 5)
      random_tune <-mlr:: makeTuneControlRandom(maxit = 1L)
      model <- mlr:: makeLearner(objective_1)
      xgBoost <- mlr::train(learner = model,task = ml_task)
      ### Set up the hyper-paramter ranges
      model_Params <- ParamHelpers::makeParamSet(
        # number of rounds
        ParamHelpers::makeIntegerParam("nrounds", lower = 50, upper = 100),
        # number of splits in each tree
        ParamHelpers::makeIntegerParam("max_depth", lower = 1, upper = 5),
        #"learing_rate" - prevents overfitting
        ParamHelpers::makeNumericParam("eta", lower = 0.01, upper = 1),
        # Regularization to prevents overfitting
        makeNumericParam("lambda", lower = 1, upper = 10),
        makeNumericParam("alpha", lower = 1, upper = 10),
        # Feature subsamplling to prevent over-fittig
        makeNumericParam("colsample_bytree",lower=0.3, upper=1))
       ### Turned model
      if(objective_2 =="binary:logistic")
      {tuned_model <- mlr::tuneParams(learner = model,
                                       task = ml_task,
                                       resampling = cv_folds,
                                       measures =mlr::mmce,
                                       par.set = model_Params,
                                       control = random_tune,
                                       show.info = FALSE)
      }
      if(objective_2 =="reg:linear"){
        tuned_model <- mlr::tuneParams(learner = model,
                                       task = ml_task,
                                       resampling = cv_folds,
                                       measures =mlr::mse,
                                       par.set = model_Params,
                                       control = random_tune,
                                       show.info = FALSE)}
      ### Turned hyper-paramters
      #nrounds<-tuned_model[[3]]$nrounds
      max_depth<-tuned_model[[3]]$max_depth
      eta<-tuned_model[[3]]$eta
      alpha<-tuned_model[[3]]$alpha
      lambda<-tuned_model[[3]]$lambda
      #colsample_bytree<-tuned_model[[3]]$colsample_bytree
      ### Model fitting
      if(objective_2 =="binary:logistic")
      {xgb_model <- xgboost::xgboost(data=x,
                                     label=ifelse(as.numeric(y)==1,0,1),
                                     nrounds=nrounds,
                                     objective=objective_2,
                                     alpha=alpha,
                                     learning_rate=eta,
                                     lambda=lambda,
                                     learning_rate=eta,
                                     colsample_bytree=colsample_bytree,
                                     max_depth=max_depth)}
      if(objective_2 =="reg:linear")
       {xgb_model <- xgboost::xgboost(data=x,
                                      label=y,
                                      nrounds=200,
                                      objective=objective_2,
                                      alpha=2,
                                      lambda=lambda,
                                      colsample_bytree=colsample_bytree,
                                      learning_rate=eta,
                                      max_depth=5)}
      ### Importance Measure
      imp_measure <-xgboost::xgb.importance(feature_names=xgb_model$feature_names,model=xgb_model)
      var_gain<-round(as.vector(imp_measure$Gain),4)### Gain as importace measure
      names(var_gain)<-imp_measure$Feature
      store_value<-as.vector(rep(0,dim(x)[2]))
      x<-as.data.frame(x)
      names(store_value)<-names(x)
      for (i in names(store_value))
      {store_value[i]<-fill_value(i,var_gain)}
      store_value<- round(store_value,4)
      return(store_value)}

### Permutation importance for boosted tree

BTVS_perm_trad <- function(x, y, perms = 40,pctl = .95,
                           objective_1="classif.xgboost",
                           objective_2="binary:logistic")
{ ###  Permutation based variable importances for xgboost
  ###  Perms is the number of permutations
  if(is.null(x) | length(x) == 0) return(NULL)
  if(dim(x)[2] == 0) return(NULL)
  if(sum(is.na(x), is.na(y)) > 0) stop("No missing values are allowed.")
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(length(p) == 0) return(NULL)
  perm_dat_list <- vector(mode = "list", length = (perms + 1))
  perm_dat_list[[1]] <- x
  perm_dat <- function(dat) { pdat <- dat[sample(x = 1:dim(dat)[1], replace = FALSE), , drop = FALSE] }
  perm_dat_list[2:length(perm_dat_list)] <- replicate(n = perms, expr = perm_dat(x), simplify = FALSE)
  imps <- lapply(X = perm_dat_list,
                 FUN =BST_vimp, y = y,
                 objective_1=objective_1,
                 objective_2=objective_2)
  imp <- imps[[1]]
  imps <- matrix(unlist(imps[2:length(perm_dat_list)]), nrow = length(imps[[2]]), ncol = perms, byrow = FALSE)
  pctls <- apply(X = imps, MARGIN = 1, FUN = quantile, probs = pctl)
  ord <- order(imp, decreasing = TRUE)
  imp_ord <- imp[ord]; pctls_ord <- pctls[ord]
  out <- cbind(imp_ord, pctls_ord)
  out <- data.frame(out)
  rownames(out) <- names(imp_ord)
  colnames(out) <- c("Imp", "Pctl")
  return(out)
}




###  (SHapley Additive exPlanation)



###########################################################################################
####                         Xgboost  covariate selection                             ####
###########################################################################################

### XY algorthems for xgboost
BT_sel_XY <- function(df,perms = 20,pctl = .95) {
  cov_index<-grep("X",names(df))
  X<-df[, cov_index]
  Y<-df[,"Y"]
  T <- factor(df$T)
  if(is.null(X)) return(NULL)
  if(dim(X)[2] == 0) return(NULL)
  if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
  if(length(table(T)) != 2) stop("Exposure variable T must have exactly two categories.")
  T <- factor(T); levels(T) <- 0:1
  Xsubs0 <- X[T == 0, , drop = FALSE]; Xsubs1 <- X[T == 1, , drop = FALSE]
  Y0 <- Y[T == 0]; Y1 <- Y[T == 1]
  X0 <- BTVS_perm_trad(x = Xsubs0, y = Y0,perms = perms, pctl = pctl,
                       objective_1="regr.xgboost",
                       objective_2="reg:linear")
  X1 <- BTVS_perm_trad(x = Xsubs1, y = Y1,perms = perms, pctl = pctl,
                       objective_1="regr.xgboost",
                       objective_2="reg:linear")
  X0 <- rownames(X0)[which(X0[,1] > 0 & X0[,1] >= X0[,2])]
  X1 <- rownames(X1)[which(X1[,1] > 0 & X1[,1] >= X1[,2])]
  XY <- sort(union(X0, X1))
  if(is.null(XY)) return(NULL)
  if(length(XY) == 0) return(NULL)
  else return (XY)}



### XT algorthems for xgboost
BT_sel_XT <- function(df,perms = 20, pctl = .95)
    {cov_index<-grep("X",names(df))
      X<-df[, cov_index]
      T <- df$T
      if(is.null(X)) return(NULL)
      if(dim(X)[2] == 0) return(NULL)
      if(sum(is.na(X), is.na(T)) > 0) stop("No missing values are allowed.")
      if(length(table(T)) != 2) stop("Exposure variable T must have exactly two categories.")
      X<-as.data.frame(X)
      perms=20
      pctl=0.95
      XT<-BTVS_perm_trad(x=X,y=T,perms = perms,
                         pctl = pctl,
                         objective_1="classif.xgboost",
                         objective_2="binary:logistic")
      XT <- rownames(XT)[which(XT[,1] > 0)]
      if(is.null(XT)) return(NULL)
      if(length(XT) == 0) return(NULL)
      else return (XT)}


