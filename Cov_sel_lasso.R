#############################################################################################
##                      Covariate selection based on lasso family                         ###
#############################################################################################
### Vanilla lasso variable importance
### Lasso families for covariate selection include : Lasso,Elastic Net, Adaptive lasso
### Typical lasso, ridge or elastic net
### Note: MLGL (2018) package could be used as alternative

cov.sel.lasso<-function (X,Y,minscreen=2,alpha=alpha){
  ## alpha=1 lasso, alpha=0 ridge, alpha=0.5 elastic net
  family<-ifelse(length(unique(Y))>2,"gaussian","binomial")
  if (!is.matrix(X)) { X<- model.matrix(~-1 + ., X)}
      fitCV <- glmnet::cv.glmnet(x = X, y = Y,alpha=alpha,type.measure = "deviance", family = family)
      whichVariable <- (as.numeric(coef(fitCV$glmnet.fit,s = fitCV$lambda.1se))[-1] !=0)
  if (sum(whichVariable) < minscreen) {
      warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
      sumCoef <- apply(as.matrix(fitCV$glmnet.fit$beta),2,function(x) sum((x != 0)))
      newCut <- which.max(sumCoef >= minscreen)
      whichVariable <- (as.matrix(fitCV$glmnet.fit$beta)[,newCut] != 0)}
  return(whichVariable)}

##########################################################################################
###                            Outcome Apdative Lasso                                  ###
###########################################################################################

### Adaptive lasso variable importance for covariate selection
## Other way for ada-lasso
#tau=1
#first.step.coef=coef(lasso)[-1]
#penalty.factor=abs(first.step.coef+1/sqrt(nrow(x)))^(-tau)
#adalasso=ic.glmnet(x.in,y.in,crit="bic",penalty.factor=penalty.factor)
#pred.adalasso=predict(adalasso,newdata=x.out)


cov.sel.ada_lasso<-function(X,Y,minscreen=2,type.measure="auc"){
  ###: type.measure : auc for catagorical and deviance for continous
  Y<-as.numeric(Y)
  family<-ifelse(length(unique(Y))>2,"gaussian","binomial")
  if (!is.matrix(X)) {X<- model.matrix(~-1 + ., X)}
                      fit_ridge1_cv <- glmnet::cv.glmnet(x = X,y =Y,type.measure =type.measure,family=family, alpha = 0) ## default nfold=10
                      best_ridge_coef <- as.numeric(coef(fit_ridge1_cv,s =fit_ridge1_cv$lambda.min))[-1]
                      fit_alasso1_cv <- glmnet::cv.glmnet(x = X, y = Y,type.measure =type.measure,family=family,alpha = 1,penalty.factor = 1 / abs(best_ridge_coef),
                                      keep = TRUE)
  whichVariable <- (as.numeric(coef(fit_alasso1_cv$glmnet.fit,family=family,s =fit_alasso1_cv$lambda.1se))[-1] !=0)
  if (sum(whichVariable) < minscreen) {
                     warning("fewer than minscreen variables passed the glmnet screen, increased lambda to allow minscreen variables")
                     sumCoef <- apply(as.matrix(fit_alasso1_cv$glmnet.fit$beta),2,function(x) sum((x != 0)))
                     newCut <- which.max(sumCoef >= minscreen)
                     whichVariable <-(as.matrix(fit_alasso1_cv$glmnet.fit$beta)[,newCut] != 0)
}
  return(whichVariable)}

###########################################################################################
###                       Hierarchical Group-Lasso Regularization                       ###
###########################################################################################
## Reference : Lim and Hastie, 2014

detect_factor<-function(vec){
  ## This function will help to detect the catagorical variables in covariate sets
  ## It helps to change catagorical variables as factor
  index<-ifelse(length(unique(vec))>10,FALSE,TRUE)
  return(index)}

cov.sel.G_lasso<-function(X,Y,type="numeric",screenLimit=2)
{ ## Note: change catagorical variables as factor use detect factor function
  ## type: numeric means variable for outcome
  ## catagorical means variable for treatment
  catagorical_index<-apply(X,2,detect_factor)
  X[,catagorical_index]<-lapply(X[,catagorical_index],as.factor)
  family<-ifelse(length(unique(Y))>2,"gaussian","binomial")
  if(type=="catagorical"){
  Y<-ifelse(as.numeric(Y)==1,0,1)}
  i_num <- sapply(X, is.numeric) ## detect numeric variable
  X[, !i_num] <- as.data.frame(apply(X[, !i_num], 2, factor))
  numLevels<- sapply(X,nlevels)
  ## recode the catagorical variable and let them start from 0 and change to them back to interger
  X[, !i_num] <- apply(X[, !i_num], 2, function(col) as.integer(as.factor(col)) - 1)
  numLevels[numLevels==0] <- 1
  ## fit G-lasso
  cv_fit <- glinternet::glinternet.cv(X,Y, numLevels,family=family,screenLimit=screenLimit)
  i_1Std <- which(cv_fit$lambdaHat1Std ==cv_fit$lambda) ## Identify the best lambda by using cross validation
  coefs <- coef(cv_fit$glinternetFit)[[i_1Std]]         ## 10 fold default
  ## identify selected covariates
  sel_main<-names(X)[unique(c(coefs$mainEffects$cat,coefs$mainEffects$cont))]
  sel_main_index<-names(X)%in%sel_main
  sel_inter<-names(X)[unique(c(coefs$interactions$catcat,coefs$interactions$contcont))]
  sel_inter_index<-names(X)%in%sel_inter
  whichVariable<- union(sel_main_index,sel_inter_index)
  return(whichVariable)}

cov.sel.G_lasso(test[,c(1:34)],test$Y)

###########################################################################################
#####                       Lasso covariate selection                                  ####
###########################################################################################

lasso_sel_XY <- function(df,alpha=0.05,minscreen=2,method='adalasso')
      { ## lasso for covariate selection using XY algorithm
        cov_index<-grep("X",names(df))
        X<-df[, cov_index]
        Y<-df[,"Y"]

        T <- factor(df$Treat)
        if(is.null(X)) return(NULL)
        if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
        dat <- data.frame(cbind(X,Y,T))
        Xsubs0 <- X[dat$T == 0,]; Xsubs1 <- X[dat$T == 1,]
        Y0 <- Y[dat$T == 0]; Y1 <- Y[dat$T == 1]
        if(method=='adalasso')
        {index_control <- cov.sel.ada_lasso(X=cbind(Xsubs0),Y=Y0,minscreen=minscreen)
         index_treat <- cov.sel.ada_lasso(X=cbind(Xsubs1),Y=Y1,minscreen=minscreen)}
        else
        {index_control <-cov.sel.G_lasso(X=cbind(Xsubs0),Y=Y0)
        index_treat <- cov.sel.G_lasso(X=cbind(Xsubs1),Y=Y1)}
        X0<-names(X)[index_control]
        X1<-names(X)[index_treat]
        XY <- sort(union(X0, X1))
        if(is.null(XY)) return(NULL)
        if(length(XY) == 0) return(NULL)
        return(sort(XY))}


lasso_sel_XT <-function(df,alpha=0.05,minscreen=2,method='adalasso')
    {## lasso family for variable selection using XT algorithm
    cov_index<-grep("X",names(df))
    X<-df[, cov_index]
    T <-df$Treat
    if(is.null(X)) return(NULL)
    if(sum(is.na(X),is.na(T)) > 0) stop("No missing values are allowed.")
    if (method=='adalasso')
    {index_expose<-cov.sel.ada_lasso(X,Y=T,minscreen=minscreen)}
    else {T<-ifelse(as.numeric(T)==1,0,1)
     index_expose<-cov.sel.G_lasso(X,Y=T,type="catagorical")}
    XT<-names(X)[index_expose]
    if(is.null(XT)) return(NULL)
    if(length(XT) == 0) return(NULL)
    return(XT)}


lasso_sel_double<- function(df,alpha=0.05,minscreen=2,method='adalasso')
{## Double lasso algorithem
  XY<-lasso_sel_XY(df,alpha=0.05,minscreen=2,method=method)
  XT<-lasso_sel_XT(df,alpha=0.05,minscreen=2,method=method)
  X_double<- union(XT,XY)
  if(is.null(X_double)) return(NULL)
  if(length(X_double) == 0) return(NULL)
  return(X_double)}

###########################################################################################
###                          Other: Permutaion based importance                         ###
###########################################################################################

### Lasso with permutation-based variable omportance
LASSO_vimp<-function (x,y,lambda,minscreen = 2, ...){
  nms <- names(x)
  family<-ifelse(length(unique(y))>2,"gaussian","binomial")
  if (!is.matrix(x)) {x<- model.matrix(~-1 + ., x)}
  fitCV <- glmnet::cv.glmnet(x = x, y = y, lambda =lambda,type.measure = "deviance", family = family)
  out<-round(coef(fitCV$glmnet.fit,s= fitCV$lambda.1se)[-1],4)
  names(out) <- nms
  return(out)
}
### Penalty terms that used for covariate selection
lambda<-seq(0,10,0.01)

###Variable importance based on lasso
LASSO_var_sel<-function(x,y,perms = 20, pctl = 0.5,lambda=lambda,minscreen=2,mc.cores = 4){
  if(is.null(x) | length(x) == 0) return(NULL)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(length(p) == 0) return(NULL)
  fac <- which(sapply(x, is.factor))
  not_fac <- setdiff(1:p, fac)
  if(length(fac) > 0) {
    levs <- sapply(x[, fac, drop = FALSE], function(f) length(levels(f)))
    if(any(levs != 2)) stop("Only binary factors allowed.") }
  perm_dat_list <- vector(mode = "list", length = (perms + 1))
  perm_dat_list[[1]] <- x
  perm_dat <- function(dat) {pdat <- dat[sample(x = 1:dim(dat)[1], replace = FALSE), ] }
  perm_dat_list[2:length(perm_dat_list)] <- replicate(n = perms, expr = perm_dat(x), simplify = FALSE)
  imps <- parallel::mclapply(X = perm_dat_list, FUN = LASSO_vimp, y = y,
                             lambda=lambda,minscreen=minscreen,mc.cores = mc.cores)
  imp <- imps[[1]]
  imps<-matrix(unlist(imps[2:length(perm_dat_list)]), nrow = length(imps[[2]]), ncol = perms, byrow = FALSE)
  qtls <- apply(X = imps, MARGIN = 1, FUN = quantile, probs = pctl)
  if(length(which(imp > qtls)) > 0)
  {out <- names(imp)[which(imp > qtls)] }
  else{return(NULL) }
  out <- sort(out)
  out
}
