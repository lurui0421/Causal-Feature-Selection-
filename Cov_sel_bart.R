######################################################################################################
###              Bayesian Additive Regression Tree (BART) used for feature selection               ###
######################################################################################################
### Author: Rui Lu 
### Reference:Bleich et al, 2014
### Reference: Linero et al,2018
### Covariate seleciton based on BART feature incluing importance 

BART_vimp<-function(X,y,ntree=20,mc.cores =4,type="wbart",sparse=TRUE){
  ## This function estimate the propotion of the variables that are used in bart fitting
  ## Wbart is for comtimous bart, lbart for logit bart (treatment)
  ## sparse=True will enable the use of sparse Dirichlet for variable selection 
  if (type=='lbart')
  {y<-ifelse(as.numeric(y)==1,0,1)
  }
  model<-BART::mc.gbart(X,y,type=type,ndpost=200L,ntree=ntree, mc.cores = mc.cores,sparse=sparse)
  Freq<-as.data.frame(model$varcount)
  names(Freq)<-names(X)
  ##proportion of times that a split is using certain variables appear 
  # in all split 
  Props<-round(apply(Freq,2,sum)/sum(Freq),4)
  return(Props)
}

### Permutaion based variable importance 
BART_var_sel<-function(x,y,perms = 20,pctl = 0.05,type=type,mc.cores = 4){
  ### Permutation based BART method for covariate selection 
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
  imps <- parallel::mclapply(X=perm_dat_list,FUN=BART_vimp,y=y,type=type,mc.cores=mc.cores)
  imp <- imps[[1]]
  imps<-matrix(unlist(imps[2:length(perm_dat_list)]), nrow = length(imps[[2]]), ncol = perms, byrow = FALSE)
  qtls <- apply(X = imps, MARGIN = 1, FUN = quantile, probs = pctl)
  if(length(which(imp > qtls)) > 0) 
  {out <- names(imp)[which(imp > qtls)] } 
  else{return(NULL) }
  out <- sort(out)
  out
}

### Covariate selection function 
bart_sel <- function(df,perms = perms, pctl = pctl, algorithem="1") {
  ### Change the interger value into numeric for further analysis 
  ### algorithem 1 is the XY for DeLuna
  ### algorithem 2 is the VanderVee and Shipster 
  cov_index<-grep("Treat",names(df))-2
  X<-df[,1:cov_index]
  Y<-df[,"Y"]
  T <- factor(df$Treat)
  if(is.null(X)) return(NULL)
  if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
  dat <- data.frame(cbind(X,Y,T))
  XsubsT<-data.frame(X,T)
  perms = 20
  pctl = 0.95
  algorithem="1"
  XT<-BART_var_sel(X,y=XsubsT$T,type='lbart',perms = perms, pctl = pctl)
  Xsubs0 <- X[dat$T == 0,]; Xsubs1 <- X[dat$T == 1,]
  Y0 <- Y[dat$T == 0]; Y1 <- Y[dat$T == 1]
  X0 <- BART_var_sel(Xsubs0,y=Y0,type='wbart',perms = perms, pctl = pctl)
  X1 <- BART_var_sel(Xsubs1,y=Y1,type='wbart', perms = perms, pctl = pctl)
  XY <- sort(union(X0, X1))
  if (algorithem=="1"){return(XY=XY)}
  else
  { XTY<-sort(union(XT, XY))
    if(is.null(XTY)) return(NULL)
    if(length(XTY) == 0) return(NULL)
    return(XTY)}
  }







