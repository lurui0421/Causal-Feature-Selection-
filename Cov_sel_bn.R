############################################################################################
##                      Covariate selection based on Bayesian Networks                    ##
############################################################################################
### Author: Rui Lu
### Reference: Jenny Haggstorm 2018
### Chpater 2: Pratical approaches to causal relationship exploration (2015) Perason
###: There are two algorthems provided : learn pc and learn mb
### Discriticize function


detect_factor<-function(vec){
### This function will help to detect the catagorical variables in covariate sets
  index<-ifelse(length(unique(vec))>20,FALSE,TRUE)
  return(index)}


discret<-function(df){
###: This function dicreticize the continous variable into intervals.
  catagorical_index<-apply(df,2,detect_factor)
  df[,catagorical_index]<-lapply(df[,catagorical_index],as.factor)
  cov_index<-grep("Treat",names(df))-2
  X<-df[,1:cov_index]
  Y<-df[,"Y"]
  T <- factor(df$Treat)
  ### reconstruct the data set
  dat<-data.frame(X,T,Y)
  covarcol<-1:(dim(X)[2])
  ycol<-(dim(X)[2])+2
  covars<-colnames(dat[,covarcol])
  Tcol<-ycol-1
  if(class(dat[,Tcol])!="factor"){dat[,Tcol]<-factor(dat[,Tcol])}
  numcol<-which(lapply(dat[covarcol],class)=="numeric")
  if(length(unique(dat[,ycol]))>2){
    tried<-try(bnlearn::discretize(data.frame(dat[,ycol]),method="quantile"),silent=TRUE)
    if(class(tried)=="data.frame"){
       dat[,ycol]<-bnlearn::discretize(data.frame(dat[,ycol]),method="quantile")}
      else{
      print(tried)
    stop("the numeric outcome could not be discretized, recode it into factor")}}
  dat[,numcol]<-bnlearn::discretize(data.frame(dat[,numcol]),method="quantile")
  if(length(numcol)>0){
    if(length(numcol)==1){
      tried<-try(bnlearn::discretize(data.frame(dat[,numcol]),method="quantile"),silent=TRUE)
      if(class(tried)=="data.frame"){
        dat[,numcol]<-bnlearn::discretize(data.frame(dat[,numcol]),method="quantile")
      }else{ print(tried)
       stop("the numeric covariate could not be discretized, recode it into factor")}}
       else{
       tried<-try(bnlearn::discretize(dat[,numcol],method="quantile"),silent=TRUE)
      if(class(tried)=="data.frame"){
        dat[,numcol]<-bnlearn::discretize(dat[,numcol],method="quantile")
      }
      else{ print(tried)
      stop("at least one numeric covariate could not be discretized, recode it into factor")}
    }
  }
     return (dat)
}



## classic covariate selection functions
## MMPC
bn_sel_mmpc_XT <- function(df,alpha=0.05,discretize=FALSE) {
  ### Change the interger value into numeric for further analysis
  ### alpha value is for mutual information based tests
  if(discretize==TRUE)
  {df<-discret(df)}
   cg_index<-names(which(unlist(lapply(df,class))=="integer"))
   for (i in cg_index){
       df[,i]<-as.numeric(df[,i])
  }
  cov_index<-grep("X",names(df))
  X<-df[, cov_index]
  T <- as.numeric(df$T)
  if(is.null(X)) return(NULL)
  if(sum(is.na(X), is.na(T)) > 0) stop("No missing values are allowed.")
  XsubsT<-data.frame(X,T)
  XT<-bnlearn::learn.nbr(XsubsT, "T", "mmpc",alpha=alpha)
  if(is.null(XT)) return(NULL)
  if(length(XT) == 0) return(NULL)
  return(XT)}


bn_sel_mmpc_XY <- function(df,alpha=0.05,discretize=TRUE) {
  ### Change the interger value into numeric for further analysis
  ### alpha value is for mutual information based tests
  if(discretize==TRUE)
  {df<-discret(df)}
  cg_index<-names(which(unlist(lapply(df,class))=="integer"))
  for (i in cg_index){
    df[,i]<-as.numeric(df[,i])
  }
  cov_index<-grep("X",names(df))
  X<-df[, cov_index]
  Y<-df[,"Y"]
  T <- df$T
  if(is.null(X)) return(NULL)
  if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
  dat <- data.frame(cbind(X,Y,T))
  Xsubs0 <- X[dat$T == 0,]; Xsubs1 <- X[dat$T == 1,]
  Y0 <- Y[dat$T == 0]; Y1 <- Y[dat$T == 1]
  X0 <- bnlearn::learn.nbr(x = cbind(Xsubs0, Y0), node = "Y0", method = "mmpc", alpha = alpha)
  X1 <- bnlearn::learn.nbr(x = cbind(Xsubs1, Y1), node = "Y1", method = "mmpc", alpha = alpha)
  XY <- sort(union(X0, X1))
  if(is.null(XY)) return(NULL)
  if(length(XY) == 0) return(NULL)
  else return (XY)}

### IAMB
### Reference: Tsamardinos, I., Aliferis, C. F., Statnikov, A. R., & Statnikov, E. (2003)
bn_sel_mb_XY <- function(df,alpha=0.05,discretize=TRUE) {
  ### Change the interger value into numeric for further analysis
  ### alpha value is for mutual information based tests
  ### algorithem 1: XY for DeLuna
  if(discretize==TRUE)
  {df<-discret(df)}
  cg_index<-names(which(unlist(lapply(df,class))=="integer"))
  for (i in cg_index){
    df[,i]<-as.numeric(df[,i])
  }
  cov_index<-grep("X",names(df))
  X<-df[, cov_index]
  Y<-df[,"Y"]
  T <- df$T
  if(is.null(X)) return(NULL)
  if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
  dat <- data.frame(cbind(X,Y,T))
  Xsubs0 <- X[dat$T == 0,]; Xsubs1 <- X[dat$T == 1,]
  Y0 <- Y[dat$T == 0]; Y1 <- Y[dat$T == 1]
  X0 <- bnlearn::learn.mb(x = cbind(Xsubs0, Y0), node = "Y0", method = "iamb", alpha = alpha)
  X1 <- bnlearn::learn.mb(x = cbind(Xsubs1, Y1), node = "Y1", method = "iamb", alpha = alpha)
  XY <- sort(union(X0, X1))
  if(is.null(XY)) return(NULL)
  if(length(XY) == 0) return(NULL)
  else return (XY)}


bn_sel_mb_XT <- function(df,alpha=0.25,discretize=FALSE) {
  ### Change the interger value into numeric for further analysis
  ### alpha value is for mutual information based tests
  ### algorithem 2: VanderVee and Shipster "
  if(discretize==TRUE)
  {df<-discret(df)}
  cg_index<-names(which(unlist(lapply(df,class))=="integer"))
  for (i in cg_index){
    df[,i]<-as.numeric(df[,i])
  }
  cov_index<-grep("X",names(df))
  X<-df[, cov_index]
  T <- as.numeric(df$T)
  if(is.null(X)) return(NULL)
  if(sum(is.na(X),is.na(T)) > 0) stop("No missing values are allowed.")
  XsubsT<-data.frame(X,T)
  XT<-bnlearn::learn.mb(XsubsT, node ="T", method ="iamb",alpha=alpha)
  if(is.null(XT)) return(NULL)
  if(length(XT) == 0) return(NULL)
  return(XT)}


#################################################################################
##                         Other covariate selecrion based on BN              ##
#################################################################################

bn_sel_hi_pc <- function(df,alpha=0.05,discretize=FALSE, algorithem="1") {
  ### Change the interger value into numeric for further analysis
  ### alpha value is for mutual information based tests
  ### algorithem 1: XY for DeLuna
  ### algorithem 2: VanderVee and Shipster "
  if(discretize==TRUE)
  {df<-discret(df)}
  cg_index<-names(which(unlist(lapply(df,class))=="integer"))
  for (i in cg_index){
    df[,i]<-as.numeric(df[,i])
  }
  cov_index<-grep("X",names(df))
  X<-df[, cov_index]
  Y<-df[,"Y"]
  T <- df$T
  if(is.null(X)) return(NULL)
  if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
  dat <- data.frame(cbind(X,Y,T))
  XsubsT<-data.frame(X,T)
  XT<-bnlearn::learn.nbr(XsubsT, "T", "si.hiton.pc",alpha=alpha)
  Xsubs0 <- X[dat$T == 0,]; Xsubs1 <- X[dat$T == 1,]
  Y0 <- Y[dat$T == 0]; Y1 <- Y[dat$T == 1]
  X0 <- bnlearn::learn.nbr(x = cbind(Xsubs0, Y0), node = "Y0", method = "si.hiton.pc", alpha = alpha)
  X1 <- bnlearn::learn.nbr(x = cbind(Xsubs1, Y1), node = "Y1", method = "si.hiton.pc", alpha = alpha)
  XY <- sort(union(X0, X1))
  if (algorithem=="1"){
    if(is.null(XY)) return(NULL)
    if(length(XY) == 0) return(NULL)
    else return (XY)}
  else{
    XTY<-sort(union(XT, XY))
    if(is.null(XTY)) return(NULL)
    if(length(XTY) == 0) return(NULL)
    return(XTY)}
}

## Reference paper: Kalisch et al. (2012)
## Reference Book: Pratical approaches to causal relationship exploration (2015) Perason
## Author:Rui Lu

bn_sel_pc <- function(df,alpha=0.05,corMethod = "standard", algorithem="1") {
  ### Change the interger value into numeric for further analysis
  ### alpha value is for mutual information based tests
  ### algorithem 1: XY for DeLuna
  ### algorithem 2: VanderVee and Shipster
  ### corMethod: a string determining the method for correlation estimation
  cov_index<-grep("X",names(df))
  X<-as.matrix(df[, cov_index])
  Y<-df[,"Y"]
  T <- ifelse(as.numeric(df$T)==1,0,1)
  if(is.null(X)) return(NULL)
  if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
  dat <- data.frame(cbind(X,Y,T))
  XsubsT<-data.frame(X,T)
  XT_model<- pcalg::pcSelect(T,X,corMethod =corMethod, alpha=alpha)
  XT<-names(XT_model$G)[ which(XT_model$G==TRUE)]
  Xsubs0 <- X[dat$T == 0,]; Xsubs1 <- X[dat$T == 1,]
  Y0 <- Y[dat$T == 0]; Y1 <- Y[dat$T == 1]
  X0_model <-pcalg::pcSelect(Y0,Xsubs0,corMethod=corMethod,alpha=alpha)
  X0<-names(X0_model$G)[ which(X0_model$G==TRUE)]
  X1_model <-pcalg::pcSelect(Y1, Xsubs1,corMethod=corMethod,alpha=alpha)
  X1<-names(X1_model$G)[ which(X1_model$G==TRUE)]
  XY <- sort(union(X0, X1))
  if (algorithem=="1"){
    if(is.null(XY)) return(NULL)
    if(length(XY) == 0) return(NULL)
    else return (XY)}
  else{
    XTY<-sort(union(XT, XY))
    if(is.null(XTY)) return(NULL)
    if(length(XTY) == 0) return(NULL)
    return(XTY)}
}


## MMHC for covariate selection
bn_sel_mmhc<-function(df, alpha=alpha,discretize=FALSE,test="mi-g",algorithem="1"){
  ## Continous : mi-g
  # Change the interger value into numeric for further analysis
  ### alpha value is for mutual information based tests
  ### algorithem 1: XY for DeLuna
  ### algorithem 2: VanderVee and Shipster
  ### corMethod: a string determining the method for correlation estimation
if(discretize==TRUE)
{df<-discret(df)}
cg_index<-names(which(unlist(lapply(df,class))=="integer"))
for (i in cg_index){
df[,i]<-as.numeric(df[,i])}
cov_index<-grep("X",names(df))
X<-df[, cov_index]
Y<-df[,"Y"]
T <- df$Treat
if(is.null(X)) return(NULL)
if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
bmT<-matrix(c(rep("T",dim(X)[[2]]),names(X)),ncol=2)
blacklistT<-data.frame(bmT)
names(blacklistT)<-c("from","to")
XsubsT<-data.frame(X,T)
XsubsT$T<-as.numeric(XsubsT$T)
model_t<-bnlearn::mmhc(XsubsT,blacklist=blacklistT,restrict.args=list(test=test ,alpha=alpha))
XT<-model_t$nodes$T$mb
bmY<-matrix(c(rep("Y",dim(X)[[2]]),names(X)),ncol=2)
blacklistY<-data.frame(bmY)
names(blacklistY)<-c("from","to")
dat_all<-data.frame(X,T,Y)
dat_y<-data.frame(X,Y)
Xsubs0 <- dat_y[dat_all$T == 0,]
Xsubs1 <- dat_y[dat_all$T == 1,]
model_y0<-bnlearn::mmhc(Xsubs0,blacklist=blacklistY, restrict.args=list(test=test,alpha=alpha))
X0<-model_y0$nodes$Y$mb
model_y1<-bnlearn::mmhc(Xsubs1,blacklist=blacklistY, restrict.args=list(test=test,alpha=alpha))
X1<-model_y1$nodes$Y$mb
XY<-unique(c(X1,X0))
if (algorithem=="1"){
  if(is.null(XY)) return(NULL)
  if(length(XY) == 0) return(NULL)
  else return (XY)}
else{
  XTY<-sort(union(XT, XY))
  if(is.null(XTY)) return(NULL)
  if(length(XTY) == 0) return(NULL)
  return(XTY)}
 }






