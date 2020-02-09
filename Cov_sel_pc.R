############################################################################################
##                  Covariate selection based on fast PC algorithm                     ##
############################################################################################
## Reference paper: Kalisch et al. (2012)
## Reference Book: Pratical approaches to causal relationship exploration (2015) Perason
## Author:Rui Lu 

bn_sel <- function(df,alpha=0.05,corMethod = "standard", algorithem="1") {
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



