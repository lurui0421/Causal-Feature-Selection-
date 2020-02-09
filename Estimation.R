######################################################################################
##                         Treatment effect estimation                              ##
#######################################################################################
##: This is the treatment effect estimatio based on various algirthems
##: These are embedded covariate selection method for causal inference
###: Author : Rui Lu

### 1: Propensity score matching
##Pre_request functions and packages for esitmation of propensity score

install.packages("Matching")
library(Matching)

logit<-function(PS){
  logitPS<- log(PS /(1 - PS))
  return(logitPS)
}


est_ps_match<-function(df,M=1,cov_sel,preprocessing=TRUE)
{
  ## Data is particular data set used for estimation
  ## cov_sel is the selected covariates
  ## M is the 1:1 or 1:2 matching

  ## Step 1: Estimate the propensity score using GBM
if (preprocessing==TRUE)
{df<-df[,c(cov_sel,"Treat","Y")]}
measure.var<-"Treat"
sel.var<-names(df[,grep("X",names(df))])
formula<-as.formula(paste(measure.var, paste(sel.var, collapse=" + "), sep=" ~ "))
twang.fit<- twang::ps(formula = formula,data=df,
                      n.trees=1000, interaction.depth=2,
                      shrinkage=0.01,stop.method="es.mean",
                      estimand = "ATE")
p.score<-twang.fit$ps
names(p.score)<-NULL
trt<-df$Treat
### Step 2: Propensity score matching
psmatch_ATT<-Match(Tr=trt,M=M,estimand ="ATT",X=sapply(p.score,logit),replace=FALSE,caliper=.2)
## Estimate ATT
matched_ATT<-df[unlist(psmatch_ATT[c("index.treated","index.control")]), ]
fit_ATT<-lm(Y~Treat,data=matched_ATT)
result_ATT<-summary(fit_ATT)
ATT<-result_ATT$coefficients[[2]]
psmatch_ATC<-Match(Tr=trt,M=M,estimand ="ATC",X=sapply(p.score,logit),replace=FALSE,caliper=.2)
## Estimate the ATC
matched_ATC<-df[unlist(psmatch_ATC[c("index.treated","index.control")]), ]
fit_ATC<-lm(Y~Treat,data=matched_ATC)
result_ATC<-summary(fit_ATC)
ATC<-result_ATC$coefficients[[2]]
ATE<-mean(ATC,ATT)
return(ATE)}

data_lin.100<-readRDS(file=file.choose())
cov_sel.dlasso.100<-readRDS(file=file.choose())
ps_100.lin.dlasso<-list()
for (i in 1:100)
{
ps_100.lin.dlasso[[i]]<-est_ps_match(data_lin.100[[i]][[1]],cov_sel=cov_sel.dlasso.100[[i]])
}


### 2: Genetic Matching
##Pre_request functions and packages for esitmation of propensity score

install.packages('rgenoud')
library(rgenoud)
install.packages("Matching")
library(Matching)

est_gen_match<-function(df,M=1,cov_sel,pre_processing=TRUE)
{ ### df represent the data set
  ## M is the number of matched samples
  ## cov_sel is the name of covariate selected
  ## pre_processing means wheather covariate selection is conducted or not
  if (pre_processing==TRUE)
    {df<-df[,c(cov_sel,"Treat","Y")]}
    ### Step 1: Genatic Matching to get weights
    X<-df[,grep("X",names(df))]
    BalanceMat <-X
    Treat<-df$Treat
    Y<-df$Y
    ### Estimate ATC
    genout_ATC <- GenMatch(Tr=Treat,X=X, BalanceMatrix=BalanceMat, estimand="ATC", M=M,
                       pop.size=16, max.generations=10, wait.generations=1)
    mout_ATC <- Match(Y=Y, Tr=Treat, X=X, estimand="ATC", Weight.matrix=genout_ATC)
    ATC<-as.vector(mout_ATC$est)

    ### Estimate ATT
    genout_ATT <- GenMatch(Tr=Treat,X=X, BalanceMatrix=BalanceMat, estimand="ATT", M=M,
                           pop.size=16, max.generations=10, wait.generations=1)
    mout_ATT <- Match(Y=Y, Tr=Treat, X=X, estimand="ATT", Weight.matrix=genout_ATT)
    ATT<-as.vector(mout_ATT$est)

    ### Step 3: Estimate ATE
    ATE<-mean(ATT,ATC)
    return(ATE)}


## 3: BART : with propensity score justification
## Reference: Hill,2012

install.packages("devtools")
devtools::install_github("vdorie/bartCause")
library(bartCause)
XY=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")
est_bart<-function(df,cov_sel=c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10"),
                   method.rsp='bart',
                   method.trt='none',
                   pre_processing=TRUE,
                   estimand='ate')
{### Vinilla BART function
  if (pre_processing==TRUE)
  {df<-df[,c(cov_sel,"Treat","Y")]}
  cov_index<-grep("X",names(df))
  X<-as.matrix(df[,cov_index])
  Y<-df[,"Y"]
  T<-df$Treat
  fit<-bartc(Y,T,X,method.rsp=method.rsp,method.trt= method.trt,estimand='ate')
  result<-summary(fit)
  return(result$estimates)}

### 4: C-TMLE_scalable,Lasso-CTMLE
### : Reference : Ju et.al, 2018
install.packages('ctmle')
library(ctmle)
est_ctmle<-function(df,method="ctmle",cov_sel=c("X3","X4","X5","X6","X7","X8","X9","X10"),pre_processing=TRUE)
{ ### Method =lasso means using CTMLE with lasso for propensity score estimation.other wise scalable CTMLE
  ##  Pre_process means if pre_rpocessing of covarite estimation is included
  ##  cov_set is the selected covariate set
  if(pre_processing==TRUE)
  {
  cov_set<-c(cov_sel,"Treat","Y")
  df<-df[,cov_set]}
  N<-dim(df)[1]
  ## Q matrix with initialization
  Q <- cbind(rep(mean(df$Y[df$Treat == 0]), N),
             rep(mean(df$Y[df$Treat == 1]), N))
  cov_index<-grep("X",names(df))
  ## W matrix
  W<-as.matrix(df[,cov_index])
  Treat<-df$Treat
  Y<-df$Y
  if (method=="lasso")
  {glmnet_fit<-cv.glmnet(Treat,x=W,family='binomial',nlambda=20)
  lambdas <- glmnet_fit$lambda[(which(glmnet_fit$lambda==glmnet_fit$lambda.min)):length(glmnet_fit$lambda)]
  ctmle_fit<- ctmleGlmnet(Y = Y,
                          A = Treat,
                          W = W,
                          Q = Q,
                          lambdas=lambdas,
                          ctmletype=1,
                          family="gaussian",
                          gbound=0.025,
                          V=5)}
  else {ctmle_fit<- ctmleDiscrete(Y = Y,
                            A = Treat,
                            W = data.frame(W),
                            Q = Q,
                            preOrder = TRUE,
                            order = rev(1:length(cov_index)),
                            detailed = TRUE)}
  result<-summary(ctmle_fit)
  return(result[[1]])
}

### 5: Super learning
install.packages("gam")
library(gam)
install.packages("tmle")
library(tmle)
install.packages("ranger")
library(ranger)
install.packages("dbarts")
library(dbarts)

est_super_learn<-function(df,cov_sel=c("X3","X4","X5","X6","X7","X8","X9","X10"),pre_processing=TRUE)
{
  if(pre_processing==TRUE)
  {cov_set<-c(cov_sel,"Treat","Y")
  df<-df[,cov_set]}
  A<-df$Treat
  Y<-df$Y
  W<-df[,grep("X",names(df))]
  SL.TMLER <- tmle(Y=Y, A=A, W=W,
                       family="gaussian",
                       Q.SL.library = c("SL.glm",
                                        "tmle.SL.dbarts2",
                                        "SL.glmnet"),
                       g.SL.library = c("SL.glm",
                                        "tmle.SL.dbarts.k.5",
                                        "SL.gam",
                                        "SL.ranger"))

  ATE<-SL.TMLER$estimates$ATE$psi
  return(ATE)}





