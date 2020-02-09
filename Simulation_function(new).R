################################################################################################
#####                            Data Generation Functions                                   ###
################################################################################################
##: This is the data generation process
###:  Author : Rui Lu
### 1: Preprations and Functions
scale_var<-function(var){
## This function is used to perform standard scalling for variables
  var<-as.numeric(var)
  new_var<-as.vector(scale(var,center = TRUE, scale = TRUE))
  return (new_var)}


## Get covariate and noise set
covariates_set<-readRDS(file=file.choose())
names(covariates_set)<-paste("X",1:34,sep="")
noise_set<-readRDS(file=file.choose())
names(noise_set)<-paste("X",35:3908,sep="")

### 2: Data generation function
DGP<-function(sample_index,
              dimension="low",
              treatment="linear",
              response="linear",
              redudant=FALSE,
              covariates_set=covariates_set,
              noise_set=noise_set,
              heterogenity=FALSE,
              scale=TRUE,
              noise=0.5)
{ ## Scale the variable
  if(scale==TRUE)
  {
  covariates_set<-data.frame(apply(covariates_set,2,scale_var))
  noise_set<-data.frame(apply(noise_set,2,scale_var))}
  ## Factor one : sample size
  n_size<-length(sample_index)
  sub_sample<-covariates_set[sample_index,]
  ## Factor two : dimesionality
  if (dimension=="low")
  {sub_sample_data<-data.frame(sub_sample)}
  else if(dimension=="medium")
  {
  p_size<-0.5*n_size
  noise_sample<-noise_set[sample_index,1:(p_size-34)]
  sub_sample_data<-data.frame(sub_sample,noise_sample)}
  else{
    p_size<-1.5*n_size
    noise_sample<-noise_set[sample_index,1:(p_size-34)]
    sub_sample_data<-data.frame(sub_sample,noise_sample)}
  ## Factor three: add redudant variables (muti-collinerity)
  if (redudant==TRUE){
    sub_sample_data$X_re_1<-2*sub_sample_data$X5+rnorm(n_size)
    sub_sample_data$X_re_2<-2*sub_sample_data$X6+rnorm(n_size)
    sub_sample_data$X_re_3<-2*sub_sample_data$X7+rnorm(n_size)
  }
  else
  {sub_sample_data<-sub_sample_data}
  attach(sub_sample_data)
  ## Factor four: Response and treatment surfaces
  if(treatment=='linear'){
    Z<-X1+X2+1.5*X3+1.5*X4+X5+X6+0.5*X7+0.5*X8+rnorm(n_size,0,noise)
    p<-abs(round(round(1/(1+exp(0.05+Z)),3)-0.05,2))
    Treat<-rbinom(n_size,1,prob=p)}
  else
  {
  Z<-X1+X2+1.5*I(X3>mean(X3))+X4*X5+rnorm(n_size,0,noise)
  p<-abs(round(round(1/(1+exp(-1.2+Z)),3)-0.1,2))
  Treat<-rbinom(n_size,1,prob=p)}
  if(response=='linear'){
    f<-0.5*X3+X4+X5+X6+1.5*X7+1.5*X8+X9+X10
    Y0<-f+rnorm(n_size,0,noise)
    Y1<-5+f+rnorm(n_size,0,noise)}
  else {
    f<-X3*X4+X5+1.5*X6+1.5*sin(X7)+2*I(X8>mean(X8))+2*X9*X10
    Y0<-f+rnorm(n_size)
    Y1<-5+f+rnorm(n_size)}
  if(heterogenity==TRUE){
    f<-0.5*X3+0.5*X4+X5+X6+1.5*X7+1.5*X8+X9+X10
    Y0<-f+I(X11>mean(X11))+rnorm(n_size,0,noise)
    Y1<-5+f+X12*X13+I(X14>mean(X14))+X15*X16+rnorm(n_size,0,noise)
  }
  Y<-Treat*Y1+(1-Treat)*Y0
  sample_data<-data.frame(sub_sample_data,Treat=Treat,Y0=Y0,Y1=Y1,Y=Y)
  detach(sub_sample_data)
  ###: Measurements
  ## (1):Nonlinearity: response surface and assignments
  fit_res<-lm(Y0~X3+X4+X5+X6+X7+X8+X9+X10,data=sample_data)
  fit_tre<-lm(Z~X1+X2+X3+X4+X5+X6+X7+X8,data=sample_data)
  R_2_response<-round(summary(fit_res)$r.squared,4)
  R_2_treat<-round(summary(fit_tre)$r.squared,4)
  ##(2):Treatment vs. control rate
  treat_control_rate<-mean(sample_data$Treat)
  ##(3):Alignments
  Alignment<-abs(cor(p,Y))
  ##(4):Signal to noise ratio
  SNR<-var(f)/noise
  PVE<-(SNR)/(1+SNR)
  ##(5):Treatment heterogenity
  THM<-sd(Y1-Y0)/sd(Y)
  ##(6):Overall measurements
  measurements<-list(R_2_response= R_2_response, R_2_treat= R_2_treat,
                     treat_control_rate=treat_control_rate,
                     Alignment=Alignment,SNR=SNR,PVE=PVE,THM=THM)
  result<-list(data=sample_data,measurements=measurements)
  return(result)}


### sub-sampling function
Sample_data<-function(  n_size,
                        covariates_set,
                        noise_set,
                        dimension="low",
                        treatment="linear",
                        response="linear",
                        redudant=FALSE,
                        heterogenity=FALSE,
                        scale=TRUE,
                        noise=1)
{## Step 1: Sample Size
sample_index<-sample(1:nrow(covariates_set),n_size,replace=FALSE)
data<-DGP(sample_index,
          dimension=dimension,
          treatment=treatment,
          response= response,
          redudant=redudant,
          covariates_set=covariates_set,
          noise_set=noise_set,
          heterogenity=heterogenity,
          scale=scale,
          noise=noise)
### Step 2: Generate data set
if(dimension=="medium")
{ n_size<-n_size/2
}
else if (dimension=="high")
{n_size<-1.5*n_size}
else
  {n_size<-34}
n_x<-n_size+1
n_y<-n_size+4
data_set<-data[[1]][,c(1:n_x,n_y)]
## Step 3: Quality control measurements
measure<-data[[2]]
return(list(data_set,measure))}

### Generate data set with various sample sizes with linear

data_100.lin_low<-replicate(100,Sample_data(100,dimension="low",covariates_set,noise_set),simplify = FALSE)
setwd("/Users/lurui/Desktop")
saveRDS(data_100.lin_low,"data_100.lin_low.RDS")
data_500.lin_low<-replicate(100,Sample_data(500,dimension="low",covariates_set,noise_set),simplify = FALSE)
saveRDS(data_500.lin_low,"data_500.lin_low.RDS")
data_2000.lin_low<-replicate(100,Sample_data(2000,dimension="low",covariates_set,noise_set),simplify = FALSE)
saveRDS(data_2000.lin_low,"data_2000.lin_low.RDS")

data_100.lin_medium<-replicate(100,Sample_data(100,dimension="medium",covariates_set,noise_set),simplify = FALSE)
saveRDS(data_100.lin_medium,"data_100.lin_medium.RDS")

data_500.lin_medium<-replicate(100,Sample_data(500,dimension="medium",covariates_set,noise_set),simplify = FALSE)
saveRDS(data_500.lin_medium,"data_500.lin_medium.RDS")

data_2000.lin_medium<-replicate(100,Sample_data(2000,dimension="medium",covariates_set,noise_set),simplify = FALSE)
saveRDS(data_2000.lin_medium,"data_2000.lin_medium.RDS")

data_100.lin_high<-replicate(100,Sample_data(100,dimension="high",covariates_set,noise_set),simplify = FALSE)
saveRDS(data_100.lin_high,"data_100.lin_high.RDS")

data_500.lin_high<-replicate(100,Sample_data(500,dimension="high",covariates_set,noise_set),simplify = FALSE)
saveRDS(data_500.lin_high,"data_500.lin_high.RDS")

data_2000.lin_high<-replicate(100,Sample_data(2000,dimension="high",covariates_set,noise_set),simplify = FALSE)
saveRDS(data_2000.lin_high,"data_2000.lin_high.RDS")

data_2000.lin_high[[1]][[1]]

################################################################################################
#####                            Simulation Quality Measurements                          ######
################################################################################################
## Measurement to ensure the quality of the simulation
#1: R square for nonlinearity
#2: Treat vs. control rate
#3: Alignment (contrl the range of the alignment)(as measurements)
#4: Signal to noise ratio. (as maro-measurements)
#5: Measure of the treatment heterogenity

################################################################################################
#####                        Simple Simulation function                                     ####
################################################################################################

### Keller Test function interaction
dgpINT <- function(n, noise = 12) {
  sig <- diag(8)
  dich <- function(vec) { vec >= mean(vec) }
  X <- mvtnorm::rmvnorm(n = n, sigma = sig)
  X[,1] <- dich(X[,1]); X[,6] <- dich(X[,6]); X[,7] <- dich(X[,7])
  coefs_ps <- t(t(c(log(2), 0, 0, log(2), 0, log(2), 0, 0)))
  coefs_oc <- t(t(c(1, 0, 0, 0, 1, 0, 2, 0)))
  ps <- (1 + exp( -(0 + X %*% coefs_ps + 1*X[,2]*X[,3]) ) )^-1
  T <- rbinom(n = n, size = 1, prob = ps)
  Y <- 0 + X %*% coefs_oc + 1*X[,2]*X[,3] + 0.5*T + rnorm(n, sd = 1)
  dch <- floor(noise/2)
  Xaug <- mvtnorm::rmvnorm(n = n, sigma = diag(noise))
  Xaug[,1:dch] <- apply(Xaug[,1:dch], 2, dich)
  dat <- data.frame(cbind(X, Xaug, Y, T))
  dat$T <- factor(dat$T)
  names(dat) <- c(paste0("X", 1:(8 + noise)), "Y", "Treat")
  dat
}

test_data<-dgpINT(1000)
### Keller Test Function linear
dgpLINEAR <- function(n, noise = 12) {
  sig <- diag(8)
  dich <- function(vec) {vec >= mean(vec) }
  X <- mvtnorm::rmvnorm(n = n, sigma = sig)
  X[,1] <- dich(X[,1]); X[,6] <- dich(X[,6]); X[,7] <- dich(X[,7])
  coefs_ps <- t(t(c(log(2), log(2), log(2), log(2), 0, log(2), 0, 0)))
  coefs_oc <- t(t(c(1, 1, 1, 0, 1, 0, 2, 0)))
  ps <- (1 + exp( -(0 + X %*% coefs_ps) ) )^-1
  T <- rbinom(n = n, size = 1, prob = ps)
  Y <- 0 + X %*% coefs_oc + 0.5*T + rnorm(n, sd = 1)
  dch <- floor(noise/2)
  Xaug <- mvtnorm::rmvnorm(n = n, sigma = diag(noise))
  Xaug[,1:dch] <- apply(Xaug[,1:dch], 2, dich)
  dat <- data.frame(cbind(X, Xaug, Y, T))
  dat$T <- factor(dat$T)
  names(dat) <- c(paste0("X", 1:(8 + noise)), "Y", "Treat")
  dat
}

test_data<-dgpLINEAR(100)





