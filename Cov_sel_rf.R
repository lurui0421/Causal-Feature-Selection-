############################################################################################
##                      Covariate selection based on Random Forest                       ##
############################################################################################
### Author: Rui Lu
### Reference : Bryan Keller 2019
### Relevant packages
install.packages("ranger")
library(ranger)
install.packages("partykit")
library(partykit)

### Common variable importance function
RF_vimp <- function(x, y,
                    ntree = 200,
                    scale = FALSE,
                    mtry = if (!is.null(y) && !is.factor(y))
                    max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x)))) {
  if(is.null(x)) return(NULL)
  if(dim(x)[2] == 0) return(NULL)
  if(sum(is.na(x), is.na(y)) > 0) stop("No missing values are allowed.")
  rf <- ranger::ranger(formula = Y ~ ., data = data.frame(cbind(x,"Y" = y)),
                       num.trees = ntree,
                       mtry = mtry,
                       importance = "permutation",
                       scale.permutation.importance = scale)
  return(rf$variable.importance)
}


# Hidden functions for conditional importance

formula_func <- function(var_names, index) {
  fmla <- paste0(var_names[index], " ~ ", paste0(var_names[-index], collapse = " + "))
  fmla }

subs_func <- function(vec, oob_ind) { vec[which(oob_ind == 1)] }
seqs_tab <- function(tab) { lapply(X = tab, FUN = seq, from = 1)}
seqs_sample <- function(sqs) { lapply(X = sqs, FUN = sample, replace = FALSE)}

mspe_func <- function(orig, newd) { mean((orig - newd)^2) }

perm_func_cond_help <- function(perm, ind, dat, nm, tree_fits_oob) {
  dat[which(tree_fits_oob[[ind]] == nm), ind] <- dat[which(tree_fits_oob[[ind]] == nm), ind][perm]
  dat }

perm_func_cond <- function(perm_lst, index, dat, tree_fits_oob) {
  nms <- names(perm_lst)
  for(i in 1:length(nms)) {
    dat <- perm_func_cond_help(perm = perm_lst[[i]],
                               ind = index,
                               dat = dat,
                               nm = nms[i],
                               tree_fits_oob = tree_fits_oob) }
  dat }

rfpred <- function(rfobj, newdat) {
  preds <- randomForest:::predict.randomForest(
    object = rfobj,
    newdata = newdat)
  preds }

tree_func <- function(fmla, dattta) {
  dattta <<- dattta
  if(!is.data.frame(dattta)) dattta <<- as.data.frame(dattta)
  nn <- dim(dattta)[1]
  treee <- rpart::rpart(formula = fmla, data = dattta,
                        control = rpart::rpart.control(maxsurrogate = 0,
                                                       xval = 0, cp = .01,
                                                       minbucket = max(10, round(nn*.05)),
                                                       maxcompete = 0))
  if(is.null(treee)) {nodes <- rep(1, times = dim(dattta)[1])} else {
    nodes <- predict(partykit::as.party(treee), type = "node") }
  nodes }

stiC <- function(rfobj, dat) { # single tree importance - conditional; outcome name in dat is assumed to be "Y"
  roc <- FALSE
  if(rfobj$type == "classification") {if(length(table(outcome)) == 2) {assign("roc", TRUE)} else {
    stop("Multicategory outcome not supported")}}
  pred_names <- names(dat[,-which(names(dat) == "Y")])
  oobs <- rfobj$oob.times
  noobs <- sum(oobs)
  dat_oobs <- dat[which(oobs == 1), -which(names(dat) == "Y")]
  oobs_otcms <- dat[which(oobs == 1), "Y"]
  npreds <- dim(dat)[2] - 1
  preds <- rfpred(rfobj = rfobj, newdat = dat_oobs)
  if(roc) {ERR <- roc_func(orig = oobs_otcms, newd = preds)} else {
    ERR <- mspe_func(orig = oobs_otcms, newd = preds)}
  forms <- lapply(X = 1:length(pred_names), FUN = formula_func, var_names = pred_names)
  # Get group classification by terminal node membership
  tree_fits <- lapply(X = forms, FUN = tree_func, dattta = dat[, pred_names])
  # Subset oob
  tree_fits_oob <- lapply(X = tree_fits, FUN = subs_func, oob_ind = oobs)
  # Tables of terminal node frequencies
  tabs <- lapply(X = tree_fits_oob, FUN = table)
  seqs <- lapply(X = tabs, FUN = seqs_tab)
  perms <- lapply(X = seqs, FUN = seqs_sample)

  perm_dat <- dat_oobs
  lngth <- dim(perm_dat)[2]
  for(j in 1:lngth) {
    perm_dat <- perm_func_cond(perm_lst = perms[[j]],
                               index = j,
                               dat = perm_dat,
                               tree_fits_oob = tree_fits_oob) }
  perm_dat_list <- vector(length = lngth, mode = "list")
  for(i in 1:lngth) {
    perm_dat_list[[i]] <- dat_oobs
    perm_dat_list[[i]][,i] <- perm_dat[,i] }

  imp_preds <- lapply(X = perm_dat_list, FUN = rfpred, rfobj = rfobj)
  if(roc) {ERRs <- sapply(X = imp_preds, FUN = roc_func, orig = oobs_otcms)} else {
    ERRs <- sapply(X = imp_preds, FUN = mspe_func, orig = oobs_otcms) }
  if(roc) {imp_raw <- ERR - ERRs} else {
    imp_raw <- ERRs - ERR }
  imp_raw
}

oaiC <- function(rf_list, dat) { # overall importance - Conditional; outcome name must be "Y"
  imps <- lapply(X = rf_list, FUN = stiC, dat = dat)
  imps <- matrix(unlist(imps), length(imps), length(imps[[1]]), byrow = TRUE)
  raw_imps <- apply(imps, 2, mean)
  scaled_imps <- apply(imps, 2, mean)/(apply(imps, 2, sd) / sqrt(dim(imps)[1]))
  out <- cbind(raw_imps, scaled_imps)
  colnames(out) <- c("raw", "scaled")
  out }


### Conditional random forest variable importance
CRF_vimp <- function(x, y,
                     scale = FALSE,
                     ntree = 50,
                     mtry = if (!is.null(y) && !is.factor(y))
                       max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x)))) {
  if(is.null(x) | length(x) == 0) return(NULL)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(p == 1) {
    crf <- RF_vimp(x = x, y = y, ntree = ntree)
    names(crf) <- names(x); return(crf) }
  dat_whole <- data.frame(cbind(x,"Y" = y))
  RFs <- replicate(n = ntree, expr = randomForest::randomForest(x = x, y = y, keep.forest = TRUE, ntree = 1), simplify = FALSE)
  imps <- oaiC(rf_list = RFs, dat = dat_whole)
  rownames(imps) <- names(x)
  if(scale) return(imps[,2])
  if(!scale) return(imps[,1])
}


### Permutation based covariate importance

RFVS_perm_trad <- function(x, y, perms = 20,
                           ntree = 200,
                           pctl = .95,
                           scale = FALSE,
                           plot = FALSE) {
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
  imps <- lapply(X = perm_dat_list, FUN = RF_vimp, y = y, ntree = ntree, scale = scale)
  imp <- imps[[1]]
  imps <- matrix(unlist(imps[2:length(perm_dat_list)]), nrow = length(imps[[2]]), ncol = perms, byrow = FALSE)
  pctls <- apply(X = imps, MARGIN = 1, FUN = quantile, probs = pctl)
  ord <- order(imp, decreasing = TRUE)
  imp_ord <- imp[ord]; pctls_ord <- pctls[ord]
  if(plot) {
    par(mar = c(5.1, 6.5, 1.1, 2.1))
    hts <- imp_ord[length(imp_ord):1]
    hts[which(hts < 0)] <- 0
    barplot(height = hts, horiz = TRUE,
            las = 1, xlab = "Random Forest Permutation Importance",
            main = "", col = "skyblue",
            cex.names=.8, xlim = c(0, max(c(imp,pctls))))
    barplot(height = pctls_ord[length(pctls_ord):1], horiz = TRUE,
            col = rgb(1, 0, 0, alpha=0.2),
            add = TRUE) }
  out <- cbind(imp_ord, pctls_ord)
  out <- data.frame(out)
  rownames(out) <- names(imp_ord)
  colnames(out) <- c("Imp", "Pctl")
  return(out)
}

### Conditional based permutation covariate importance

RFVS_perm_cond <- function(x, y, perms = 10,
                           ntree = 50,
                           pctl = .95,
                           scale = FALSE,
                           plot = FALSE,
                           permcores = 1) {
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
  imps <- parallel::mclapply(X = perm_dat_list, FUN = CRF_vimp, y = y, ntree = ntree, scale = scale, mc.cores = permcores)
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



### Variable selection use random forest XT

RF_sel_XT <- function(df,perms = 20,ntree = 200,pctl = .95,scale = FALSE)
{cov_index<-grep("X",names(df))
  X<-df[, cov_index]
  Y<-df[,"Y"]
  T <- factor(df$Treat)
  if(is.null(X)) return(NULL)
  if(dim(X)[2] == 0) return(NULL)
  if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
  if(length(table(T)) != 2) stop("Exposure variable T must have exactly two categories.")
  levels(T) <- 0:1
  Xsubs0 <- X[df$T == 0,]; Xsubs1 <- X[df$T == 1,]
  XT<-RFVS_perm_trad(x=X,y=T, ntree = ntree,perms = perms, pctl = pctl,scale = scale)
  XT <- rownames(XT)[which(XT[,1] > 0)]
  if(is.null(XT)) return(NULL)
  if(length(XT) == 0) return(NULL)
  else return (XT)}

### Variable selection use random forest XY

RF_sel_XY <- function(df,perms = 20,ntree = 200,pctl = .95,scale = FALSE) {
  cov_index<-grep("X",names(df))
  X<-df[, cov_index]
  Y<-df[,"Y"]
  T <- factor(df$Treat)
  if(is.null(X)) return(NULL)
  if(dim(X)[2] == 0) return(NULL)
  if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
  if(length(table(T)) != 2) stop("Exposure variable T must have exactly two categories.")
  levels(T) <- 0:1
  Xsubs0 <- X[df$T == 0,]; Xsubs1 <- X[df$T == 1,]
  Y0 <- Y[T == 0]; Y1 <- Y[T == 1]
  X0 <- RFVS_perm_trad(x = Xsubs0, y = Y0, ntree = ntree,
                       perms = perms, pctl = pctl,
                       scale = scale)
  X1 <- RFVS_perm_trad(x = Xsubs1, y = Y1, ntree = ntree,
                       perms = perms, pctl = pctl,
                       scale = scale)
  X0 <- rownames(X0)[which(X0[,1] > 0 & X0[,1] >= X0[,2])]
  X1 <- rownames(X1)[which(X1[,1] > 0 & X1[,1] >= X1[,2])]
  XY<-union(X1,X0)
  if(is.null(XY)) return(NULL)
  if(length(XY) == 0) return(NULL)
  else return (XY)}

### variable selection use conditional random forest XY

RF_sel_XY_CON <- function(df,
                          perms = 10,
                          ntree = 50,
                          pctl = .95,
                          scale = FALSE,
                          permcores = 1) {
  cov_index<-grep("X",names(df))
  X<-df[, cov_index]
  Y<-df[,"Y"]
  T <- factor(df$Treat)
  if(is.null(X)) return(NULL)
  if(dim(X)[2] == 0) return(NULL)
  if(sum(is.na(X), is.na(Y), is.na(T)) > 0) stop("No missing values are allowed.")
  if(length(table(T)) != 2) stop("Exposure variable T must have exactly two categories.")
  T <- factor(T); levels(T) <- 0:1
  Xsubs0 <- X[T == 0, , drop = FALSE]; Xsubs1 <- X[T == 1, , drop = FALSE]
  Y0 <- Y[T == 0]; Y1 <- Y[T == 1]
  X0 <- RFVS_perm_cond(x = Xsubs0, y = Y0, ntree = ntree,
                       perms = perms, pctl = pctl,
                       scale = scale, permcores = permcores)
  X1 <- RFVS_perm_cond(x = Xsubs1, y = Y1, ntree = ntree,
                       perms = perms, pctl = pctl,
                       scale = scale, permcores = permcores)

  X0 <- rownames(X0)[which(X0[,1] > 0 & X0[,1] >= X0[,2])]
  X1 <- rownames(X1)[which(X1[,1] > 0 & X1[,1] >= X1[,2])]
  XY <- sort(union(X0, X1))
  if(is.null(XY)) return(NULL)
  if(length(XY) == 0) return(NULL)
  return(XY) }




