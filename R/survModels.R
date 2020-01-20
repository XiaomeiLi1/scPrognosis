## This file is structured as follows:
## 1) required libraries are loaded,
## 2) required learners are wrapped in the R package mlr,
## 3) measures are defined,
## 4) benchmark experiments are defined for each dataset (estimated mean perfor
##    mances are computed after benchmarking),
## 5) and boxplot is generated as a .tex file, which needs to be compiled to ob-
##    tain the final figure presented in the paper.

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(mlr)
library(caret)
library(irace)# for iterated F-racing
library(Hmisc)# for rcorr.cens
library(Rmisc)
library(BBmisc)
library(checkmate)
library(survival)
library(ggplot2)
library(mlr)
library(plyr)
require(tikzDevice)
# ------------------------------------------------------------------------------
# Wrap learners into the mlr package
# ------------------------------------------------------------------------------

#--- TuneWrapper for cox model
wrapperPHModell <- function(center = TRUE, scale = TRUE,
                            ...) {
  #task <- makeSurvTask(data = data, target = target)
  lrn <- makePreprocWrapperScale(makeLearner(cl = "surv.coxph", ...),
                                 center = center, scale = scale)
  return(lrn)
}

# ------------------------------------------------------------------------------
# Wrap measures into the mlr package
# ------------------------------------------------------------------------------
#' my.ci.fun computes the C-index for "survivalsvmprediction" objects.
#'
#' @param task survival task used to fit the model.
#'        This will be later used by "makeMeasure" from the pakage "mlr".
#' @param model model of class "survivalsvm".
#'        This will be later used by "makeMeasure" from the pakage "mlr".
#' @param pred The predictions.
#' @param feats will be later used by "makeMeasure" from the pakage "mlr".
#' @param extra.args will be later used by "makeMeasure" from the pakage "mlr".
#'
#' @return The C-Index
#'
my.ci.fun <- function(task, model, pred, feats, extra.args) {
  myci = rcorr.cens(x = getPredictionResponse(pred),
                    S = getPredictionTruth(pred))
  return(myci["C Index"])
}
# constructs the C-index measure for survivalsvmprediction objects.
c.i <- makeMeasure(
  id = "ci", name = "C-Index",
  properties = c("surv"),
  minimize = FALSE, best = 1, worst = 0,
  fun = my.ci.fun
)

# constructs the C-index measure for PH Model.
ph.cindex = makeMeasure(id = "ci", minimize = FALSE, best = 1, worst = 0,
                        properties = c("surv", "req.pred", "req.truth"),
                        name = "C-Index",
                        fun = function(task, model, pred, feats, extra.args) {
                          requirePackages("Hmisc", default.method = "load")
                          resp = pred$data$response
                          if (anyMissing(resp))
                            return(NA_real_)
                          s = Surv(pred$data$truth.time, pred$data$truth.event)
                          Hmisc::rcorr.cens(-1 * resp, s)[["C Index"]]
                        }
)

# --- Log-rank statistic ######
#
#' @param t1 Observations in group 1
#' @param d1 Predictions in group 1
#' @param t2 Observations in group 2
#' @param d2 Predictions in group 2
#'
#' @return log-rank statistik
logrank <- function(t1, d1, t2, d2){
  t1.ord <- order(t1)
  t1 <- t1[t1.ord]
  d1 <-  d1[t1.ord]
  t2.ord <- order(t2)
  t2 <- t2[t2.ord]
  d2 <-  d2[t2.ord]
  t <- c(t1, t2)
  d <- c(d1, d2)
  times <- unique(t[d == 1])
  n <- length(times)
  tel.ner <- sapply(times, function(i){
    o <- sum(t == i) # failures
    r <- sum(t >= i) # at risk
    o1 <- sum(t1 == i) # failures
    r1 <- sum(t1 >= i) # at risk
    o2 <- sum(t2 == i) # failures
    r2 <- sum(t2 >= i) # at risk
    teller <- o1 - r1*o/r
    noemer <- r2*r1*o*(r-o) / (r^2 *(r - 1))
    return(c(teller, noemer))
  })
  search.na <- colSums(tel.ner)
  na.index <- which(is.na(search.na))
  if(length(na.index) > 0){
    tel.ner <- tel.ner[, -na.index] 
  }
  chi <- sum(tel.ner[1,])^2 / sum(tel.ner[2,])
  return(chi_sq = chi)
}

#' @param task a survival task
#' @param model survival model
#' @param pred survival predictions
#' @param feats features
#' @param extra.args further arguments
#'
#' @return log-rank statistik
lr.fun <- function(task, model, pred, feats, extra.args){
  t <- getPredictionTruth(pred)[,1]
  delta <- getPredictionTruth(pred)[,2]
  u <- getPredictionResponse(pred)
  part2 <- which(u > mean(u))
  part1 <- setdiff(1:length(t), part2)
  return(logrank(t1 = t[part1], d1 = delta[part1], t2 = t[part2],
                 d2 = delta[part2]))
}

lgrk <- makeMeasure(
  id = "logrank", name = "Log-Rank",
  properties = c("surv"),
  minimize = FALSE, best = Inf, worst = 0,
  fun = lr.fun
)

# --- Hazardratio
#' hr.fun computes the hazard rate for survivalsvm objects given:
#'
#' @param task a survival task
#' @param model a survivalsvm model
#' @param pred object of class survivalsvm
#' @param feats the features
#' @param extra.args further arguments
#'
#' @return hazard rate
hr.fun <- function(task, model, pred, feats, extra.args){
  u <- getPredictionResponse(pred)
  a <- min(u)
  b <- max(u)
  u <- (u-a) / (b-a)
  mod.ph <- try(coxph(getPredictionTruth(pred) ~ u), silent = TRUE)
  if(!("list" %in% is(mod.ph))) return(1)
  return(exp(mod.ph$coef))
}

# construction of the measure
hr <- makeMeasure(
  id = "hr", name = "Hazardrate",
  properties = c("surv"),
  minimize = FALSE, best = Inf, worst = 0,
  fun = hr.fun
)

# --- Hazardrate
#' computes the hazard rate for reference models given:
#'
#' @param task a task
#' @param model a model
#' @param pred predictions
#' @param feats features
#' @param extra.args further arguments
#'
#' @return hazard rate
hrph.fun <- function(task, model, pred, feats, extra.args){
  u <- getPredictionResponse(pred)
  a <- min(u)
  b <- max(u)
  u <- (u-a) / (b-a)
  mod.ph <- try(coxph(getPredictionTruth(pred) ~ u), silent = TRUE)
  return(exp(-mod.ph$coef))
}

hrph <- makeMeasure(
  id = "hr", name = "Hazardrate",
  properties = c("surv"),
  minimize = FALSE, best = Inf, worst = 0,
  fun = hrph.fun
)

#my_seed <- .Random.seed

# #---------- Some pre-processing steps ----------------------------------------
makeTrainAndTestDataSets <- function(train.dataset = NULL, test.dataset = NULL, list = NULL){
  idx = which(rownames(train.dataset) %in% list)
  temp=train.dataset[idx,]
  expd = t(exprs(temp))
  clind = pData(temp)[,c('t.rfs',"e.rfs")]
  t1 = as.numeric(clind[,1])
  t1[t1<0] <- 0
  t2 = as.numeric(clind[,2])
  clind = data.frame(cbind(t1,t2))
  colnames(clind) = c("time", "status")
  data = cbind(expd, clind)
  train.task <- makeSurvTask(data = data, target = c("time", "status"))
  
  ####Build testing task
  temp = t(exprs(test.dataset))
  train.genes = colnames(expd)
  expd = NULL
  for (i in train.genes) {
    idx = which(colnames(temp) == i)
    if (length(idx)) expd = cbind(expd, temp[,idx])
    else expd = cbind(expd, rep(0, nrow(expd)))
  }
  colnames(expd) = train.genes
  clind = pData(test.dataset)[,c('t.rfs',"e.rfs")]
  t1 = as.numeric(clind[,1])
  t2 = as.numeric(clind[,2])
  clind = data.frame(cbind(t1,t2))
  colnames(clind) = c("time", "status")
  data = cbind(expd, clind)
  test.task <- makeSurvTask(data = data, target = c("time", "status"))
  
  return(list("train.task" = train.task, "test.task" = test.task))
}

makeCVDataSets <- function(dd = NULL, pd = NULL, list = NULL){
  idx = which(rownames(dd) %in% list)
  #print(length(idx))
  expd=t(dd[idx,])
  t1 = as.numeric(pd[,1])
  t1[is.na(t1)] <- 0
  t1[t1<0] <- 0
  t2 = as.numeric(pd[,2])
  pd = data.frame(cbind(t1,t2))
  colnames(pd) = c("time", "status")
  if(nrow(expd) == 1) {data = cbind(t(expd), pd)
  colnames(data) = c("G1","time", "status")
  }
  else data = cbind(expd, pd)
  datasets <- makeSurvTask(data = data, target = c("time", "status"))
  return(datasets)
}

makeCVData <- function(mRNA = NULL, pd = NULL){
  #dd = get(dn) 
  expd=t(mRNA)
  t1 = as.numeric(pd[,1])
  t1[is.na(t1)] <- 0
  t1[t1<0] <- 0
  t2 = as.numeric(pd[,2])
  pd = data.frame(cbind(t1,t2))
  colnames(pd) = c("time", "status")
  if(nrow(expd) == 1) {data = cbind(t(expd), pd)
  colnames(data) = c("G1","time", "status")
  }
  else data = cbind(expd, pd)
  return(data)
}

# doPHModel <- function(train.task = NULL, test.task = NULL, n.iter = 100) {
#   ci.vec <- vector("double",n.iter)
#   lgrk.vec <- vector("double",n.iter)
#   #hrph.vec <- vector("double",n.iter)
#   for (i in 1:n.iter) {
#     print(i)
#     set.seed(seed[i])
#     lrn <- makeLearner(cl = "surv.coxph")
#     rdesc  <-  makeResampleDesc(method = "CV", iters = 10L)
#     res = resample(lrn, train.task, rdesc, measures = list(ph.cindex, lgrk),show.info = FALSE)
#     temp = res$aggr
#     ####test join from ci
#     # resp = res$pred
#     # a = Hmisc::rcorr.cens(x = getPredictionResponse(resp),
#     #                       S = getPredictionTruth(resp))
#     # res.list[[i]] = list("ci" = a[["C Index"]], "mean.aggr" = temp)
#     ci.vec[i] = temp[1]
#     lgrk.vec[i] = temp[2]
#     #hrph.vec[i] = temp[3]
#   }
#   means = NULL
#   sds = NULL
#   means = c(means,mean(ci.vec,na.rm = TRUE))
#   sds = c(sds,sd(ci.vec,na.rm = TRUE))
#   means = c(means,mean(lgrk.vec,na.rm = TRUE))
#   sds = c(sds,sd(lgrk.vec,na.rm = TRUE))
#   # means = c(means,mean(hrph.vec,na.rm = TRUE))
#   # sds = c(sds,sd(hrph.vec,na.rm = TRUE))
#   #return(list("cis" = ci.vec, "lrs" = lgrk.vec, "hrs" = hrph.vec))
#   #return(list("means" = means, "sds" = sds))
#   return (means)
# }

doPHModel <- function(train.task = NULL, test.task = NULL, n.iter = 100) {
  res.aggr = matrix(data = NA, nrow = n.iter, ncol = 2)
  colnames(res.aggr) = c("ci", "lgrk")
  for (i in 1:n.iter) {
    tryCatch({
      print(i)
      set.seed(seed[i])
      lrn <- makeLearner(cl = "surv.coxph")
      rdesc  <-  makeResampleDesc(method = "CV", iters = 10L)
      res = resample(lrn, train.task, rdesc, measures = list(ph.cindex, lgrk),show.info = FALSE)
      temp = res$aggr
      ####test join from ci
      # resp = res$pred
      # a = Hmisc::rcorr.cens(x = getPredictionResponse(resp),
      #                       S = getPredictionTruth(resp))
      # res.list[[i]] = list("ci" = a[["C Index"]], "mean.aggr" = temp)
      res.aggr[i,"ci"] = temp[1]
      res.aggr[i,"lgrk"] = temp[2]
      #hrph.vec[i] = temp[3]
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  means = colMeans(res.aggr,na.rm = TRUE)
  return (means)
}

PHPred <- function(train.task = NULL, test.task = NULL, n.iter = 1) {
  x = list("vector", n.iter)
  #hrph.vec <- vector("double",n.iter)
  for (i in 1:n.iter) {
    set.seed(seed[i])
    lrn <- makeLearner(cl = "surv.coxph")
    rdesc  <-  makeResampleDesc(method = "CV", iters = 10L)
    res = resample(lrn, train.task, rdesc, measures = list(ph.cindex, lgrk),show.info = FALSE)
    resp = res$pred
    x[[i]] = getPredictionResponse(resp)
    print(i)
  }
  return(x)
}