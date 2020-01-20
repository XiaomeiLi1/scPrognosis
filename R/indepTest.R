library(genefu)

makeDS <- function(dd = NULL, pd = NULL, list = NULL){
  temp = t(dd)
  expd = NULL
  for (i in list) {
    idx = which(colnames(temp) == i)
    if (length(idx)) expd = cbind(expd, temp[,idx])
    else expd = cbind(expd, rep(0, nrow(temp)))
  }
  colnames(expd) = list
  # t1 = as.numeric(pd[,1])
  # t1[is.na(t1)] <- 0
  # t1[t1<0] <- 0
  # t2 = as.numeric(pd[,2])
  # pd = data.frame(cbind(t1,t2))
  colnames(pd) = c("time", "status")
  if(nrow(expd) == 1) {data = cbind(t(expd), pd)
  colnames(data) = c("G1","time", "status")
  }
  else data = cbind(pd,expd)
  datasets <- data.frame(data,stringsAsFactors = F)
  return(datasets)
}

####training on TCGA
indepTest <- function(signatures = NULL) {
  load("/home/xiaomei/Topic1/data/TCGA753.rda")
  #pd = survival_data[,c(5,7)] # rf
  pd = survival_data[,c(6,7)] #os
  train_data = makeDS(RNASeq_data,pd,signatures)
  coxph1 <- coxph(Surv(time,status)~.,data=train_data)
  
  ###OS METABRIC
  pd = pData(METABRIC)[,c(3,4)]
  pd[,2]= as.numeric(pd[,2])
  pd[which(pd[,2] == 2),2] <- 0
  meta_data = makeDS(exprs(METABRIC),pd,signatures)
  predict1 <- predict(coxph1,newdata=meta_data)
  s = Surv(pd[,1], pd[,2])
  a = Hmisc::rcorr.cens(-1 * predict1, s)[["C Index"]]
  
  ####RF
  pd = pData(METABRIC)[,c(5,6)]
  meta_data = makeDS(exprs(METABRIC),pd,signatures)
  predict1 <- predict(coxph1,newdata=meta_data)
  s = Surv(pd[,1], pd[,2])
  a = c(a, Hmisc::rcorr.cens(-1 * predict1, s)[["C Index"]])
  
  ###GEO RF
  pd = pData(GEO)[,c(11,12)]
  pd = apply(pd, 2, as.numeric)
  geo_data = makeDS(exprs(GEO),pd,signatures)
  predict1 <- predict(coxph1,newdata=geo_data)
  s = Surv(pd[,1], pd[,2])
  a = c(a, Hmisc::rcorr.cens(-1 * predict1, s)[["C Index"]])
  
  ###UK DFS
  pd = pData(UK)[,c(20,19)] #rfs
  uk_data = makeDS(exprs(UK),pd,signatures)
  predict1 <- predict(coxph1,newdata=uk_data)
  s = Surv(pd[,1], pd[,2])
  a = c(a, Hmisc::rcorr.cens(-1 * predict1, s)[["C Index"]])
  return(a)
}

####training on METABRIC
indepTest1 <- function(signatures = NULL) {
  ###OS METABRIC
  # pd = pData(METABRIC)[,c(3,4)]
  # pd[,2]= as.numeric(pd[,2])
  # pd[which(pd[,2] == 2),2] <- 0
  # train_data = makeDS(exprs(METABRIC),pd,signatures)
  # coxph1 <- coxph(Surv(time,status)~.,data=train_data)
  
  ####RF
  pd = pData(METABRIC)[,c(5,6)]
  train_data = makeDS(exprs(METABRIC),pd,signatures)
  coxph1 <- coxph(Surv(time,status)~.,data=train_data)
  
  
  ###OS METABRIC
  pd = pData(METABRIC)[,c(3,4)]
  pd[,2]= as.numeric(pd[,2])
  pd[which(pd[,2] == 2),2] <- 0
  meta_data = makeDS(exprs(METABRIC),pd,signatures)
  predict1 <- predict(coxph1,newdata=meta_data)
  s = Surv(pd[,1], pd[,2])
  a = Hmisc::rcorr.cens(-1 * predict1, s)[["C Index"]]
  save(predict1, file = "META_OS_pred.rda")
  
  ####TCGA
  load("C:/unisa/RstudioProjects/Topic1/data/TCGA753.rda")
  pd = survival_data[,c(5,7)] # rf
  tcga_data = makeDS(RNASeq_data,pd,signatures)
  predict1 <- predict(coxph1,newdata=tcga_data)
  s = Surv(pd[,1], pd[,2])
  a = c(a, Hmisc::rcorr.cens(-1 * predict1, s)[["C Index"]])
  save(predict1, file = "TCGA_OS_pred.rda")
  
  pd = survival_data[,c(6,7)] #os
  predict1 <- predict(coxph1,newdata=tcga_data)
  s = Surv(pd[,1], pd[,2])
  a = c(a, Hmisc::rcorr.cens(-1 * predict1, s)[["C Index"]])
  save(predict1, file = "TCGA_RF_pred.rda")
  
  ###GEO RF
  pd = pData(GEO)[,c(11,12)]
  pd = apply(pd, 2, as.numeric)
  geo_data = makeDS(exprs(GEO),pd,signatures)
  predict1 <- predict(coxph1,newdata=geo_data)
  s = Surv(pd[,1], pd[,2])
  a = c(a, Hmisc::rcorr.cens(-1 * predict1, s)[["C Index"]])
  save(predict1, file = "GEO_pred.rda")
  
  ###UK DFS
  pd = pData(UK)[,c(20,19)] #rfs
  uk_data = makeDS(exprs(UK),pd,signatures)
  predict1 <- predict(coxph1,newdata=uk_data)
  s = Surv(pd[,1], pd[,2])
  a = c(a, Hmisc::rcorr.cens(-1 * predict1, s)[["C Index"]])
  save(predict1, file = "UK_pred.rda")
  return(a)
}
