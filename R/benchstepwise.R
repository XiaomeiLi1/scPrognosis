benchstepwise <- function(ranking1 = NULL, ranking2 = NULL, ranking3 = NULL, params = NULL) {
  i = 1
  load("/home/xiaomei/Topic1/data/TCGA753.rda")
  pd = survival_data[,c(5,7)] # os
  a1 = params[i,"a1"]
  a2 = params[i,"a2"]
  a3 = 1 - (a1+a2)
  ranking = geneRank(ranking1, ranking2, ranking3,a1,a2,a3)
  A = NULL
  for(N in 1:100) {
    print(N)
    tryCatch({
      A = rbind(A,run_topN(mRNA = RNASeq_data, pd = pd, ranking = ranking, N = N))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  write.csv(A, file = paste("tcga.vim.os.",a1,a2,".csv",sep = ""))
  
  i = i+1
  pd = survival_data[,c(6,7)] #rf
  a1 = params[i,"a1"]
  a2 = params[i,"a2"]
  a3 = 1 - (a1+a2)
  ranking = geneRank(ranking1, ranking2, ranking3,a1,a2,a3)
  A = NULL
  for(N in 1:100) {
    print(N)
    tryCatch({
      A = rbind(A,run_topN(mRNA = RNASeq_data, pd = pd, ranking = ranking, N = N))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  write.csv(A, file = paste("tcga.vim.rf.",a1,a2,".csv",sep = ""))
  
  ###OS METABRIC
  survival_data = pData(METABRIC)
  pd = survival_data[,c(3,4)]
  pd[,2]= as.numeric(pd[,2])
  pd[which(pd[,2] == 2),2] <- 0
  i = i+1
  a1 = params[i,"a1"]
  a2 = params[i,"a2"]
  a3 = 1 - (a1+a2)
  ranking = geneRank(ranking1, ranking2, ranking3,a1,a2,a3)
  A = NULL
  for(N in 1:100) {
    print(N)
    tryCatch({
      A = rbind(A,run_topN(mRNA = exprs(METABRIC), pd = pd, ranking = ranking, N = N))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  write.csv(A, file = paste("meta.vim.os.",a1,a2,".csv",sep = ""))
  
  ####RF
  pd = survival_data[,c(5,6)]
  i = i+1
  a1 = params[i,"a1"]
  a2 = params[i,"a2"]
  a3 = 1 - (a1+a2)
  ranking = geneRank(ranking1, ranking2, ranking3,a1,a2,a3)
  A = NULL
  for(N in 1:100) {
    print(N)
    tryCatch({
      A = rbind(A,run_topN(mRNA = exprs(METABRIC), pd = pd, ranking = ranking, N = N))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  write.csv(A, file = paste("meta.vim.rf.",a1,a2,".csv",sep = ""))

  ###GEO RF
  survival_data = pData(GEO)
  pd = survival_data[,c(11,12)]
  pd = apply(pd, 2, as.numeric)
  i = i+1
  a1 = params[i,"a1"]
  a2 = params[i,"a2"]
  a3 = 1 - (a1+a2)
  ranking = geneRank(ranking1, ranking2, ranking3,a1,a2,a3)
  A = NULL
  for(N in 1:100) {
    print(N)
    tryCatch({
      A = rbind(A,run_topN(mRNA = exprs(GEO), pd = pd, ranking = ranking, N = N))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  write.csv(A, file = paste("geo.vim.rf.",a1,a2,".csv",sep = ""))

  ###UK DFS
  survival_data = pData(UK)
  pd = survival_data[,c(20,19)] #rfs
  i = i+1
  a1 = params[i,"a1"]
  a2 = params[i,"a2"]
  a3 = 1 - (a1+a2)
  ranking = geneRank(ranking1, ranking2, ranking3,a1,a2,a3)
  A = NULL
  for(N in 1:100) {
    print(N)
    tryCatch({
      A = rbind(A,run_topN(mRNA = exprs(UK), pd = pd, ranking = ranking, N = N))
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  write.csv(A, file = paste("uk.vim.rf.",a1,a2,".csv",sep = ""))
}