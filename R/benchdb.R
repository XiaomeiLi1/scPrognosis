benchdb <- function(ranking = NULL, N = 50) {
  pd = survival_data[,c(5,7)] # os
  A = run_topN(mRNA = RNASeq_data, pd = pd, ranking = ranking,N = N)
  
  pd = survival_data[,c(6,7)] #rf
  A = rbind(A,run_topN(mRNA = RNASeq_data, pd = pd, ranking = ranking,N = N))
  
  ###OS METABRIC
  survival_data = pData(METABRIC)
  pd = survival_data[,c(3,4)]
  pd[,2]= as.numeric(pd[,2])
  pd[which(pd[,2] == 2),2] <- 0
  A = rbind(A,run_topN(mRNA = exprs(METABRIC), pd = pd, ranking = ranking,N = N))
  
  ####RF
  pd = survival_data[,c(5,6)]
  A = rbind(A,run_topN(mRNA = exprs(METABRIC), pd = pd, ranking = ranking,N = N))
  
  ###GEO RF
  survival_data = pData(GEO)
  pd = survival_data[,c(11,12)]
  pd = apply(pd, 2, as.numeric)
  A = rbind(A,run_topN(mRNA = exprs(GEO), pd = pd, ranking = ranking,N = N))
  
  ###UK DFS
  survival_data = pData(UK)
  pd = survival_data[,c(20,19)] #rfs
  A = rbind(A,run_topN(mRNA = exprs(UK), pd = pd, ranking = ranking,N = N))
  return(A)
}
