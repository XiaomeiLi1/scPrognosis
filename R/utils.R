run_topN <- function(mRNA = NULL, pd = NULL, ranking = NULL, N = 50){
  ranking = names(ranking) ###notice this line
  index = which(ranking %in% rownames(mRNA))
  t_rank = ranking[index]
  ddata = mRNA[t_rank,]
  genes = t_rank[1:N]
  dataSets = makeCVDataSets(ddata,pd,genes)
  res = doPHModel(dataSets)
  res
}