grid.search <- function(ranking1 = NULL, ranking2 = NULL, ranking3 = NULL, N = 50){
  res = NULL
  for(a1 in seq(0,1,0.1)){
    for(a2 in seq(0,1-a1,0.1)){
      a3 = 1 - (a1+a2)
      ranking = geneRank(ranking1, ranking2, ranking3,a1,a2,a3)
      temp = benchdb(ranking, N)
      if(is.null(res)) res = temp
      else res = cbind(res, temp)
    }
  }
  res
}
