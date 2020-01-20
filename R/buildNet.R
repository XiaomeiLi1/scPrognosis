library(LEAP)

# data: data matrix for which the rows are genes and the columns are experiments, sorted by their pseudotime.
netLeap <- function(data = NULL, cutoff = 0.2, info = FALSE) {
  gn = rownames(data)
  start.time <- Sys.time()
  MAC = MAC_counter(data = data, MAC_cutoff = cutoff)
  if(info) save(MAC, file = "MAC.rda")
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  #cat("Running time: ", time.taken, "\n\n")
  res = MAC[,c(4,3,1)]
  res[,1] = gn[as.numeric(res[,1])]
  res[,2] = gn[as.numeric(res[,2])]
  colnames(res) = c("child","parent","weight")
  res
}