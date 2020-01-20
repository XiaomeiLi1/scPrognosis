getPseudotime <- function(PTime = NULL){
# Scale pseudotime to 0 to 1.
  flag = is.numeric(PTime)
  if(flag) PTime = (PTime-min(PTime))/(max(PTime)-min(PTime))
  else
    stop("Please input numeric time vector!\n\n")
}