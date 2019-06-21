get_ticks <- function(data, col_name) {
  n1 <- floor(log10(range(data[col_name])))
  if (n1[1] == -Inf) n1[1] = 0
  pow <- seq(n1[1], n1[2]+1)
  ticks <- as.vector(sapply(pow, function(p) (c(1,5)*10^p)))
  return(ticks)
}
