f1 <- function(n){
  i1 <- sample.int(n, 1)
  i2 <- sample.int(n-1, 1)
  i3 <- sample.int(n-2, 1)
  elems <- 0:(n-1)
  elems <- elems[-i1]
  j2 <- elems[i2]
  elems <- elems[-i2]
  j3 <- elems[i3]
  return(c(i1, j2, j3))
}

f2 <- function(n){
  i1 <- sample.int(n, 1)
  i2 <- sample.int(n-1, 1)
  i3 <- sample.int(n-2, 1)
  elems <- 1:n
  elems[i1] = elems[n]
  j2 <- elems[i2]
  elems[i2] = elems[n-1]
  j3 <- elems[i3]
  return(c(i1, j2, j3)-1)
}

f3 <- function(n){
  i1 <- sample.int(n, 1) - 1
  i2 <- sample.int(n-1, 1) - 1
  i3 <- sample.int(n-2, 1) - 1
  if(i3 == i2) i3 = n-2
  if(i3 == i1) i3 = n-1
  if(i2 == i1) i2 = n-1
  return(c(i1, i2, i3))
}

table(
  replicate(480000, paste0(f1(4), collapse = "-"))
)

table(
  replicate(120000, paste0(f2(4), collapse = "-"))
)

table(
  replicate(480000, paste0(f3(5), collapse = "-"))
)