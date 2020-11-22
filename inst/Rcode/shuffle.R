
f_sample_k <- function(k, nsims){
  replicate(nsims, paste0(sample.int(k), collapse = "-"))
}

f_sample_n_k <- function(n, k, nsims){
  replicate(nsims, paste0(sample.int(n, k), collapse = "-"))
}

f_shuffle <- function(n, k, nsims){
  replicate(nsims, paste0(sample.int(n)[1:k], collapse = "-"))
}

table(f_sample_n_k(5, 3, 50000))/50000
table(f_shuffle(5, 3, 50000))/50000