convolve = function(...) {
  l = list(...)
  
  a = l[[1]]
  for (i in seq_len(length(l))[-1]) {
    b = l[[i]]
    
    prod = outer(a, b, "*")
    inds = row(prod) + col(prod)
    
    convolved = sapply(unique(c(inds)), function(x) {
      sum(prod[inds == x])
    })
    a = convolved/max(convolved)
  }
  a
}


cbind.fill <- function(...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}


clustgroup = function(mat, scales, factor) {
  if (nrow(mat) < 2) return(rep(1, nrow(mat)))
  
  cutree(hclust(dist(mat/rep(scales, each = nrow(mat))), method="single"), h = factor)
}

breakgroup = function(x, breaksize = 1) {
  o = order(x)
  ass = c(0, which(diff(x[o]) > breaksize), length(x))
  
  as.integer(cut(seq_along(x), breaks = ass))[order(o)]
}