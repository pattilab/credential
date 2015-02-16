possComb = function(poss) {
  row.poss = alply(poss, 2, function(x){
    if (!any(x, na.rm=T)) { return(NA) }
    c(which(x), NA)
  })
  poss.grid = expand.grid(row.poss)
  poss.grid.m = as.matrix(poss.grid)
  poss.grid.rep = adply(poss.grid.m, 1, testRow)
  poss.grid.rep[,-1]
}

testRow = function(x) {
  x = matrix(x, nrow=1)
  if(any(count(na.omit(unlist(x)))[,"freq"] > 1)) return(NULL)
  if(any(count(na.omit(unlist(possToI(x))))[,"freq"] > 1)) return(NULL)
  return(x)
}

possToI = function(poss) { # Helper function to convert rows of possibilities to row, col
  pairs = which(!is.na(poss), arr.ind=T)
  data.frame(row = poss[pairs], col = pairs[,"col"])
}

credPairs = function(x, ppm.lim, mzdiff, mpc) {
  x = x[order(x$mz),]
  
  mz.m2 = array(x$mz, dim=c(length(x$mz), length(x$mz)))
  dmz = outer(x$mz, x$mz, "-")
  cn = round(dmz * x[,"charge"] / mzdiff)
  ppm = (cn / x[,"charge"] * mzdiff - dmz) / x$mz * 1E6
  
  # Build mpc tf matrix
  mpc.m2 = mz.m2*x[,"charge"]/cn
  mpc.max = mpc.m2; mpc.min = mpc.m2;
  mpc.max[] = mpc[as.character(cn),"max_mpc"]
  mpc.min[] = mpc[as.character(cn),"min_mpc"]
  mpc.tf = mpc.max > mpc.m2 & mpc.m2 > mpc.min
  
  #maxo.r.a = outer(x$maxo.a, x$maxo.a, "/")
  #maxo.r.b = outer(x$maxo.b, x$maxo.b, "/")
  
  poss.m2 = abs(ppm) < ppm.lim & mpc.tf
  if (sum(poss.m2, na.rm=T) < 3) {
    poss = possComb(poss.m2)
    pr = possToI(poss)
    pair = cbind(
      c12.g = x[pr$col,"group"], 
      c13.g= x[pr$row,"group"], 
      mpc = mpc.m2[cbind(pr$row,pr$col)], 
      step=1, 
      ppm = ppm[cbind(pr$row,pr$col)], 
      mult = nrow(pr),
      detected.12 = sum(is.na(x[pr$col,c("pn.a", "pn.b")])),
      detected.13 = sum(is.na(x[pr$row,c("pn.a", "pn.b")]))
    )
    return(pair) 
  }
  
  # Find seq direction
  seq.dir = ddply(x, "seq", function(y) {
    cbind(seq = y[1,"seq"], dir = coefficients(lm(maxo ~ mz, data=y))[["mz"]])
  })
  seq.dir.v = merge(x, seq.dir, by="seq")[,c("dir","mz")]
  seq.dir.v = seq.dir.v[order(seq.dir.v[,"mz"]),"dir"]
  seq.dir.comp = outer(seq.dir.v, seq.dir.v, "*")
  seq.dir.tf = seq.dir.comp < 0 | is.na(seq.dir.comp)
  
  
  # Find base of seq
  seq.base = unlist(alply(x, 1, function(y) {
    z = subset(x, seq == y[["seq"]])
    c = coefficients(lm(maxo ~ mz, data=z))[["mz"]]
    if (is.na(c)) {return(T)}
    if (c < 1) {return(y["pn.a"] == z[1,"pn.a"])}
    if (c > 1) {return(z[nrow(z), "pn.a"] == y["pn.a"])}
  }))
  seq.base.aperm = aperm(array(seq.base, dim = rep(length(seq.base),2)),c(2,1))
  
  poss.m2 = abs(ppm) < ppm.lim & mpc.tf & seq.base.aperm & seq.base & seq.dir.tf
  if (sum(poss.m2, na.rm=T) < 3) {
    poss = possComb(poss.m2)
    pr = possToI(poss)
    pair = cbind(
      c12.g = x[pr$col,"group"], 
      c13.g= x[pr$row,"group"], 
      mpc = mpc.m2[cbind(pr$row,pr$col)], 
      step=2, 
      ppm = ppm[cbind(pr$row,pr$col)], 
      mult = nrow(pr),
      detected.12 = sum(is.na(x[pr$col,c("pn.a", "pn.b")])),
      detected.13 = sum(is.na(x[pr$row,c("pn.a", "pn.b")]))
    )
    return(pair) 
  }
  
  if (sum(poss.m2, na.rm=T) > 25) { 
    poss.m2[] = NA; poss = possComb(poss.m2)
  } else { # Any more and this gets really slow.  Bug.
    poss = possComb(poss.m2)
  }
  carbon.dist = unlist(alply(poss, 1, function(x) {
    t = laply(which(!is.na(x)), function(i) {
      cn[x[[i]],i]
    })
    prod(t)
  }))
  c.dist.best = which(carbon.dist == max(carbon.dist))
  
  pr = possToI(poss[c.dist.best,,drop=F])
  pair = cbind(
    c12.g = x[pr$col,"group"], 
    c13.g= x[pr$row,"group"], 
    mpc = mpc.m2[cbind(pr$row,pr$col)], 
    step=3, 
    ppm = ppm[cbind(pr$row,pr$col)], 
    mult = nrow(pr),
    detected.12 = sum(!is.na(x[pr$col,c("pn.a", "pn.b")])),
    detected.13 = sum(!is.na(x[pr$row,c("pn.a", "pn.b")]))
  )
  return(pair)
}