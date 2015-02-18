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

testReturn = function(poss.tf) {
  if (sum(poss.tf, na.rm=T) < 10) {
    poss = possComb(poss.tf)
    pr = possToI(poss)
    if (length(unlist(pr)) == length(unique(unlist(pr)))) {
      return(pr)
    }
  }
  return(F)
}

formatPair = function(x, mpc.m2, pr, step=1, ppm) {
  if (nrow(pr) < 1) return(noMatchPair(step))
  cbind(
    c12.g = x[pr[,"col"],"group"], 
    c13.g= x[pr[,"row"],"group"], 
    mpc = mpc.m2[cbind(pr[,"row"],pr[,"col"])], 
    step=1, 
    ppm = ppm[cbind(pr[,"row"],pr[,"col"])], 
    mult = nrow(pr),
    detected.12 = rowSums(!is.na(x[pr[,"col"],c("pn.a", "pn.b")])),
    detected.13 = rowSums(!is.na(x[pr[,"row"],c("pn.a", "pn.b")]))
  )
}

noMatchPair = function(step=NA) {
  cbind(
    c12.g = NA, 
    c13.g= NA, 
    mpc = NA, 
    step=step, 
    ppm = NA, 
    mult = 0,
    detected.12 = 0,
    detected.13 = 0
  )
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
  
  poss.m2 = abs(ppm) < ppm.lim & mpc.tf  
  pr = testReturn(poss.m2)
  if (!identical(pr, F)) return(formatPair(x, mpc.m2, pr, step=1, ppm))
  
  # Find base of seq
  seq.base = unlist(llply(seq_along(x[,1]), function(i) {
    y = x[i,]
    this.seq= which(x[,"seq"] == y[["seq"]])
    this.cn = cn[,1]
    z = x[this.seq,,drop=F]
    
    c = coefficients(lm(maxo ~ mz, data=z))[["mz"]]
    if (is.na(c)) {return(T)}
    if (c < 1) {
      return(any(i == which(min(this.cn[this.seq]) == this.cn)))
    }
    if (c > 1) {
      return(any(i == which(max(this.cn[this.seq]) == this.cn)))
    }
  }))
  seq.base.aperm = aperm(array(seq.base, dim = rep(length(seq.base),2)),c(2,1))
    
  # Find seq direction
  seq.dir = ddply(x, "seq", function(y) {
    cbind(seq = y[1,"seq"], dir = coefficients(lm(maxo ~ mz, data=y))[["mz"]])
  })
  seq.dir.v = merge(x, seq.dir, by="seq")[,c("dir","mz")]
  seq.dir.v = seq.dir.v[order(seq.dir.v[,"mz"]),"dir"]
  seq.dir.comp = outer(seq.dir.v, seq.dir.v, "*")
  seq.dir.tf = seq.dir.comp < 0 | is.na(seq.dir.comp)
  
  
  poss.m2 = abs(ppm) < ppm.lim & mpc.tf & seq.base.aperm & seq.base & seq.dir.tf
  pr = testReturn(poss.m2)
  if (!identical(pr, F)) return(formatPair(x, mpc.m2, pr, step=2, ppm))
  
  if (sum(poss.m2, na.rm=T) > 25) { # Any more and this gets really slow.  Bug but doesn't happen often.
    return(noMatchPair(step=0))
  }
  
  poss = possComb(poss.m2)
  carbon.dist = unlist(alply(poss, 1, function(x) {
    t = laply(which(!is.na(x)), function(i) {
      cn[x[[i]],i]
    })
    prod(t)
  }))
  c.dist.best = which(carbon.dist == max(carbon.dist))
  pr = possToI(poss[c.dist.best,,drop=F])
  
  return(formatPair(x, mpc.m2, pr, step=3, ppm))
}