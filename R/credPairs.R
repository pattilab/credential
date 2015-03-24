fillMat = function(x, byrow=F) {
  matrix(x, nrow=length(x), ncol=length(x), byrow=byrow)
}

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
    c12.pn.a = x[pr[,"col"],"pn.a"], 
    c12.pn.b = x[pr[,"col"],"pn.b"], 
    c13.g= x[pr[,"row"],"group"], 
    c13.pn.a = x[pr[,"row"],"pn.a"], 
    c13.pn.b = x[pr[,"row"],"pn.b"], 
    mpc = mpc.m2[cbind(pr[,"row"],pr[,"col"])], 
    step=step, 
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
    c12.pn.a = NA, 
    c13.pn.a= NA, 
    c12.pn.b = NA, 
    c13.pn.b= NA, 
    mpc = NA, 
    step=step, 
    ppm = NA, 
    mult = 0,
    detected.12 = 0,
    detected.13 = 0
  )
}

credPairs = function(x, ppm.lim, mzdiff, mpc, min.c = 3, log.maxor.r = c(-.9, .9)) {
  x = x[order(x$mz),]
  
  mz.m2 = array(x$mz, dim=c(length(x$mz), length(x$mz)))
  dmz = outer(x$mz, x$mz, "-")
  cn = round(dmz * x[,"charge"] / mzdiff)
  ppm = (cn / x[,"charge"] * mzdiff - dmz) / x$mz * 1E6
  
  #Min carbons
  cn.tf = cn >= min.c
  
  # Maxor
  maxo.m2 = outer(x$maxo, x$maxo, "/")
  maxo.tf = log10(maxo.m2) > log.maxor.r[1] & log10(maxo.m2) < log.maxor.r[2]
  
  # Build mpc tf matrix
  mpc.m2 = mz.m2*x[,"charge"]/cn
  mpc.max = mpc.m2; mpc.min = mpc.m2;
  mpc.max[] = mpc[as.character(cn),"max_mpc"]
  mpc.min[] = mpc[as.character(cn),"min_mpc"]
  mpc.tf = mpc.max > mpc.m2 & mpc.m2 > mpc.min
  
  poss.m2 = abs(ppm) < ppm.lim & mpc.tf & cn.tf & maxo.tf
  pr = testReturn(poss.m2)
  if (!identical(pr, F)) return(formatPair(x, mpc.m2, pr, step=1, ppm))
  
  
  for (step in 2:3) {
    
    #Find edges of the seqs
    seq.edges = unlist(llply(seq_along(x[,1]), function(i) {
      y = x[i,]
      my.cn= cn[i,1]
      this.seq= which(x[,"seq"] == y[["seq"]])
      this.cn = cn[this.seq,1]
      z = x[this.seq,,drop=F]
      #find fuzz
      fuzz = z[,"maxo"]/mean(z[,"maxo"]) < .1
      if (my.cn == min(this.cn[!fuzz])) {
        return(1)
      } else if (my.cn == max(this.cn[!fuzz])) {
        return(-1)
      } else {
        return(0)
      }
    }))
    
    # Find seq direction
    seq.dir = ddply(x, "seq", function(z) {
      fuzz = z[,"maxo"]/mean(z[,"maxo"]) < .1
      c(seq = z[1,"seq"], dir = coefficients(lm(maxo ~ mz, data=z[!fuzz,]))[["mz"]])
    })
    
    #Find is a peak is the base of a seq based on edge/direction
    seq.base.1 = merge(seq.dir, cbind(seq.edges, seq = x[,"seq"], n = seq(seq.edges)), by="seq")
    seq.base.1 = seq.base.1[order(seq.base.1[,"n"]),]
    seq.base.2 = seq.base.1[,2] * seq.base.1[,3]
    seq.base = is.na(seq.base.2) | sign(seq.base.2) == -1
    
    #Find if peaks seqs are opposite facing
    seq.dir.comp = outer(sign(seq.base.1[,"dir"]), sign(seq.base.1[,"dir"]), "*")
    seq.dir.tf = seq.dir.comp == -1 | is.na(seq.dir.comp)
    poss.m2 = abs(ppm) < ppm.lim & mpc.tf & cn.tf & maxo.tf &
      (fillMat(seq.edges) != 0 |  fillMat(seq.edges, T) != 0) &
      outer(seq.base, seq.base, "&") &
      seq.dir.tf
    pr = testReturn(poss.m2)
    if (!identical(pr, F)) if(nrow(pr) > 0 | step == 3) return(formatPair(x, mpc.m2, pr, step=step, ppm))
    
    if (!identical(pr, F)) if(nrow(pr) == 0) {
      x = ddply(x, "seq", function(z) {
        fuzz = z[,"maxo"]/mean(z[,"maxo"]) < .1
        mzs = round(z$mz, 2)
        mzs.u = unique(mzs)
        mid = round(sum(!fuzz)/2)
        z[mzs %in% mzs.u[(mid+1):length(mzs)],"seq"] = 100*z[1,"seq"]
        z
      })
      x = x[order(x$mz),]
    }
  }
      
  
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
  pr = unique(pr)
  
  return(formatPair(x, mpc.m2, pr, step=4, ppm))
}