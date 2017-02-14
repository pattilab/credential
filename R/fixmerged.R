checkmerged = function(knot) {
  if (sum((!knot$tail))<4) { return(rep(1, nrow(knot))) }
  
  cc.o = knot$cc
  o = order(knot$mz) 
  
  half = floor(sum(!knot$tail)/2)
  full = sum(!knot$tail)
  
  front = knot[o][tail == F][1:half]
  back = knot[o][tail == F][((half+1):full)]
  
  s1 = mean(diff(head(front$i,3)))/max(head(front$i,3))
  s2 = mean(diff(tail(back$i,3)))/max(tail(back$i,3))
  
  if (s1 < -0.2 & s2 > 0.2) {
    return(c(rep(1, half), rep(2, nrow(knot)))[match(knot$cc, cc.o[o])])
  }
  else return(rep(1, nrow(knot)))
}

calcdir = function(ps) {
  ps = subset(ps, tail == F)
  
  f = tail(seq_len(nrow(ps)), n=3)
  b = head(seq_len(nrow(ps)), n=3)
  
  ds = c(
    mean(diff(ps$i[order(ps$mz)][f]))/max(ps$i),
    mean(diff(ps$i[order(ps$mz)][b]))/max(ps$i)
  )
  
  val = ds[which.max(abs(ds))]
  if (length(val) < 1) return(as.numeric(NA))
  val
}

knotstats = function(knot) {
  
  o = order(knot$mz)
  
  data.frame(
    meanr = mean(knot$mz %% (13.00335-12)),
    meanmz = mean(knot$mz),
    mainmz = knot$mz[which.max(knot$i)],
    rt = mean(knot$rt),
    maxi = max(knot$i),
    n = nrow(knot),
    z = knot$z[1],
    dir = calcdir(knot[,.(mz, i, tail)])
  )
}


fixmergedquipu = function(Knot, features) {
  
  merged = Knot$cc_knot[features,,on="cc", nomatch=0][,merged := checkmerged(.SD), by="knot"][,.(cc,merged)]
  Knot$cc_knot[merged,merged := merged,on="cc"]
  
  temp = Knot$knot[Knot$cc_knot,,on="knot"]
  
  temp[,':='(knot=as.integer(factor(paste(knot, merged))), merged = NULL)]
  
  knot_new = features[,.(mz, i, rt, cc)][temp,,on="cc", nomatch = 0][,knotstats(.SD),by="knot"]
  
  Knot$knot = knot_new
  
  Knot$cc_knot = temp[,.(cc, knot, tail)]
  
  Knot
}