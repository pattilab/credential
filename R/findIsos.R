buildIsoGroups = function(
  peaks,
  xr,
  rt.corr,
  ppm.lim,
  rcorr.lim,
  rt.lim,
  mzdiff,
  charges,
  g=100
) { cat("Building isotope groups within a sample. Total: ", nrow(peaks), "\n"); #Misses peaks when they arent within rt.lim of the proper group.  Possible considering centwave's rt method.
    carryover.pns = c()
    last.pn = 0
    isogs = list()
    
    repeat {
      end = last.pn+g
      if (end > nrow(peaks)) {end = nrow(peaks)}
      if ((last.pn+1) >= end) {break;}
      
      p.now = subset(peaks, pn %in% c(peaks[(last.pn+1):end,"pn"], carryover.pns)) 
      cat(paste0(last.pn, "+", nrow(p.now), ";"))
      if(nrow(p.now) < 1) { next;}
      carryover.pns = p.now[which(abs(p.now[,"rt"] - max(p.now$rt)) < rt.lim), "pn"]
      
      isogs.now = findIsos(
        peaks = p.now,
        xr = xr,
        rt.corr = rt.corr,
        ppm.lim = ppm.lim,
        rcorr.lim = rcorr.lim,
        mzdiff = mzdiff,
        charges = charges
      )
      
      co.g = laply(isogs.now, function(x) {
        any(carryover.pns %in% x[,"pn"])
      })
      more.co = unlist(llply(isogs.now[which(co.g)], function(x) x[,"pn"]))
      carryover.pns = unique(c(carryover.pns, more.co))
      
      isogs = c(
        isogs, 
        isogs.now[-which(co.g)]
      )
      
      last.pn = end;
    }
    isogs
}



findIsos = function(
  peaks,
  xr,
  rt.corr,
  ppm.lim,
  rcorr.lim,
  mzdiff = atm$c13 - atm$c12,
  charges = c(1,2,3,4,5,6)
) { # Takes a data.frame of peaks and returns a list, each entry is a different isotope channel with members assigned to sequences.
  
  # Two m3 peak matrices, the first has peaks$mz along the z axis, the second has peaks$mz along the y axis
  peaks.m3 = outer(outer(peaks$mz, rep(1, length(peaks$mz)), "*"), rep(1, length(charges)), "*")
  peaks2.m3 = outer(outer(rep(1, length(peaks$mz)), peaks$mz, "*"), rep(1, length(charges)), "*")
  
  #Pairwise differences and carbon number calculation
  mzd.m2 = outer(peaks$mz, peaks$mz, "-")
  mzd.m2[lower.tri(mzd.m2, F)] = NA
  cn.m3 = outer(mzd.m2, charges, function(mzd, chg) {
    round(mzd * chg / mzdiff)
  })
  pmz.m3 = outer(mzd.m2, charges, function(mzd, chg) {
    mzdiff/chg
  }) * cn.m3 + peaks2.m3
  
  ppm.m3 = (pmz.m3 - peaks.m3) / pmz.m3 * 1E6
  
  # Which pairs exist within the mass error tolerance
  isos = which(abs(ppm.m3) < ppm.lim, arr.ind=T)
  dimnames(isos) = list(NULL, c("base", "iso", "charge"))

  # Finding eics of putative matches only once. names(eics) correspond to rownumber of peak or index in peaks$mz
  needed.eics = unique(as.vector(isos[,c("base","iso")]))
  eics.l = llply(needed.eics, function(x) nmEIC.rt(data.frame(peaks[x,,drop=F]), xr, rt.corr))
  
  scans = unlist(llply(eics.l, function(x) x[,"scan"]))
  scan.min = min(scans)
  scan.max = max(scans)
  eic.mat = matrix(NA, ncol=length(needed.eics), nrow=scan.max-scan.min+1, dimnames=list(NULL, as.character(needed.eics)))
  
  for(i in seq_along(eics.l)) {
    eic = eics.l[[i]]
    eic.mat[eic$scan-scan.min+1, as.character(i)] = eic$int
  }
  
  #Finding corr between putative matches
  eic.corr = ddply(data.frame(isos), c("base", "iso"), function(x) {
    eic.merge = eic.mat[,as.character(x[1,c("base", "iso")])]
    if (nrow(eic.merge) < 5) {return(data.frame(rcorr = 0))}
    data.frame(rcorr = rcorr(eic.merge[,1], eic.merge[,2])$r[1,2])
  })
  
  #Convert this 3 column df into an array conformable with ppm.m3
  missing.base = which(!(1:length(peaks$mz) %in% eic.corr$base))
  missing.iso = which(!(1:length(peaks$mz) %in% eic.corr$iso))
  if (length(missing.base) > 0 & length(missing.iso) > 0) eic.corr = rbind(eic.corr, cbind(base = missing.base, iso=missing.iso, rcorr=NA))
  rcorr.m2 = acast(eic.corr, base~iso, value.var="rcorr", fill=NA)
  rcorr.m3 = outer(rcorr.m2, rep(T, length(charges), "&"))
  
  #Determine charge state and iso packets
  # Remove cn of all disqualified peaks
  not.isos = which(abs(ppm.m3) > ppm.lim | rcorr.m3 < rcorr.lim)
  cn.isos = cn.m3
  cn.isos[not.isos] = NA
  
  #Find 1C spacing series
  cn.spacings = aaply(cn.isos, c(2,3), function(x) {
    c(NA,diff(x[order(x)]))
  })
  sequential.n = aaply(cn.spacings, c(1,2), function(x) {
    sum(x==1,na.rm=T)
  })
  
  #Pick best charge state or 1 temporarily to define groupings
  #possible.isos = !is.na(cn.isos)
  #has.iso = aaply(possible.isos, c(1), any)
  best.charge = aaply(sequential.n, 1, function(x) { order(x, decreasing=T)[1] })
  
  ppm.m3.b.m2 = aaply(1:dim(ppm.m3)[2], 1, function(i) {
    ppm.m3[i,,best.charge[[i]]]
  })
  rcorr.m3.b.m2 = aaply(1:dim(rcorr.m3)[2], 1, function(i) {
    rcorr.m3[i,,best.charge[[i]]]
  })
  
  # partition isotope groups
  isos2 = which(abs(ppm.m3.b.m2) < ppm.lim & rcorr.m3.b.m2 > rcorr.lim, arr.ind=T)
  dimnames(isos2) = list(NULL, c("base", "iso"))
  isos2 = subset(data.frame(isos2), base != iso)
  
  toag = unique(unlist(isos2))
  iso.groups = list()
  for (ag in toag) {
    if(ag %in% unlist(iso.groups)) { next; }
    repeat {
      y = subset(data.frame(isos2), base %in% ag | iso %in% ag)
      ag.n = unique(unlist(y))
      
      if (all(ag.n %in% ag)) { iso.groups[[length(iso.groups)+1]] = ag; break;}
      ag = ag.n
    }
  }
  
  # Choose best charge state based on other groups.
  best.charge.g = laply(iso.groups, function(x) {
    order(aaply(sequential.n[x,], 2, sum), decreasing=T)[1]
  })
  
  # With charge state specified and isotope group assigned we can find sequences.
  cn.g = llply(seq_along(iso.groups), function(i) {
    x = iso.groups[[i]]
    mzs = peaks$mz[x]
    base = peaks$mz[x][1]
    
    cbind(cn=round((base-mzs) * best.charge[[i]] / mzdiff), pn=peaks$pn[x])
  })
  cn.spacings.g = llply(cn.g, function(x) {
    ord = order(x[,"cn"])
    cns = x[ord,"cn"]
    cbind(cn = c(diff(cns), NA), p = x[ord,"pn"])
  })
  iso.seqs = llply(cn.spacings.g, function(x) {
    i = 0
    j = 0
    y = cbind(x, seq=0)
    for (z in c(which(x[,"cn"]>1), (nrow(x)))) {
      j = j+1
      y[(i+1):(z),"seq"] = j
      i = z
      if(i == nrow(x)) {break;}
    }
    y
  })
  
  isos.end= llply(seq_along(iso.groups), function(i) {
    x = iso.groups[[i]]
    i.s = iso.seqs[[i]]
    ord = match(peaks[x,"pn"], i.s[,"p"])
    
    cbind(pn = peaks[x,"pn"], charge = best.charge.g[[i]], seq = i.s[ord,"seq"])
  })
  
  isos.end
}