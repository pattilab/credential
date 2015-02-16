mov.av = function(x, n=3) {
  temp = filter(x,rep(1/n,n), sides=2)
  temp[is.na(temp)] = x[is.na(temp)]
  as.numeric(temp)
}

nmEIC.an = function(pn, an, xr.l) {
  peak = an@xcmsSet@peaks[pn,]
  xr = xr.l[[peak["sample"]]]
  rt.corr =an@xcmsSet@rt$corrected[[peak["sample"]]]
  
  nmEIC.rt(peak, xr, rt.corr)
}

nmEIC.rt = function(peak, xr, rt.corr) {# Peak has rtmin, rtmax, mzmin, mzmax
  sc.range = c(order(abs(rt.corr-min(peak["rtmin"])), decreasing=F)[1],
               order(abs(rt.corr-max(peak["rtmax"])), decreasing=F)[1])
  
  eic = rawEIC(xr, mzrange=peak[c("mzmin", "mzmax")], scanrange=sc.range)$intensity
  eic.r = cbind(scan = sc.range[1]:sc.range[length(sc.range)], int = eic, pn = peak[["pn"]])
  data.frame(eic.r)
}