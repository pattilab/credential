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
  sc.range = c(which.min(abs(rt.corr - peak[["rtmin"]])), which.min(abs(rt.corr - peak[["rtmax"]])))
  
  do.call("cbind", rawEIC(xr, mzrange=peak[c("mzmin", "mzmax")], scanrange=sc.range))
}
