formatPeaks = function(
  mf, 
  ps
) {
  c("mz", "rt", "rtmax", "rtmin", "match_mz", "n_carbons", "peaknum_a", "adduct", "iso", "charge")
  
  if (nrow(mf) < 1) {return("")}
  
  tmp = lapply(1:nrow(mf), function(i) {
    peak = mf[i,,drop=F]
    peak_a = ps[which(ps[,"peaknum"] == peak[,"master_peaknum_a"]),,drop=F]
    match_a = ps[which(ps[,"peaknum"] == peak[,"peaknum_a"]),,drop=F]

    cbind(
      peaknum_a = peak_a[,"peaknum"],
      peaknum_a13 = match_a[,"peaknum"],
      mz= peak_a[,"mz"], 
      rt= peak_a[,"rt"], 
      rtmin= peak_a[,"rtmin"], 
      rtmax= peak_a[,"rtmax"],
      match_mz = match_a[,"mz"],
      n_carbons = peak[,"p_carbons_a"],
      adduct = peak[,"carule"],
      iso = peak[,"ciso"],
      charge = peak[,"p_charge_a"],
      maxo = peak_a[,"maxo"]
    )
  })
  do.call("rbind", tmp)
}
