mov_av <- function(x, n=3){
  temp = filter(x,rep(1/n,n), sides=2)
  nas = is.na(temp)
  temp[which(nas)] = x[which(nas)]
  as.numeric(temp)
}

mPeaks = function(
  peak, 
  peaks, 
  ppm_for_isotopes,
  mpc,
  mpc_f
) {
  p_rt_diff = abs(peaks[,"rt"] - peak["rt"])
  p_maxo_r = peak["maxo"]/peaks[,"maxo"]
  peaks = cbind(peaks, p_maxo_r = p_maxo_r)
  
  tmp = lapply(c(1,2,3), function(charge) { 
    
    p_carbons = round((peaks[,"mz"] - peak["mz"]) * charge / (aC13-aC12))
    p_mass = p_carbons*(aC13-aC12)/charge + peak["mz"]
    p_mpc = peak["mz"] * charge / p_carbons
    p_iso_ppm = abs((peaks[,"mz"] - p_mass)/p_mass) * 1E6
    
    peaks = cbind(peaks, p_charge=rep(charge, nrow(peaks)), p_carbons=p_carbons, p_mass=p_mass, p_mpc = p_mpc, p_iso_ppm = p_iso_ppm)
    
    cs = abs(peaks[,"p_carbons"])
    cs[cs==0] = 1
    maxmpc = mpc[cs, "max_mpc"]*2
    minmpc = mpc[cs, "min_mpc"]
    maxmpc[is.na(maxmpc)] = mpc[nrow(mpc),"max_mpc"]
    minmpc[is.na(minmpc)] = mpc[nrow(mpc),"min_mpc"]
    
    ppm_mz = ppm_for_isotopes(peaks[,"mz"])
    
    x = complete.cases(peaks[,"mz"]) & #No filtering - all match
      #peaks$p_rt_diff < isotope_rt_delta_s & #Exclude all peaks which do not elute within x seconds
      #peaks$mz > peak$mz & #Implied
      peaks[,"p_iso_ppm"] < ppm_mz & #Prospective isotope must be n*(C13-C12) m/z away
      peaks[,"p_mpc"] <= maxmpc * mpc_f &
      peaks[,"p_mpc"] >= minmpc / mpc_f &
      peaks[,"p_carbons"] > 0
      #peaks[,"p_maxo_r"] < max_maxo_r  & #Single Sample maxo Ratio
      #peaks[,"p_maxo_r"] > min_maxo_r  #Single Sample maxo Ratio
    
    peaks[x, c("peaknum", "p_charge", "p_carbons", "p_mpc", "p_iso_ppm")]
  })
  
  do.call("rbind", tmp)
} # Returns data.frame("peaknum", "p_charge", "p_carbons", "p_mpc", "p_rt_diff", "p_iso_ppm") of matches


pwms = function(
  peak_table,
  isotope_rt_delta_s, 
  ppm_for_isotopes,
  mpc,
  mpc_f
) {
  cat("\nLooking for pairwise 13C isotopic peak matches within", nrow(peak_table), "peaks:\n")
  pairwise_matches = alply(peak_table, 1, .progress="text", function(peak) {
    peaks = peak_table[abs(peak_table[,"rt"] - peak["rt"]) < isotope_rt_delta_s, ,drop=F]
    
    mPeaks(peak, peaks, ppm_for_isotopes, mpc, mpc_f)
})
  
  count = sum(sapply(pairwise_matches, function(x) { nrow(x) }), na.rm=T)
  names(pairwise_matches) = peak_table[,"peaknum"]
  cat("\nPeaks:",nrow(peak_table), "Pairwise matches: ",count)
  
  pairwise_matches
}

filterCorrs = function(pwmsx, peak_table, an, xr, sample) {
  cat("\nCalculating EIC correlations between putative isotopes:\n")
  rt.raw = an@xcmsSet@rt$raw[[sample]]
  rt.corr = an@xcmsSet@rt$corrected[[sample]]
  
  tmp = llply(names(pwmsx), .progress="text", function(i) {
    pns = unique(pwmsx[[i]][,"peaknum"])
    pn = as.numeric(i)
    if(length(pns) < 1) {return(pwmsx[[i]])}
    
    peak = peak_table[peak_table[,"peaknum"] %in% pn,,drop=F]
    peaks = peak_table[peak_table[,"peaknum"] %in% pns,,drop=F]
    
    min = min(c(peak[,"rtmin"], peaks[,"rtmin"]))
    max = max(c(peak[,"rtmax"], peaks[,"rtmax"]))
    
    min = order(abs(min - rt.corr), decreasing = F)[1]
    max = order(abs(max - rt.raw), decreasing = F)[1]
    
    min  = rt.raw[min]
    max  = rt.raw[max]
    
    eic_peak = rawEIC(xr, mzrange=peak[,c("mzmin", "mzmax")], rtrange=c(min, max))$intensity
    if(length(eic_peak) > 6) {eic_peak = mov_av(eic_peak)}
    
    eic_peaks = alply(peaks, 1, function(x) {
      temp = rawEIC(xr, mzrange=x[c("mzmin", "mzmax")], rtrange=c(min, max))$intensity
      if(length(temp) > 6) {temp = mov_av(temp)}
      temp
    })
    
    corrs_peaks = laply(eic_peaks, function(x) {
      if(length(eic_peak) < 5 || length(x) < 5) {return(1)}
      rcorr(eic_peak, x)[[1]][1,2]
    })
    
    pns_passed = pns[corrs_peaks > 0.7]
    
    pwmsx[[i]][
      pwmsx[[i]][,"peaknum"] %in% pns_passed,
      ,drop=F
      ]
  })
  names(tmp) = names(pwmsx)
  tmp
}

filterMpc = function(matches, mpc_f) {
  cat("\nFiltering matches based on mass per carbon limits.\n")
  
  minmpc = mpc[matches[, "p_carbons_a"],"min_mpc"]/mpc_f
  maxmpc = mpc[matches[, "p_carbons_a"],"max_mpc"]*mpc_f

  #foo = matches[
  #  (matches[,"p_mpc_a"] < maxmpc & matches[,"p_mpc_a"] > minmpc)
  #  ,,drop=F]
  idx = which(matches[,"p_mpc_a"] < maxmpc & matches[,"p_mpc_a"] > minmpc)
  foo = matches[idx,]
}

#Looks for credentialed features which are a series of 1 carbon spacings.  This indicates natural abundance isotopes and the more isotopes found the more likely we have identified the correct charge state.
buildChargeSupport = function(
  pwmsx
) { 
  cat("\nTrying to determine charge state based on isotopic spacings on", length(pwmsx), "peaks:\n")
  pairwise_matches = llply(1:length(pwmsx), .progress="text", function(i) {
   ms = pwmsx[[i]]                               
   ms = ms[order(ms[,"p_carbons"]),,drop=F][order(ms[,"p_charge"]),,drop=F]
   
   srch = lapply(c(1,2,3), function(c) {
     peaks = ms[which(ms[,"p_charge"] == c),,drop=F]
     
     spacing = c(peaks[,"p_carbons"][-1] - peaks[,"p_carbons"][-length(peaks[,"p_carbons"])], 0)
     
     for (i in 1:length(spacing)) {
       if(spacing[i] > 1) { spacing[i] = 0 }
     }
     count = 0
     supporting_peaks = c()
     for (i in 1:length(spacing)) {
       supporting_peaks = c(supporting_peaks, count)
       if(spacing[i] == 0) { count = 0 }
       if(spacing[i] == 1) { count = count + 1}
     }
     
     if(nrow(peaks) == 0) { spacing = numeric(); supporting_peaks = numeric() }
     if(nrow(peaks) == 1) { spacing = 0 }
     
     cbind(peaks, charge_support = spacing, total_charge_support = rep(length(which(spacing==1)), nrow(peaks)), supporting_peaks = supporting_peaks)
   })
   
   do.call("rbind", srch)
})

  names(pairwise_matches) = names(pwmsx)
  pairwise_matches
}


buildAlignIndex = function(
  xs
  ) { #TODO: Fix up grouping
  npeaks = sum(xs@peaks[,"sample"]==1)
  hum = ldply(xs@groupidx, .progress="text", function(x) {
    peaks_from_1 = x <= npeaks
    
    if(!any(peaks_from_1) || !any(!peaks_from_1)) { return(NULL) }
    
    ldply(which(peaks_from_1), function(y) {
      cbind(peaknum_a = x[y], peaknum_b = x[!peaks_from_1] - npeaks)
    })
  })
  
  cat("\nUnique peaknum_a aligns: ", length(unique(hum[,"peaknum_a"])), ". Total aligns: ", nrow(hum))
  return(hum)
}

matchAlignsCarbonCharge = function(
  aligns, 
  pwma, 
  pwmb
  ) {
  cat("\nBefore combination and carbon/charge filtering. Unique peaknum_a aligns: ", length(unique(aligns[,"peaknum_a"])), ". Total aligns: ", nrow(aligns))
  cat("\nWorking on peak of", nrow(aligns), "total:\n")
  all_matches = llply(1:nrow(aligns), .progress="text", function(i) { #Looks for matches in a, that have matches in b with the same carbon number
    x = aligns[i,,drop=F]
    
    ma = pwma[[as.character(x[,"peaknum_a"])]]
    mb = pwmb[[as.character(x[,"peaknum_b"])]]
        
    tmp = lapply(seq_along(ma[,1]), function(j) {
      y = ma[j,,drop=F]
      
      possible_validation = 
        mb[,"p_carbons"] == y[,"p_carbons"] &
        mb[,"p_charge"] == y[,"p_charge"]
      
      if (sum(possible_validation) > 0) {
        ys = y
        #for (z in 1:(sum(possible_validation)-1)) { ys = rbind(y, ys) }
        
        rows = nrow(mb[possible_validation,,drop=F])
        
        possibilities = cbind(
          rep(as.numeric(x[,"peaknum_a"]), rows), 
          matrix(rep(ys, rows), nrow = rows, byrow=T), 
          rep(as.numeric(x[,"peaknum_b"]), rows),
          mb[possible_validation,,drop=F]
          )
        colnames(possibilities) = c("master_peaknum_a", paste(sep="", colnames(ys), "_a"), "master_peaknum_b", paste(sep="", colnames(mb), "_b"))
        
        return(possibilities)
      }
    })
    do.call("rbind",tmp,quote=T)    
  })
  hum = do.call("rbind", all_matches, quote=T)
  cat("\nAfter combination and carbon/charge filtering between samples. Unique peaknum_a aligns: ", length(unique(hum[,"master_peaknum_a"])), ". Total possibilities: ", nrow(hum))
  
  hum
  # list(A12 ->
  #   list(A13 ->
  #       list(B12 ->
}


pickLowestCharge = function(
  ams
  ) {
  cat("\nBefore selecting lowest charge. Unique: ", length(unique(ams[,"master_peaknum_a"])), "Total: ", length(ams[,"master_peaknum_a"]))
  pairs = paste(sep=".", ams[,"master_peaknum_a"], ams[,"master_peaknum_b"], ams[,"peaknum_a"], ams[,"peaknum_b"])
  unique_pairs = unique(pairs)
  cat("\nWorking on ",length(unique_pairs), "total:\n")
  tmp = llply(1:length(unique_pairs), .progress="text", function(i) {
    
    pair = unique_pairs[i]
    
    peaks = ams[pairs == pair,,drop=F]
    
    charge = min(peaks[,"p_charge_a"])
    
    return(peaks[charge == peaks[,"p_charge_a"],,drop=F])
    })
  
  foo = matrix(unlist(tmp), ncol=ncol(tmp[[1]]), byrow=T, dimnames=list(NULL, colnames(tmp[[1]])))
  cat("\nAfter selecting lowest charge. Unique: ", length(unique(foo[,"master_peaknum_a"])), "Total: ", length(foo[,"master_peaknum_a"]), "\n")
  
  foo
}

filterMaxo = function(
  ams, 
  peaks_a, 
  peaks_b, 
  maxo_r_a, 
  maxo_r_b, 
  maxo_r_fs, 
  maxo_r_fc
  ) {
  ams = as.matrix(ams)
  peaks_a = as.matrix(peaks_a)
  peaks_b = as.matrix(peaks_b)
  
  lookup_a = rep(NA, max(peaks_a[,"peaknum"])); for(i in 1:nrow(peaks_a)) { lookup_a[peaks_a[i,"peaknum"]] = i }
  lookup_b = rep(NA, max(peaks_b[,"peaknum"])); for(i in 1:nrow(peaks_b)) { lookup_b[peaks_b[i,"peaknum"]] = i }
  
  cat("\nBefore maxo filter. Unique: ", length(unique(ams[,"master_peaknum_a"])), "Total: ", length(ams[,"master_peaknum_a"]))
  cat("\nWorking on quad of", nrow(ams), "total:\n")
  tmp = lapply(1:nrow(ams), function(j) { # Each Peak Quad
    cat("\r", j); flush.console();
    
    quad = ams[j,,drop=F]
    
    a12 = peaks_a[lookup_a[quad[,"master_peaknum_a"]],,drop=F]
    a13 = peaks_a[lookup_a[quad[,"peaknum_a"]],,drop=F]
    b12 = peaks_b[lookup_b[quad[,"master_peaknum_b"]],,drop=F]
    b13 = peaks_b[lookup_b[quad[,"peaknum_b"]],,drop=F]
    
    a_maxo_r = a12[,"maxo"] / a13[,"maxo"]
    b_maxo_r = b12[,"maxo"] / b13[,"maxo"] 
    
    bar = list()
    bar$maxo_r_a = a_maxo_r
    bar$maxo_r_b = b_maxo_r
      
    test_a = maxo_r_a/maxo_r_fs < a_maxo_r & a_maxo_r < maxo_r_a*maxo_r_fs
    test_b = maxo_r_b/maxo_r_fs < b_maxo_r & b_maxo_r < maxo_r_b*maxo_r_fs
    test_ab = maxo_r_a/maxo_r_b/maxo_r_fc < a_maxo_r/b_maxo_r & a_maxo_r/b_maxo_r < maxo_r_a/maxo_r_b*maxo_r_fc
    
    foo = test_a & test_b & test_ab
    
    bar$tf = foo
    
    return(bar)
    })
  
  yn = sapply(tmp, function(x) {x$tf})
  a_maxo_r = sapply(tmp, function(x) {x$maxo_r_a})
  b_maxo_r = sapply(tmp, function(x) {x$maxo_r_b})
  
  hum = ams[yn, ,drop=F]
  cat("\nAfter maxo filter. Unique: ", length(unique(hum[,"master_peaknum_a"])), "Total: ", length(hum[,"master_peaknum_a"]), "\n")
  
  cbind(hum, a_maxo_r = a_maxo_r[yn], b_maxo_r = b_maxo_r[yn])
}

prep_peaktable = function(
  peak_table
  ) {
  peak_table = as.matrix(peak_table)
  peak_table = cbind(peaknum = 1:nrow(peak_table), peak_table)
  peak_table = peak_table[order(peak_table[,"rt"]),,drop=F]
  peak_table
}

filterIsos = function(
  pms
  ) {
  cat("\nBefore CAMERA isotope removal. Unique: ", length(unique(pms[,"master_peaknum_a"])), "Total: ", length(pms[,"master_peaknum_a"]), "\n")
  pms = pms[((is.na(pms[,"ciso"]) |  pms[,"ciso"] < 1)),]
  cat("\nAfter CAMERA isotope removal. Unique: ", length(unique(pms[,"master_peaknum_a"])), "Total: ", length(pms[,"master_peaknum_a"]), "\n")
  pms
}


addIsoAdd = function(
  mf_s, 
  an
) { cat("\nCompiling adduct and isotopic data from CAMERA analysis.\n")
    pns = mf_s[,"master_peaknum_a"]
    ns = matrix(NA, nrow = nrow(mf_s), ncol=5, dimnames = list(NULL,c("psg", "cacharge", "carule", "ciso", "ccharge")))
    
    iso_add = l_ply(which(!sapply(an@isotopes,is.null)), function(group) {
      possible_peaks = an@xcmsSet@groupidx[[group]]
      peaks = possible_peaks[possible_peaks %in% pns]
      if (length(which(pns %in% peaks)) < 1) { return(NULL) }
      
      i = an@isotopes[[group]]
      
      n = length(which(pns %in% peaks))
      ns[which(pns %in% peaks),
         c("ciso", "ccharge")
         ] = c(rep(i$iso,n), rep(i$charge, n))
      ns <<- ns
    })
    
    add_add = l_ply(which(!sapply(an@derivativeIons,is.null)), function(group) {
      possible_peaks = an@xcmsSet@groupidx[[group]]
      peaks = possible_peaks[possible_peaks %in% pns]
      if (length(which(pns %in% peaks)) < 1) { return(NULL) }
      
      d = an@derivativeIons[[group]]
      
      n = length(which(pns %in% peaks))
      ns[which(pns %in% peaks),
         c("carule", "cacharge")
         ] = c(rep(d[[1]]$rule_id, n), rep(d[[1]]$charge, n))
      ns <<- ns
    })
    
    ps_add = l_ply(1:length(an@pspectra), function(i) {
      possible_peaks = an@pspectra[[i]]
      peaks = possible_peaks[possible_peaks %in% pns]
      if (length(which(pns %in% peaks)) < 1) { return(NULL) }
      
      ns[which(pns %in% peaks),
         c("psg")
         ] = rep(i, length(which(pns %in% peaks)))
      ns <<- ns
    })
    
    cbind(mf_s, ns)
}


removeDuplicatePeakPicks = function(
  mf_p
  ) {
  ddply(as.data.frame(mf_p), "master_peaknum_a", function(peaks) {
    
    if(length(unique(peaks[,"p_carbons_a"])) == 1 & length(unique(peaks[,"p_charge_a"])) == 1) {
      return(peaks[1,,drop=F])
    } else {
      return(peaks)
    }
    
  })
}

filterChargeSupport = function(
  matches
) { cat("\nLooking for consensus between manual charge prediction and CAMERA.\n")
  ddply(as.data.frame(matches), "master_peaknum_a", .progress="text", function(peaks) {
    if (nrow(peaks) <= 1) { return(peaks) }

    charge_supporta = sapply(1:max(peaks[,"p_charge_a"]), function(x) {
      sum(peaks[x==peaks[,"p_charge_a"],"charge_support_a"] > 0)
    })
    charge_supportb = sapply(1:max(peaks[,"p_charge_a"]), function(x) {
      sum(peaks[x==peaks[,"p_charge_a"],"charge_support_b"] > 0)
    })

    camera_votes=lapply(1:100, function(x){0})
    camera_votes[101] = 1
    for (y in c("ccharge", "cacharge")) {
      if (y %in% colnames(peaks)) { 
        if(!is.na(peaks[1,y])) { 
          if (peaks[1,y] != 0) {
            camera_votes[[abs(peaks[1,y])]] = camera_votes[[abs(peaks[1,y])]]  + 1
          }
        }
      } 
    }
    if (max(charge_supporta) != 0) { 
      for (z in which(max(charge_supporta) == charge_supporta)) {
        camera_votes[[z]] = camera_votes[[z]]  + 1
      }
    }
    if (max(charge_supportb) != 0) {
      for (z in which(max(charge_supportb) == charge_supportb)) {
        camera_votes[[z]] = camera_votes[[z]]  + 1
      }
    }
    
    camera_consensus = which(max(unlist(camera_votes)) == unlist(camera_votes))           
    
    if (length(camera_consensus) > 1 || camera_consensus != 101) {
      the_peaks = which(abs(peaks[,"p_charge_a"]) %in% camera_consensus)
    } else {
      the_peaks = complete.cases(peaks[,"p_charge_a"])
    }
    
    return(peaks[the_peaks,,drop=F])
  })
}

filterPeakSupport = function(mf_pd) {
  ddply(as.data.frame(mf_pd), "master_peaknum_a", .progress="text", function(peaks) {
    
    votes = peaks[,"supporting_peaks_a"] + peaks[,"supporting_peaks_b"]
    
    peaks[which.max.nm(votes),,drop=F]
  })
}

which.max.nm = function(x) {
 which(x == max(x, na.rm=T)) 
}

cNumMatters = function(mf_pd) {
  ddply(as.data.frame(mf_pd), "master_peaknum_a", .progress="text", function(peaks) {
    
    if(nrow(peaks) <= 1) { return(peaks) }
    if (length(unique(peaks[,"p_carbons_a"])) == 1) { return(peaks[1,,drop=F]) }
    return(peaks)
    })
}


selectBest = function(
  mf_pdc
  ) {
  #pick only well deifned peaks
  counts= sapply(1:nrow(mf_pdc), function(i) {
    sum(mf_pdc[,"master_peaknum_a"] == mf_pdc[i,"master_peaknum_a"])
  })
  
  mf = mf_pdc[counts==1,,drop=F]
  mf_bad = mf_pdc[counts>1,,drop=F]
  
  cat("\nFinal tally. Good peaks: ", nrow(mf), " Bad peaks: ", nrow(mf_bad), "Unique bad peaks: ", length(unique(mf_bad[,"master_peaknum_a"])), "\n")
  return(list(good_peaks = mf, bad_peaks = mf_bad)) 
}


heaviestMatters = function(m) {
  ddply(as.data.frame(m), "master_peaknum_a", .progress="text", function(peaks) {
    peaks[peaks[,"p_carbons_a"] == max(peaks[,"p_carbons_a"]),,drop=F]
  })
}
