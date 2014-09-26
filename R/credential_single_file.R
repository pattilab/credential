credential.one = function(
  xs_a,
  r_12t13_a,
  an_a, 
  
  isotope_rt_delta_s = 5,
  ppm_for_isotopes = 4,
  
  mixed_ratio_factor = 4,
  xs_a_file=NULL, 
  mpc_f= 1.1,
  write_files = T
) {
  
  if (!is.null(xs_a_file)) {xs_a@filepaths = xs_a_file  }
  
  data(mpc)
  #mpc = read.csv(header=T, file="dependencies/mm_mpc.csv")
  mpc = mpc[1:150,,drop=F]
  
  csum=c()
  
  peaks_a = prep_peaktable(xs_a@peaks)
  
  #Generate pairwise matches of isotopic mass within each sample
  pwms_a = pwms(peaks_a, 
                isotope_rt_delta_s = isotope_rt_delta_s, 
                ppm_for_isotopes = ppm_for_isotopes,
                mpc = mpc,
                mpc_f = mpc_f,
                mixed_ratio_12t13 = r_12t13_a,
                mixed_ratio_factor = mixed_ratio_factor)
  has_match_a = sapply(pwms_a, nrow) > 0
  csum = c(csum, paste(sep=" ","A - peaks:",nrow(peaks_a), "Pairwise matches: ",sum(sapply(pwms_a, function(x) { nrow(x) }), na.rm=T), "Has a match: ", sum(has_match_a)))
  
  if(length(pwms_a) < 1) {stop("No pairwise matches found.  Labeled samples? Proper ratio?")}  
  
  #Look for patterns of credentialed natural isotopes
  pwms_af = buildChargeSupport(pwms_a)
  
  pwms = format_pwms(pwms_af)
  
  pwms_afi = addIsoAdd.one(pwms, an_a)
  
  m_ic = filterChargeSupport.one(pwms_afi)
  
  m_icmi = filterIsos.one(m_ic)
  m_icmid = removeDuplicatePeakPicks.one(m_icmi)
  m_icmidp = filterPeakSupport.one(m_icmid)
  m_icmidpc = heaviestMatters.one(m_icmidp)
  
  cfs = formatPeaks.one(m_icmidpc, peaks_a)
  
  write.csv(cfs, "credentialed_features_1.csv")
}

filterChargeSupport.one = function(
  mf_pd
) {
  matrixLapplyUniqueRows(mf_pd, "master_peaknum", function(peaks) {
    
    if(nrow(peaks) <= 1) {return(peaks)}
    
    charge_supporta = sapply(1:max(peaks[,"p_charge"]), function(x) {
      sum(peaks[x==peaks[,"p_charge"],"charge_support"] == 1)
    })
    
    camera_votes=list(0,0,0,0,0,0,0,1)
    for (y in c("ccharge_a", "cacharge_a")) {
      if (y %in% colnames(peaks)) { if(!is.na(peaks[1,y])) { camera_votes[[peaks[1,y]]] = camera_votes[[peaks[1,y]]]  + 1}
      } }
    if (max(charge_supporta) != 0) { 
      for (z in which(max(charge_supporta) == charge_supporta)) {
        camera_votes[[z]] = camera_votes[[z]]  + 1
      }
    }
    
    camera_consensus = which(max(unlist(camera_votes)) == unlist(camera_votes))           
    
    if (length(camera_consensus) > 1 || camera_consensus != 8) {
      the_peaks = which(peaks[,"p_charge"] %in% camera_consensus)
    } else {
      the_peaks = complete.cases(peaks[,"p_charge"])
    }
    
    return(peaks[the_peaks,,drop=F])
  })
}

addIsoAdd.one = function(
  mf_s, 
  an_a
) { #TODO: Incoporate earlier in workflow
  #add_iso_info
  as = lapply(1:nrow(mf_s), function(i) {
    peak = mf_s[i,,drop=F]
    
    isos = an_a@isotopes[[peak[,"master_peaknum"]]]
    if (is.null(isos)) { return(c(NA, NA)) }
    
    return(c(isos$iso, isos$charge))
  })
  as2 = do.call("rbind", as)
  colnames(as2) = c("ciso_a", "ccharge_a")
  
  
  adas = lapply(1:nrow(mf_s), function(i) {
    peak = mf_s[i,,drop=F]
    
    isos = an_a@derivativeIons[[peak[,"master_peaknum"]]]
    if (is.null(isos)) { return(c(NA, NA)) }
    
    return(c(isos$rule_id, isos$charge))
  })
  adas2 = do.call("rbind", as)
  colnames(adas2) = c("crule_a", "cacharge_a")
  
  
  #Add PS info
  pspecs = lapply(1:length(an_a@pspectra), function(i) {
    cbind(master_peaknum = rep(i, length(an_a@pspectra[[i]])), psg = an_a@pspectra[[i]])
  })
  pspecs = do.call("rbind", pspecs)
  pspecs2 = sapply(1:nrow(mf_s), function(x) {
    peak = mf_s[x,,drop=F]
    pspecs[which(pspecs[,"master_peaknum"] == peak[,"master_peaknum"]), "psg"][1]
  })
  
  cbind(mf_s, as2, adas2, pspec = pspecs2)
}

filterIsos.one = function(
  pms
) {
  cat("\nBefore CAMERA isotope removal. Unique: ", length(unique(pms[,"master_peaknum"])), "Total: ", length(pms[,"master_peaknum"]), "\n")
  pms = pms[((is.na(pms[,"ciso_a"]) |  pms[,"ciso_a"] < 1)),]
  cat("\nAfter CAMERA isotope removal. Unique: ", length(unique(pms[,"master_peaknum"])), "Total: ", length(pms[,"master_peaknum"]), "\n")
  pms
}

removeDuplicatePeakPicks.one = function(
  mf_p
) {
  matrixLapplyUniqueRows(mf_p, "master_peaknum", function(peaks) {
    
    if(length(unique(peaks[,"p_carbons"])) == 1 & length(unique(peaks[,"p_charge"])) == 1) {
      return(peaks[1,,drop=F])
    } else {
      return(peaks)
    }
    
  })
}
filterPeakSupport.one = function(mf_pd) {
  matrixLapplyUniqueRows(mf_pd, "master_peaknum", function(peaks) {
    
    votes = peaks[,"supporting_peaks"] + peaks[,"supporting_peaks"]
    
    peaks[which.max.nm(votes),,drop=F]
  })
}
heaviestMatters.one = function(m) {
  matrixLapplyUniqueRows(m, "master_peaknum", function(peaks) {
    peaks[peaks[,"p_carbons"] == max(peaks[,"p_carbons"]),,drop=F]
  })
}

formatPeaks.one = function(
  mf, 
  ps
) {
  c("mz", "rt", "rtmax", "rtmin", "match_mz", "n_carbons", "peaknum_a", "adduct", "iso", "charge")
  
  if (nrow(mf) < 1) {return("")}
  
  tmp = lapply(1:nrow(mf), function(i) {
    peak = mf[i,,drop=F]
    peak_a = ps[which(ps[,"peaknum"] == peak[,"master_peaknum"]),,drop=F]
    match_a = ps[which(ps[,"peaknum"] == peak[,"peaknum"]),,drop=F]
    
    crule = peak[,"crule_a"]
    ciso = peak[,"ciso_a"]
    
    cbind(
      peaknum_a = peak_a[,"peaknum"], 
      mz= peak_a[,"mz"], 
      rt= peak_a[,"rt"], 
      rtmin= peak_a[,"rtmin"], 
      rtmax= peak_a[,"rtmax"],
      match_mz = match_a[,"mz"],
      n_carbons = peak[,"p_carbons"],
      adduct = crule,
      iso = ciso,
      charge = peak[,"p_charge"],
      maxo = peak_a[,"maxo"]
    )
  })
  do.call("rbind", tmp)
}

format_pwms = function(
  pwma
) {
  all_matches = lapply(1:length(pwma), function(i) {  cat("\r", i); flush.console();  #Looks for matches in a, that have matches in b with the same carbon number
    cbind(pwma[[i]], master_peaknum = rep(as.numeric(names(pwma)[i]), nrow(pwma[[i]])))
    })
  hum = do.call("rbind", all_matches, quote=T)
}