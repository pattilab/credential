credential = function(
  an,
  r_12t13_a,
  r_12t13_b, 

  isotope_rt_delta_s = 5,
  ppm_for_isotopes = 4,# A number or a function(mz); Takes a mass and returns a +/- ppm (2*ppm would be the window). Must be vectorized

  mixed_ratio_factor = 4,
  mixed_ratio_ratio_factor = 1.8,
  
  mpc_f= 1.1,
  write_files = T
  ) {
  

  if (any(!file.exists(an@xcmsSet@filepaths))) {
    stop("Raw data not found: xsAnnotate@xcmsSet@filepaths must point to the raw data.")
  }
  
  if(is.numeric(ppm_for_isotopes)) {
    ppm = ppm_for_isotopes
    ppm_for_isotopes = function(mz) { return(rep(eval(ppm),length(mz))) }
  }
  
  data(mpc)
  
  sample = an@xcmsSet@peaks[,"sample"]
  peaks_a = prep_peaktable(an@xcmsSet@peaks[sample==1,,drop=F])
  peaks_b = prep_peaktable(an@xcmsSet@peaks[sample==2,,drop=F])
  
  #Generate pairwise matches of isotopic mass within each sample
  pwms_a = pwms(peaks_a, 
                isotope_rt_delta_s = isotope_rt_delta_s, 
                ppm_for_isotopes = ppm_for_isotopes,
                mpc = mpc,
                mpc_f = mpc_f
                )

  pwms_b = pwms(peaks_b, 
                isotope_rt_delta_s =isotope_rt_delta_s, 
                ppm_for_isotopes =ppm_for_isotopes,
                mpc = mpc,
                mpc_f = mpc_f
                )

  cat("\nBuilding alignment index between sample A and sample B based on supplied grouping.\n")
  aligns = buildAlignIndex(an@xcmsSet)

  # Possible alternative, manually group and match.  Probably a good idea.
  #aligns = align(
  #  peaks_a[has_match_a,,drop=F], 
  #  peaks_b[has_match_b,,drop=F],
  #  ppm=30, 
  #  drt=30
  #  )

  
  #Look for patterns of credentialed natural isotopes
  pwms_af = buildChargeSupport(pwms_a)
  pwms_bf = buildChargeSupport(pwms_b)
  
  #Find all putative matches which were replicated aross samples.
  matches = matchAlignsCarbonCharge(aligns, pwms_af, pwms_bf)

  # Adds from xsAnnotate: ciso_a, ciso_b, ccharge_a, ccharge_b, cacharge_a, cacharge_b, crule_a, crule_b
  m_i = addIsoAdd(matches, an) 
  
  m_i2 = filterMpc(m_i, mpc_f)
  
  m_ic = filterChargeSupport(m_i2)
  
  m_icm = filterMaxo(
    m_ic, 
    peaks_a,
    peaks_b,
    maxo_r_a = r_12t13_a,
    maxo_r_b = r_12t13_b,
    maxo_r_fs = mixed_ratio_factor,
    maxo_r_fc = mixed_ratio_ratio_factor
  )
        
  #Remove camera isotopes
  m_icmi = filterIsos(m_icm)
  
  # Remove very similar peaks found by the peakfinding (errors?)
  m_icmid = removeDuplicatePeakPicks(m_icmi)
      
  #Filters based on isotope support
  m_icmidp = filterPeakSupport(m_icmid)
      
  # Last resort, picks the lowest charge
  m_icmidpc = pickLowestCharge(m_icmidp)
       
  #Last last resort. Don't care as long as C# is correct
  m_icmidpcn = cNumMatters(m_icmidpc)
       
  #Last last last resort. Pick heaviest
  m_icmidpcnh = heaviestMatters(m_icmidpcn)
        
  #########################################
  #Any filtering before this step is good - Just discard inderminant peaks
  mf_lst = selectBest(m_icmidpcnh)
  mf = mf_lst$good_peaks
  
  cfs = formatPeaks(mf, peaks_a)
       
  
  bcfs = formatPeaks(mf_lst$bad_peaks, peaks_a)

  if (write_files) {
    
    prepath = paste(sep="_", "credential", round(as.numeric(Sys.time())))  
    dir.create(prepath, showWarnings = FALSE)
    dump("ppm_for_isotopes", paste0(prepath,"/ppm_for_isotopes_dump.R"))
    
    
    csum = paste(llply(
      list( c("r_12t13_a", r_12t13_a),
            c("r_12t13_b", r_12t13_b),
            c("isotope_rt_delta_s", isotope_rt_delta_s),
            c("ppm_for_isotopes", "See ppm_for_isotopes_dump.R"),
            c("mixed_ratio_factor", mixed_ratio_factor),
            c("mixed_ratio_ratio_factor", mixed_ratio_ratio_factor),
            c("mpc_f", mpc_f)
      ), 
      paste, collapse=" - "), collapse="\n")
    
    csum = c(csum, paste(sep=" ","\n\nA - peaks:",nrow(peaks_a), "Pairwise matches: ",sum(sapply(pwms_a, function(x) { nrow(x) }), na.rm=T), "Has a match: ", sum(sapply(pwms_a, nrow) > 0)))
    csum = c(csum, paste(sep=" ","B - peaks:",nrow(peaks_b), "Pairwise matches: ",sum(sapply(pwms_b, function(x) { nrow(x) }), na.rm=T), "Has a match: ", sum(sapply(pwms_b, nrow) > 0)))
    csum = c(csum, paste(sep=" ","\n\nUnique peaknum_a aligns: ", length(unique(aligns[,"peaknum_a"])), ". Total aligns: ", nrow(aligns)))
    csum = c(csum,  paste(sep=" ","After combination and carbon/charge filtering between samples. Unique peaknum_a aligns: ", length(unique(matches[,"master_peaknum_a"])), ". Total possibilities: ", nrow(matches)))
    csum = c(csum,  paste(sep=" ","After filtering based on maxo ratios. Unique peaknum_a aligns: ", length(unique(m_icm[,"master_peaknum_a"])), ". Total possibilities: ", nrow(m_icm)))
    csum = c(csum,  paste(sep=" ","After filtering based on CAMERA isotopes. Unique peaknum_a aligns: ", length(unique(m_icmi[,"master_peaknum_a"])), ". Total possibilities: ", nrow(m_icmi)))
    csum = c(csum,  paste(sep=" ","After filtering peakpicking duplicates. Unique peaknum_a aligns: ", length(unique(m_icmid[,"master_peaknum_a"])), ". Total possibilities: ", nrow(m_icmid)))
    csum = c(csum,  paste(sep=" ","After filtering based on isotope support. Unique peaknum_a aligns: ", length(unique(m_icmidp[,"master_peaknum_a"])), ". Total possibilities: ", nrow(m_icmidp)))
    csum = c(csum,  paste(sep=" ","After filtering based on charge support. Unique peaknum_a aligns: ", length(unique(m_ic[,"master_peaknum_a"])), ". Total possibilities: ", nrow(m_ic)))
    csum = c(csum,  paste(sep=" ","After filtering for lowest charge state. Unique peaknum_a aligns: ", length(unique(m_icmidpc[,"master_peaknum_a"])), ". Total possibilities: ", nrow(m_icmidpc)))
    csum = c(csum,  paste(sep=" ","After filtering by assuming the common carbon number is correct (bad filter). Unique peaknum_a aligns: ", length(unique(m_icmidpcn[,"master_peaknum_a"])), ". Total possibilities: ", nrow(m_icmidpcn)))
    csum = c(csum,  paste(sep=" ","After filtering by picking the heaviest match (bad filter). Unique peaknum_a aligns: ", length(unique(m_icmidpcnh[,"master_peaknum_a"])), ". Total possibilities: ", nrow(m_icmidpcnh)))
    csum = c(csum,  paste(sep=" ","\n\nGood credentialed features (credentialed_features.csv): ", nrow(cfs)))
    csum = c(csum,  paste(sep=" ","Bad credentialed features (credentialed_features_bad.csv): ", nrow(bcfs)))
    
    write.csv(cfs, paste0(prepath,"/credentialed_features.csv"), row.names=F)
    write.csv(mf, paste0(prepath,"/credentialed_features_raw_data"), row.names=F)
    writeLines(csum, paste0(prepath,"/summary_of_credentialing_steps.txt"))
    write.csv(bcfs, paste0(prepath,"/credentialed_features_indeterminant.csv"), row.names=F)
    
    save("an", file=paste0(prepath,"/xsAnnotate.Rdata"))
    save("mpc", file=paste0(prepath,"/mpc.Rdata"))
    
    #Filtered hist(maxo_r)
    f_maxors = list(
      
      ggplot(mf, aes(x=log10(a_maxo_r/b_maxo_r), xintercept=log10(r_12t13_a/r_12t13_b))) +
        geom_histogram(fill="white", colour="black", binwidth=0.0025) + xlim(-0.1,0.1) +
        geom_vline() +
        ggtitle("Maxo Ratio of Sample A/B; Filtered"),
      
      ggplot(mf, aes(x=log10(a_maxo_r), xintercept=log10(r_12t13_a))) +
        geom_histogram(fill="white", colour="black", binwidth=0.015) + xlim(-0.5,0.5) +
        geom_vline() +
        ggtitle("Maxo Ratio of Sample A; Filtered"),
      
      ggplot(mf, aes(x=log10(b_maxo_r), xintercept=log10(r_12t13_b))) +
        geom_histogram(fill="white", colour="black", binwidth=0.015) + xlim(-0.5,0.5) +
        geom_vline() +
        ggtitle("Maxo Ratio of Sample B; Filtered")
      
    )
    
    
    #build unfiltered
    unf = cbind(m_ic,
                a_maxo_r = peaks_a[m_ic[,"master_peaknum_a"],"maxo"] / peaks_a[m_ic[,"peaknum_a"],"maxo"],
                b_maxo_r = peaks_b[m_ic[,"master_peaknum_b"],"maxo"] / peaks_b[m_ic[,"peaknum_b"],"maxo"]
    )
    #Unfiltered hist(maxo_r)
    unf_maxors = list(
      
      ggplot(unf, aes(x=log10(a_maxo_r/b_maxo_r), xintercept=log10(r_12t13_a/r_12t13_b))) +
        geom_histogram(fill="white", colour="black", binwidth=0.1) + xlim(-5,5) +
        geom_vline() +
        ggtitle("Maxo Ratio of Sample A/B"),
      
      ggplot(unf, aes(x=log10(a_maxo_r), xintercept=log10(r_12t13_a))) +
        geom_histogram(fill="white", colour="black", binwidth=0.1) + xlim(-5,5) +
        geom_vline() +
        ggtitle("Maxo Ratio of Sample A"),
      
      ggplot(unf, aes(x=log10(b_maxo_r), xintercept=log10(r_12t13_b))) +
        geom_histogram(fill="white", colour="black", binwidth=0.1) + xlim(-5,5) +
        geom_vline() +
        ggtitle("Maxo Ratio of Sample B")
      
    )
    
    
    raw = cbind(mf, an@xcmsSet@peaks[mf[,"master_peaknum_a"],])
    
    pdf(paste0(prepath,"/credentialed_maxo_graphic.pdf"), width=10, height = 12)
    #ppm vs maxo
    ggplot(raw, aes(x=mz, y=p_iso_ppm_a, colour=log10(maxo))) + 
      geom_point() +
      scale_colour_continuous(low="white", high="darkblue") + 
      ggtitle("ppm Error Comparing U12C to U13C")
    
    do.call("grid.arrange", unf_maxors)
    
    do.call("grid.arrange", f_maxors)
    
    dev.off()
  }
  
  return(cfs)
}