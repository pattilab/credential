.onLoad = function(libname, pkgname) {
  data("atm", package=pkgname, envir=parent.env(environment()))
}

credential = function(
  an,
  r.a,
  r.b,
  mzdiff = atm$c13-atm$c12,
  rt.lim = 10,
  rcorr.lim = 0.7,
  ppm.lim = 2.5,
  mpc.f= 1.1,
  charges = c(1,2,3,4,5,6),
  g=100
  ) {
  
  # Prepare data

  if (any(!file.exists(an@xcmsSet@filepaths))) { stop("Raw data not found: xsAnnotate@xcmsSet@filepaths must point to the raw data.") }
  
  data(mpc)
  data(atm)
  mpc = within(mpc, min_mpc <- min_mpc / mpc.f)
  mpc = within(mpc, max_mpc <- max_mpc * mpc.f)
  
  cat("Loading raw files. ")
  xr.a = xcmsRaw(an@xcmsSet@filepaths[1], profstep = 10)
  xr.b = xcmsRaw(an@xcmsSet@filepaths[2], profstep = 10)
  
  p = an@xcmsSet@peaks
  
  group = rep(NA, nrow(p))
  for (i in 1:length(an@xcmsSet@groupidx)) {
    group[an@xcmsSet@groupidx[[i]]] = i
  }
  
  psg = rep(NA, nrow(p))
  for (i in 1:length(an@pspectra)) {
    psg[an@pspectra[[i]]] = i
  }
  
  pn = 1:nrow(an@xcmsSet@peaks)
  an@xcmsSet@peaks = cbind(p, pn, group, psg)
  p = data.frame(an@xcmsSet@peaks)
  
  
  # Build groups of peaks which could be isotopes of each other
  
  iso.a = buildIsoGroups(
    peaks = p[p[,"sample"] == 1,],
    xr = xr.a,
    rt.corr = an@xcmsSet@rt$corr[[1]],
    ppm.lim = ppm.lim,
    rcorr.lim = rcorr.lim,
    rt.lim =   rt.lim,
    mzdiff = mzdiff,
    charges = charges,
    g=g
    )
  
  iso.b = buildIsoGroups(
    peaks = p[p[,"sample"] == 2,],
    xr = xr.b,
    rt.corr = an@xcmsSet@rt$corr[[2]],
    ppm.lim = ppm.lim,
    rcorr.lim = rcorr.lim,
    rt.lim = rt.lim,
    mzdiff = mzdiff,
    charges = charges,
    g=g
  )
  
  cat("Merging peaks. ")
  riso.g = replicatedGroups(an, iso.a, iso.b) # Based on supplied grouping in an@xs@groupidx
  cam.g = compileCamera(riso.g$group, an)
  
  risoc.g = merge(cam.g, riso.g, by = "group", all=T)
  peak.gs.a = merge(risoc.g, an@xcmsSet@peaks, by.x=c("pn.a", "group"), by.y=c("pn","group"), all.x=T)
  peak.gs.ab = merge(peak.gs.a, an@xcmsSet@peaks, by.x=c("pn.b", "group"), by.y=c("pn","group"), all.x=T, suffixes = c(".a", ".b"))

  pgs.f = within(peak.gs.ab, mz <- mergeValues(mz.a,mz.b))
  pgs.f = within(pgs.f, maxo <- mergeValues(maxo.a,maxo.b))
  
  cat("Finding pairs of credentialed (U12-U13) peaks in the ", length(unique(pgs.f[,"isog"]))," isotope groups. \n")
  pwms = ddply(pgs.f, "isog", .progress="text", function(x) {
    credPairs(x, ppm.lim, mzdiff, mpc)
    })

  
  # Save results
  cat("Saving results. ")
  prepath = paste(sep="_", "credential", round(as.numeric(Sys.time())))
  dir.create(prepath, showWarnings = FALSE)
  
  write.csv(pwms, paste0(prepath,"/pwms.csv"), row.names=F)
  write.csv(pgs.f, paste0(prepath,"/pgs.f.csv"), row.names=F)
  
  write.csv(iso.a, paste0(prepath,"/iso.a.csv"), row.names=F)
  write.csv(iso.b, paste0(prepath,"/iso.b.csv"), row.names=F)
  
  save("an", file=paste0(prepath,"/xsAnnotate.Rdata"))
  save("mpc", file=paste0(prepath,"/mpc.Rdata"))
  f = as.list(match.call()); write.csv(cbind(names(f), unlist(f)), paste0(prepath,"/arguments.csv"), row.names=F, col.names=F, sep="\t")
  
  plotCfOutput(prepath = prepath, pwms, pgs.f)
    
  return(pwms)
}