plotCfOutput = function(prepath = ".", pwms, pgs.f) {
  peaks = merge(pwms, pgs.f, by.x="c12.g", by.y="group")
  peaks2 = merge(peaks, pgs.f,  by.x="c13.g", by.y="group", suffixes=c(".12", ".13"))
  
  
  # Maxo Ratio Plots
  maxo.plots = list(
    ggplot(peaks2, aes(x=log10(maxo.a.12/maxo.a.13))) + 
      geom_histogram(binwidth=.1) + 
      facet_grid(. ~ mult) +
      ggtitle("Maxo Ratios of A") +
      xlim(-1.5,1.5)
    ,
    ggplot(peaks2, aes(x=log10(maxo.b.12/maxo.b.13))) + 
      geom_histogram(binwidth=.1) + 
      facet_grid(. ~ mult) +
      ggtitle("Maxo Ratios of B") +
      xlim(-1.5,1.5)
    ,
    ggplot(peaks2, aes(x=log10(maxo.a.12/maxo.a.13/(maxo.b.12/maxo.b.13)))) + 
      geom_histogram(binwidth=.1) + 
      facet_grid(. ~ mult) +
      ggtitle("Maxo Ratios of A/B") +
      xlim(-1.5,1.5)
  )
  do.call("grid.arrange", maxo.plots)
  
  # ppm Plots
  ppm.plots = list(
    ggplot(peaks2, aes(x=ppm)) + 
      geom_histogram(binwidth=.3) + 
      facet_grid(. ~ detected.12) +
      ggtitle("ppm Error from sample A")
    ,
    ggplot(peaks2, aes(x=ppm, y=mult)) + 
      geom_point() + 
      ggtitle("ppm Error of peaks which could not be determined")
  )
  do.call("grid.arrange", ppm.plots)
  
  sample.plots = llply(sample(pwms[,1], 9), function(x) plotPwms(pwms[x,], pgs.f))
  do.call("grid.arrange", sample.plots)

  qc.plots = list(
    ggplot(pwms, aes(x=mult)) +
             geom_histogram(binwidth=1)
    ,
    ggplot(pwms, aes(x=step))  +
      geom_histogram(binwidth=1)
    ,
    ggplot(pwms, aes(x=detected.12))  +
      geom_histogram(binwidth=1) +
      xlim(1,11)
    ,
    ggplot(pwms, aes(x=detected.13))  +
      geom_histogram(binwidth=1) +
      xlim(1,11)
    ,
    ggplot(pwms, aes(x=factor(is.na(c12.g)))) +
      geom_bar() +
      ggtitle("IsoGroups without a credentialed isotope.")
  )
  
  pdf(paste0(prepath, "/plotCfOutput.pdf"), width=10, height=8) 
  {
    do.call("grid.arrange", qc.plots)
    do.call("grid.arrange", maxo.plots)
    do.call("grid.arrange", ppm.plots)
    do.call("grid.arrange", sample.plots)
  } 
  dev.off()

}




plotPwms = function(pwms.1 = pwms[sample(pwms[,1],1),], pgs.f = pgs.f2) {
  isog = subset(data.frame(pgs.f), isog %in% pwms.1[,"isog"])
  pair = subset(data.frame(isog), group %in% pwms.1["c12.g"] | group %in% pwms.1["c13.g"])

  if(any(is.na(pwms.1))) {
    return(ggplot(data.frame(isog), aes(x = mz, xend=mz, yend=0, y= maxo, colour=factor(seq))) + 
      geom_segment() +
      ggtitle("No matches.")) + 
      theme(legend.position="none")
  }
  
  
  maxmaxo = max(isog[,"maxo"]) * 1.05
  pair2 = data.frame(pair, maxmaxo)
  
  
  ggplot() + 
    geom_segment(data=data.frame(isog, maxmaxo), mapping=aes(x = mz, xend=mz, yend=0, y= maxo, colour=factor(seq))) + 
    geom_point(data=pair2,mapping=aes(x=mz, y=maxmaxo)) + 
    ggtitle(paste0("Isog ", pwms.1[,"isog"], ". Step: ", pwms.1[,"step"], ". 1x", sum(!is.na(pair[1,c("pn.a", "pn.b")])), ". 2x", sum(!is.na(pair[2,c("pn.a", "pn.b")])))) +
    theme(legend.position="none")
}

plotIsog = function(isog, an) {
  ggplot(merge(isog, an@xcmsSet@peaks, by="pn"), aes(x = mz, xend=mz, yend=0, y= maxo, colour=factor(seq))) + geom_segment()
}

plotMasses = function(peaks) {
  ggplot(peaks, aes(x = mz, xend=mz, yend=0, y= maxo, colour=factor(seq))) + geom_segment()
}

agSpecs = function(pn.a = pos.matches[sample(nrow(pos.matches), 1),"pn.a"], ppm=50) {
  
  ag = pn.a
  count = 0
  repeat {
    y = subset(pos.matches, pn.a %in% ag | pn13.a %in% ag | pn.b %in% ag | pn13.b %in% ag & p.ppm.a < ppm)
    ag.n = unique(unlist(y[,c("pn.a", "pn13.a","pn.b", "pn13.b")]))
    
    if (all(ag.n %in% ag)) { print(count); break;}
    count = count + 1
    ag = ag.n
  }
  
  xr.l = vector("list", length=2)
  xr.l[[1]] = xr_a
  xr.l[[2]] = xr_b
  eic.list=llply(ag, nmEIC.an, an, xr.l)
  eics = eic.list[[1]]
  for(x in eic.list[-1]) {
    eics = merge(eics, x, by="scan", all=T)
  }
  
  
  eic.plot = ggplot(melt(eics, "scan"), aes(x = scan, y=value, colour=variable)) +
    geom_path()
  
  spec.plot = ggplot(data.frame(subset(an@xcmsSet@peaks, pn %in% ag)), aes(x = mz, xend=mz, yend=0, y= maxo, colour=factor(pn))) +
    geom_segment() +
    facet_grid(sample ~ .) + 
    ggtitle(paste0(
      paste0(collapse=", ", count(subset(an@xcmsSet@peaks, pn %in% ag)[,"sample"])[,2]),
      paste0(" pn: ", ag[1])
    ))
  
  grid.arrange(spec.plot, eic.plot)
  
  return(y)
}
