replicatedGroups = function(an, isog.a, isog.b) {
  isog.ad = dcast(melt(isog.a), ... ~ Var2)[,-1]
  names(isog.ad)[1] = "isog"
  isog.bd = dcast(melt(isog.b), ... ~ Var2)[,-1]
  names(isog.bd)[1] = "isog"
  
  m.a = merge(an@xcmsSet@peaks[,c("pn", "group")], isog.ad, by="pn")
  m.b = merge(an@xcmsSet@peaks[,c("pn", "group")], isog.bd, by="pn")
  
  risoc.g = merge(m.a, m.b, by="group", all=T, suffixes = c(".a", ".b"))
  
  risoc.g = cbind(risoc.g, isog = mergeFactors(risoc.g[,"isog.a"], risoc.g[,"isog.b"]))
  risoc.g = cbind(risoc.g, charge = mergeFactors(risoc.g[,"charge.a"], risoc.g[,"charge.b"]))
  
  ddply(risoc.g, "isog", function(x) {
    cbind(x, seq = mergeFactors(x[,"seq.a"], x[,"seq.b"]))
  })
}

mergeValues = function(v1, v2) {
  v.new = v1
  v.new[is.na(v.new)] = v2[is.na(v.new)]
  v.new
}

mergeFactors = function(v1, v2) {
  v.new = v1
  na.v1 = which(is.na(v.new))
  replacement = llply(na.v1, function(y) {
    v2.val = which(v2 %in% v2[y])
    unique(na.omit(v1[v2.val]))
  })
  for (i in seq_along(na.v1)) {
    if (length(replacement[[i]]) < 1) { next; }
    if (length(replacement[[i]]) < 1) { warning("Vectors not merged cleanly."); }
    v.new[na.v1[i]] = replacement[[i]]
  }
  na.new = which(is.na(v.new))
  v.new[na.new] = v2[na.new]*1000
  v.new = as.numeric(factor(v.new))
  v.new
}


compileCamera = function(group, an) { # pwms must contain pn12 to subset groups on
  d = ddply(data.frame(group), "group", function(x) {
    d = an@derivativeIons[[x[1,1]]]
    if (!is.null(d)) {
      return(data.frame(cam.d=d[[1]], cam.d.n=length(d)))
    }
  })
  
  i = ddply(data.frame(group), "group", function(x) {
    d = an@isotopes[[x[1,1]]]
    if (!is.null(d)) {
      return(data.frame(cam.i=d))
    }
  })
  
  merge(d, i, by="group", all=T)
}

