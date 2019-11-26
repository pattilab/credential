#' Step 1 of Credentialing: Find Isotopic Knots and Resolve Charge State
#'
#' @description This function searches for 'knot' from feature table based on retention time and mass residual analysis. 
#' A knot is a set of isotopologue-like features with defined m/z and retention time difference and same charge states.
#' Based on the fact that isotopologues have identical mass residual against \code{cd} = 13C-12C (1.00335), mass residual 
#' analysis can be used to cluster isotopologues. Mass residual of a feature is calculated as: (mz %% cd). 
#' @usage findknots(features, .zs=1:4, ppmwid=4, rtwid = 1, cd = 13.00335-12) 
#' @param features data.table A feature table including at least four exclusive columns: "cc", "mz", "rt" and "i". Each row is a feature (peak).
#' @param .zs integer A vector indicating all possible charge states to be searched against.
#' @param ppmwid numeric Mass error tolerance in ppm.
#' @param rtwid numeric The maximum retention time difference allowed in seconds.
#' @param cd numeric Unit mass difference between unlabeled and labeled atoms. Used for mass residual calculation. Defalut is 13C-12C = 13.00335 - 12.
#' @import data.table utils 
#' @return list A list with two tables "cc_knot" and "knot".  \code{cc_knot} data.table The index of feature-to-knot 
#' assignment. \code{knot} data.table Aggregated information of all knots: knot#, meanr (mean mass residual), meanmz 
#' (mean m/z), basemz (lowest mz peak), rt (mean rt), maxi (maximum intensity), n (number of peaks), z (charge state), 
#' dir (relative of the knot from left to right, positive value for increasing, negative for decreasing).
#' @keywords credential findknots
#' @seealso \code{\link{credentialing}}
#' @export
#' 
findknots = function(features, .zs=1:4, ppmwid=4, rtwid = 1, cd = 13.00335-12) { # Should search for and return isotope knots (including z with one knot per peak) representing individual compounds
  factor = 1
  Cc = copy(features)
  
  # Find putative knots in various charge states
  zknots = lapply(.zs, function(.z) {
    cat("\rFinding isotope knots. Charge state:", .z)
    findzknots(Cc, .z, ppmwid = ppmwid, rtwid = rtwid)
  })
  names(zknots) = .zs
  
  # Aggregate
  knots = do.call(what=rbind, lapply(zknots, function(x) {
    Cc[x$c13[x$c13.,,on="c13"],,on="cc"]
  }))
  
  # Resolve peak knot charge state assignment
  knots$ind = as.numeric(factor(knots$cc))
  counts = matrix(nrow = length(unique(knots$cc)), ncol = length(.zs))
  ass= cbind(knots$ind, knots$z)
  counts[ass] = knots$n
  
  keeps = which(counts == matrixStats::rowMaxs(counts), arr.ind=T)
  keeps = keeps[!duplicated(keeps[,1]), ]
  
  keepme= match(paste(keeps[,1], keeps[,2]), paste(ass[,1], ass[,2]))
  
  temp = knots[keepme][,c13 := as.numeric(factor(paste(c13, z)))]
  
  #Calculate Direction
  annotatetails = function(ps) {
    o = order(ps$mz)
    i = ps$i[o]
    term = which(diff(diff(i/max(i))) < -0.6)[1] + 1
    
    if (!is.na(term) & term > 1) {
      c(rep(F, term), rep(T, nrow(ps) - term))[order(o)]
    } else {
      rep(F, nrow(ps))
    }
  }
  
  temp[, tail := annotatetails(data.frame(mz, i)), by="c13"]
  
  calcdir = function(ps) {
    ps = subset(ps, tail == F)

    f = tail(seq_len(nrow(ps)), n=3)
    b = head(seq_len(nrow(ps)), n=3)
    
    ds = c(
      mean(diff(ps$i[order(ps$mz)][f]))/max(ps$i),
      mean(diff(ps$i[order(ps$mz)][b]))/max(ps$i)
    )
    
    ds[which.max(abs(ds))]
  }
  
  #Aggregate
  c13. = temp[,.(meanr = mean(mz %% (cd/z)), meanmz = mean(mz), basemz = min(mz), mainmz = mz[which.max(i)], rt = mean(rt), maxi=max(i),  n = length(mz), dir = calcdir(data.frame(mz, i, tail)), z= z[1]),by=c13]
  c13.[n==1, z := 0]
  c13 = temp[,.(cc, c13, tail)]
  
  names(c13)[2] = "knot"
  names(c13.)[1] = "knot"
  Knot = list(cc_knot = c13, knot = c13.)
  
  # Resolve merged knots
  Knot = fixmergedknots(Knot,features)
  
  cat("\nFound", nrow(Knot$knot), "isotope knots.")
  
  return(Knot)
}

#' Core function of findknots(): Searching knots at charge state z
#' 
#' @description This function is a subset function of findknots(). It searches isotopologue-like features at charge state z.
#' @usage findzknots(features, .z=1, ppmwid=5, rtwid = 1, cd = 13.00335-12)
#' @param features data.table A feature table including at least four exclusive columns: "cc", "mz", "rt" and "i". Each row is a feature (peak).
#' @param .z integer Charge state of the features.
#' @param ppmwid numeric Mass error tolerance in ppm.
#' @param rtwid numeric The maximum retention time difference allowed in seconds.
#' @param cd numeric Unit mass difference between unlabeled and labeled atoms. Defalut is 13C-12C = 13.00335 - 12.
#' @import igraph magrittr stats
#' @return list A list with values "cc_knot" and "knot".  \code{cc_knot} contains the features assignments to a knot. \code{knot} contains aggregate information about each knot.
#' @seealso \code{\link{findknots}}
#'

findzknots = function(features, .z=1, ppmwid=5, rtwid = 1, cd = 13.00335-12) {
  
  Cc = copy(data.table(features))
  factor = 1
  scales = c(ppmwid * 700 / 1E6, rtwid)
  
  # Cache residual
  Cc[,c13r := mz %% (cd/.z)]
  
  # Iterative splitting
  gcols = paste0("g", 1:6)
  Cc[,(gcols) := "1"]
  ngs = 1
  repeat{
    Cc[,g4 := paste(g4, clustgroup(cbind(c13r, rt), scales, factor)), by=gcols]
    Cc[,g3 := paste(g3, breakgroup((mz %/% (cd/.z)), 1)), by=gcols]
    
    l = length(unique(Cc[,do.call(paste, .SD), .SDcols = gcols])); if (l == ngs) break else ngs = l
  }
  Cc[,c13c := as.integer(factor(paste(do.call(paste, .SD)))), .SDcols = gcols]
  
  c13c_meta = Cc[,.(meanr = mean(mz %% (cd/.z)), rt = mean(rt), n = length(mz), c13 = c13c[1], z= .z), by="c13c"][, c13c:=NULL]
  
  return(list(c13 = Cc[,.(cc, c13 = c13c)], c13. = c13c_meta))
}
