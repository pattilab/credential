#' Find knots - sets of features with defined spacing and retention time.
#'
#' @param features data.frame Rows are features. Columns are "cc"
#' @param .z integer Charge state to assume
#' @param ppmwid numeric The maximum mass error in ppm.
#' @param rtwid numeric The maximum retention time difference in seconds.
#' @param cd numeric The mass spacing to search for (defaults to C13 - C12)
#' 
#'
#' @return list A list with values "cc_knot" and "knot".  \code{cc_knot} contains the features assignments to a knot. \code{knot} contains aggregate information about each knot.
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
  
  list(c13 = Cc[,.(cc, c13 = c13c)], c13. = c13c_meta)
}


#' Find knots and resolve charge state
#'
#' @param features data.frame Rows are features. Columns are "cc"
#' @param .zs integer A vector containing the absolute value of the possible charges.
#' @param ppmwid numeric The maximum mass error in ppm.
#' @param rtwid numeric The maximum retention time difference in seconds.
#' @param cd numeric The mass spacing to search for (defaults to C13 - C12)
#' @param minlength Integer. If an ROI is shorter than this it is discarded.
#'
#' @seealso \link{\code{findzknots}}
#' 
#' @return list A list with values "cc_knot" and "knot".  \code{cc_knot} contains the features assignments to a knot. \code{knot} contains aggregate information about each knot.
#'
findknots = function(Cc.in, .zs=1:4, ppmwid=4, rtwid = 1, cd = 13.00335-12) { # Should search for and return isotope knots (including z with one knot per peak) representing individual compounds
  factor = 1
  Cc = copy(Cc.in)
  
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
  c13. = temp[,.(meanr = mean(mz %% (cd/z)), meanmz = mean(mz), mainmz = mz[which.max(i)], rt = mean(rt), maxi=max(i),  n = length(mz), dir = calcdir(data.frame(mz, i, tail)), z= z[1]),by=c13]
  c13.[n==1, z := 0]
  c13 = temp[,.(cc, c13, tail)]
  
  names(c13)[2] = "knot"
  names(c13.)[1] = "knot"
  
  cat("\nFound", nrow(c13.), "isotope knots.")
  
  list(cc_knot = c13, knot = c13.)
}
