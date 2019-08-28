#' Step 3 of Credentialing: Find Isotope Knots and Resolve Charge State
#'
#' @description This function searches 'quipu' in knots. A quipu is a set of knots with same charge state and similar 
#' retention time and m/z spacing. A credentialed quipu also qualifies for theoretical mass-per-cabon range and 
#' ratio of the head-tail peak intensities
#' @usage credentialknots(Knots, ppmwid = 9, rtwid = 1, mpc = c(12, 120), Ratio = 1/1, Ratio.lim = 0.1, 
#' maxnmer = 4, cd = 13.00335-12, .zs = 1:4)
#' @param Knots list The original output of findknots()
#' @param ppmwid numeric Mass error tolerance (ppm) of credentialed knots.
#' @param rtwid numeric Retention time tolerance (seconds) of credentialed knots.
#' @param mpc numeric Mass-per-carbon filter. A vector of two values setting minimum and maximum theorectical mass per carbon.
#' @param Ratio numeric The mixing ratio of unlabled/labeled (12C/13C) samples.
#' @param Ratio.lim numeric Ratio filter. A decimal number (0,1] controling ranges of acceptable deviation of 
#' unlabeled/labeled ratio from mixing ratio.
#' @param maxnmer numeric Maximumly allowed number of knots in each credentialed groups, The default value is 4.
#' @param cd numeric Unit mass difference between unlabeled and labeled atoms. Defalut is 13C-12C = 13.00335 - 12.
#' @param .zs numeric Possible charge states of isotopologues. Default value is 1:4.
#' @keywords credentialing credentialknots 
#' @import data.table igraph magrittr stats matrixStats
#' @seealso \code{\link{findknots}} \code{\link{credentialing}}
#' @return list A list of three tables "knot_quipu", "quipu" and "quipu_stat". \code{knot_quipu} data.table The index of quipu-to-knot 
#' assignment. \code{quipu} data.table Aggregated information of all quipus: quipu#, nknot (number of knots), 
#' npeak (number of peaks), rtmean (mean retention time), basemz (lowest mz), mainmz 1 (credentialed pair - 
#' unlabeled mz). mainmz 2 (credentialed pair - labeled mz), int 1 (credentialed pair - intensity of unlabeled mz). 
#' int 2 (credentialed pair - intensity of labeled mz), ratio (ratio of unlabled/labeled intensity).
#' \code{quipu_stat} data.table Statistical information of all quipus: quipu#, minsupport (minimum number 
#' of peaks among the knots), maxsupport (maximum number of peaks among the knots), nknot
#' @export
#' 

credentialknots = function(Knots, ppmwid = 9, rtwid = 1, mpc = c(12, 120), Ratio = 1/1, Ratio.lim = 0.1, 
                           maxnmer = 4, cd = 13.00335-12, .zs = 1:4) {
  
  Knot = copy(Knots$knot)
  cat("\nCredentialing within", length(unique(Knot$knot)), "supplied knots.")
  factor = 1
  do.plot = F
  maxdimer = maxnmer
  scales = c(ppmwid * 700 / 1E6, rtwid)
  
  # Initial Grouping by rt and meanr and z
  gcols = paste0("g", 1:6)
  Knot[,(gcols) := 1L]
  
  Knot[,g1 := z]
  Knot[,g2 := as.integer(clustgroup(cbind(meanr, rt), scales, factor)), by=gcols]
  
  Knot[,rrtg := as.integer(factor(paste(do.call(paste, .SD)))), .SDcols = gcols]
  Knot[,(gcols) := NULL]
  
  cat("\nWorking with supported charge states.")
  Knot_quipu = copy(Knot[,.(knot)])
  for (.rrtg in unique(Knot[z>0]$rrtg)) {
    knots = Knot[rrtg == .rrtg]
    credential(knots, Knot_quipu, ppmwid = ppmwid, rtwid = rtwid, factor = factor, mpc = mpc, ratio = Ratio, ratio.lim = Ratio.lim, maxdimer = maxdimer, cd = cd, do.plot = do.plot)
  }
  cat("\rWorking with supported charge states.", "Found", length(unique(Knot_quipu$q)), "credentialed knots.")
  lastn = length(unique(Knot_quipu$q))
  
  for (.z in unique(c(Knot$z, .zs))) {
    cat("\nWorking with unsupported charge state:", .z)
    lastn = length(unique(Knot_quipu$q))
    if (.z == 0) next
    
    subknots = Knot[Knot_quipu[is.na(q)], , on="knot"][,gtemp := as.integer(clustgroup(cbind(meanr, rt), scales, factor))][z %in% c(0, .z)]
    for (.gtemp in unique(subknots$gtemp)) {
      knots = subknots[gtemp == .gtemp]
      knots[,meanr := meanmz %/% (cd/.z)]
      knots[1,z:=.z]
      credential(knots, Knot_quipu, ppmwid = ppmwid, rtwid = rtwid, factor = factor, mpc = mpc, ratio = Ratio, ratio.lim = Ratio.lim, maxdimer = maxdimer, cd = cd, do.plot = do.plot)
    }
    
    cat("\rWorking with unsupported charge state:", .z, "Found", length(unique(Knot_quipu$q)) - lastn, "credentialed knots.")
  }
  
  
  Knot_quipu[,':='(quipu = as.integer(factor(q)), q = NULL)]
  
  Quipu = Knot_quipu[Knot,,on="knot"][,.(minsupport = min(n), maxsupport = max(n), nknot = .N, ratio = calcratio(.SD)), by="quipu"][!is.na(quipu) & !duplicated(quipu)][ratio > Ratio*Ratio.lim & ratio < Ratio/Ratio.lim]
  
  cat("\nAfter isotope ratio check, found", nrow(Quipu), "credentialed knots (quipus).")
  
  Quipu_stat = Quipu[,.(quipu,minsupport,maxsupport,nknot)]
  
  Quipu = Knot[Knot_quipu[Quipu,,on="quipu"],,on="knot"]
  Quipu = Quipu[,.(nknot = nknot, npeak = sum(n), rtmean = mean(rt), charge=z, basemz = min(basemz), mainmz1 = min(mainmz), mainmz2 = max(mainmz), int1 = min(maxi), int2 = max(maxi), ratio = ratio), by = "quipu"]
  Quipu = Quipu[!duplicated(Quipu$quipu)]
  
  return(list(knot_quipu = Knot_quipu[which(quipu %in% unique(Quipu$quipu))], quipu = Quipu, quipu_stat = Quipu_stat))
}

#' Core algorithm of credentialknots: Heuristic search for knots that satisfy credentialing filters
#'
#' @description This functions identify credentialed knots pairs by intensity ratio check, mass-per-carbon check,
#' direction check followed by graph-clustering. 
#' @usage credential(knots, Knot_quipu, ppmwid, rtwid, factor, mpc, ratio, ratio.lim, maxdimer, cd, do.plot)
#' @param knots data.table. A subset of knots with same charge state to be credentialed.
#' @param Knot_quipu data.table The initial index of knot-to-quipu assignment.
#' @param ppmwid numeric Mass error tolerance (ppm) of credentialed knots.
#' @param rtwid numeric Retention time tolerance (seconds) of credentialed knots.
#' @param factor numeric Height cutoff of hierarchical clustering knots by (meanr, rt)
#' @param mpc numeric Mass-per-carbon filter. A vector of two values setting minimum and maximum theorectical mass per carbon.
#' @param ratio numeric The mixing ratio of unlabled/labeled (12C/13C) samples.
#' @param ratio.lim numeric Ratio filter. A decimal number (0,1] controling ranges of acceptable deviation of 
#' unlabeled/labeled ratio from mixing ratio.
#' @param maxdimer numeric Maximumly allowed number of knots in each credentialed groups, The default value is 4.
#' @param cd numeric Unit mass difference between unlabeled and labeled atoms. Defalut is 13C-12C = 13.00335 - 12.
#' @param do.plot logical If set to TRUE, plot the credentialing scores of the knots
#' @keywords credential credentialknots
#' @import igraph magrittr stats matrixStats
#' @seealso \code{\link{credentialknots}} \code{\link{findknots}}
#' @return \code{Knot_quipu} data.table The new index of knot-to-quipu assignment after credentialing. 
#'

credential = function(knots, Knot_quipu, ppmwid, rtwid, factor, mpc, ratio, ratio.lim, maxdimer, cd, do.plot) {
  
  scales = c(ppmwid * 700 / 1E6, rtwid)
  knots[z==0, z:=max(knots$z)]
  
  # Split up peaks if possible
  withinmpc = { (-1 * knots$mainmz * knots$z / cd) / round(outer(knots$mainmz, knots$mainmz, "-")) } %>% { . > mpc[1] } #& . < mpc[2]  } #Dont use max mpc at this stage - will disqualify dimers.
  
  knots[,direction := dir %>% { .[abs(.) < 0.1] = 0; . } %>% sign]
  knots[is.na(direction),direction:=0]
  dirworks = outer(knots$direction, knots$direction, "<") | outer(knots$direction, knots$direction) == 0
  
  intsanity = outer(knots$maxi, knots$maxi, "/") %>% { . < ratio/ratio.lim & . > ratio*ratio.lim}
  
  poss = which(withinmpc & dirworks & intsanity, arr.ind = T)
  
  if (nrow(poss) < 1) return(NULL)
  
  gs = graph.data.frame(poss) %>% clusters %>% '[['("membership")
  gs = gs[order(as.numeric(names(gs)))]
  
  knots[as.numeric(names(gs)),mdig := gs]
  
  for (.mdig in unique(na.omit(knots$mdig))) {
    knots2 = knots[mdig == .mdig]
    
    # Generate and assess combinations of peaks
    maxdimerx = maxdimer; maxdimerx[maxdimerx > nrow(knots2)] = nrow(knots2)
    cs = do.call(what=cbind.fill, lapply(seq_len(maxdimerx)[-1], combn, x = seq_len(nrow(knots2))))
    
    #Calculate intenisty residual for each combination
    gscore = vector(length=ncol(cs), mode="numeric")
    for (j in seq_len(ncol(cs))) {
      knots3 = knots2[cs[,j]][!is.na(maxi)]
      
      npeaks = sum(!is.na(knots3$maxi))
      setkey(knots3, "meanmz")
      
      
      spacing = diff(round((knots3$mainmz * knots3$z) %/% 1))
      if (npeaks == 2) {
        ints = convolve(c(1, ratio))
        
      } else if (npeaks == 3) {
        if (length(unique(spacing)) == 1) {
          ints = convolve(c(1, ratio), c(1, ratio))
        } else { # Cant have three peaks and unequal spacing
          ints = c(0,0,0)
        }
        
      } else if (npeaks == 4) {
        if (length(unique(spacing)) == 1) { # Unequal spacing necessitates dimer type pattern
          ints = convolve(c(1, ratio), c(1, 0, ratio))
        } else { #Equal spacing necessitates trimer type pattern
          ints = convolve(c(1, ratio), c(1, ratio), c(1,ratio))
        }
        
      } else { # Assume Homomultimer
        ints = do.call(what=convolve, lapply(seq_len(nrow(knots3)-1), function(x) c(1, ratio)))
      }
      
      gscore[j] = sum((knots3$maxi/max(knots3$maxi, na.rm=T) - ints)^2)
    }
    
    
    # Calculate scores based on rt and meanr
    distances = dist(knots2[,.(meanr, rt)]/rep(scales, each = nrow(knots2))) %>% as.matrix
    
    dscore = vector(length=ncol(cs))
    for (i in seq_len(ncol(cs))) {
      measureme = sapply(seq_len(length(na.omit(cs[,i]))-1), function(j) cs[j:(j+1),i]) %>% aperm
      dscore[i] = sapply(seq_len(nrow(measureme)), function(j) distances[measureme[j,,drop=F]]) %>% sum
    }
    
    
    #Calculate feasibility based on mass per carbon
    mzs = is = zs = ns = cs
    mzs[] = knots2$mainmz[cs]
    zs[] = knots2$z[cs]
    is[] =  knots2$maxi[cs]
    ns = colSums(!is.na(cs))
    cnums = round(matrixStats::rowDiffs(colRanges(mzs,na.rm=T))/cd * matrixStats::colMaxs(zs, na.rm=T))
    
    mpc2 = colMaxs(zs, na.rm=T) * colMins(mzs,na.rm=T) / colRanges(mzs, na.rm=T) %>% rowDiffs %>% { (. * colMaxs(zs, na.rm=T)) %/% cd }
    mpctf = c(mpc2 > mpc[1] & mpc2 < mpc[2])
    
    dirtf = vector(length=ncol(cs))
    for (i in seq_len(ncol(cs))) {
      dirs = knots2$direction[na.omit(cs[,i])]
      dirs[is.na(dirs)] = 0
      r = sum(dirs == 1)
      n = sum(dirs == 0)
      f = sum(dirs == -1)
      
      
      dirtf[i] = f+r <= 3 & tail(dirs[order(knots2$mainmz[na.omit(cs[,i])])], n = 1) != -1 & head(dirs[order(knots2$mainmz[na.omit(cs[,i])])], n = 1) != 1
    }
    
    if (sum(mpctf & dirtf) < 1) return(NULL)
    
    #Use this information to make quipu
    #gscore (ints), dscore (mz, rt), cnums (carbon number), mpctf, direction
    cs = cs[,mpctf & dirtf,drop = F]
    scores = rbind(gscore, dscore, cnums = c(cnums))[,mpctf&dirtf,drop=F]
    
    if (do.plot) plot(scores[1:2,,drop=F] %>% aperm, xlim = c(0, 4), ylim = c(0,4))
    
    
    totscore = c(colSums(scores[1:2,,drop=F]^2)^0.5)
    
    tso = order(totscore)
    creds = tso[1]
    fails = 0
    repeat {
      if (tso[length(tso)] == creds[length(creds)]) break
      
      newi = tso[length(creds)+1]
      
      if (!(any(na.omit(cs[,newi]) %in% na.omit(cs[,creds])))) {
        creds = c(creds, newi)
      } else {
        fails = fails + 1
      }
      
      if (fails > 3) break
    }
    
    
    DT = cs[,creds] %>% melt %>% data.table
    Knot_quipu[knots2[DT$value],q := paste(paste(knots$knot, collapse = " "), .mdig, DT$Var2),on="knot"]
  }
}
#' Calculation of ratio of unlabeled/labeled (12C/13C) peak intensities.
#' 
#' @usage calcratio(knot)
#' @param knot a set of knots to be calculated. The data structure are the same as
#' @return numeric Ratio of measured intensities 
#' @seealso \code{\link{findknots}} \code{\link{findknots}}
#' @import data.table

calcratio = function(knot) {
  li = which.min(knot$mainmz)
  hi = which.max(knot$mainmz)
  knot$maxi[li]/knot$maxi[hi]
}

