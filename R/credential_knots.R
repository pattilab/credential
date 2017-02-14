credential = function(knots, Knot_quipu, ppmwid, rtwid, factor, mpc, ratio, ratio.lim, maxdimer, cd, do.plot, scales) {
  #cat("\rWorking on knots:", knots$knot, "                                                      ")
  
  
  knots[z==0, z:=max(knots$z)]
  
  # Split up peaks if possible
  withinmpc = { (-1 * knots$mainmz * knots$z / cd) / round(outer(knots$mainmz, knots$mainmz, "-")) } %>% { . > mpc[1] } #& . < mpc[2]  } #Dont use max mpc at this stage - will disqualify dimers.
  
  knots[,direction := dir %>% { .[abs(.) < 0.1] = 0; . } %>% sign]
  knots[is.na(direction),direction:=0]
  dirworks = outer(knots$direction, knots$direction, "<") | outer(knots$direction, knots$direction) == 0
  
  intsanity = outer(knots$maxi, knots$maxi, "/") %>% { . < 1/ratio.lim & . > ratio.lim}
  
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

#' Heuristic search for knots which satisfy credentialing filters
#'
#' @param Knots data.frame Rows are knots.
#' @param ppmwid numeric The maximum mass error in ppm.
#' @param rtwid numeric The maximum retention time difference in seconds.
#' @param mpc numeric Vector containing the minimum and maximum allowable mass per carbon.
#' @param ratio numeric The ratio of U12C/U13C.
#' @param ratio.lim numeric The maximum deviation from that ratio.
#' @param maxnmer integer The biggest n-mer to search for.
#' @param cd numeric The mass spacing to search for (defaults to C13 - C12)
#' 
#' @seealso \link{\code{findknots}}
#' 
#' @return list A list with values "knot_quipu" and "quipu".  \code{knot_quipu} assigns feature knots to quipu - credentialed groups. \code{quipu} contains aggregate information about each credentialed group.
#'
credentialknots = function(Knots, ppmwid = 9, rtwid = 1, mpc = c(12, 120), ratio = 1/1, ratio.lim = 0.1, maxnmer = 4, cd = 13.00335-12, .zs = 1:4) {
  cat("\nCredentialing within", length(unique(Knots$knot)), "supplied knots.")
  factor = 1
  maxdimer = maxnmer
  do.plot = F
  Knot = copy(Knots)
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
    credential(knots, Knot_quipu, ppmwid = ppmwid, rtwid = rtwid, factor = factor, mpc = mpc, ratio = ratio, ratio.lim = ratio.lim, maxdimer = maxdimer, cd = cd, do.plot = do.plot, scales = scales)
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
      credential(knots, Knot_quipu, ppmwid = ppmwid, rtwid = rtwid, factor = factor, mpc = mpc, ratio = ratio, ratio.lim = ratio.lim, maxdimer = maxdimer, cd = cd, do.plot = do.plot,scales = scales)
    }
    
    cat("\rWorking with unsupported charge state:", .z, "Found", length(unique(Knot_quipu$q)) - lastn, "credentialed knots.")
  }
  
  
  Knot_quipu[,':='(quipu = as.integer(factor(q)), q = NULL)]
  
  Quipu = Knot_quipu[Knot,,on="knot"][,.(minsupport = min(n), maxsupport = max(n), nknot = .N, ratio = calcratio(.SD)), by="quipu"][!is.na(quipu) & !duplicated(quipu)]
  
  cat("\nFound", nrow(Quipu)-1, "credentialed knots.")
  
  list(knot_quipu = Knot_quipu, quipu = Quipu)
}


calcratio = function(knot) {
  li = which.min(knot$mainmz)
  hi = which.max(knot$mainmz)
  
  knot$maxi[li]/knot$maxi[hi]
}

