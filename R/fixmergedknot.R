#' Step 2 of Credentialing: Split Merged knots with Multiple Isotopologue Patterns
#'
#' @description This function breaks knots with multiple isotopologue patterns for credentialing step 3. If certain
#' descending and ascending pattern are found at the beginning (head) and ending (tail) region of the konts, the knot
#' are broken half as two new knots.
#' @usage fixmergedknots(Knot, features) 
#' @param Knot list The output of findknots() function.
#' @param features data.table Feature table corresponding to the Knot.
#' @import data.table utils 
#' @return list A list with identical data sturcture as the input \code{Knot}.
#' @keywords credentialing, fixmergedknot, metabolomics
#' @seealso \code{\link{findknots}} \code{\link{credentialing}}
#' @export

fixmergedknots = function(Knot, features) {
  
  merged = Knot$cc_knot[features,,on="cc", nomatch=0][,merged := checkmerged(.SD), by="knot"][,.(cc,merged)]
  Knot$cc_knot[merged,merged := merged,on="cc"]
  
  temp = Knot$knot[Knot$cc_knot,,on="knot"]
  
  temp[,':='(knot=as.integer(factor(paste(knot, merged))), merged = NULL)]
  
  knot_new = features[,.(mz, i, rt, cc)][temp,,on="cc", nomatch = 0][,knotstats(.SD),by="knot"]
  
  Knot$knot = knot_new
  
  Knot$cc_knot = temp[,.(cc, knot, tail)]
  
  return(Knot)
}

#' Search and label if certain increase and decrease are found at the head and tail region respectively.
#' @param knot data.table A single knot with all aggregated information
#' @seealso \code{\link{findknots}}
#' @return numeric a vector of values (1 or 2) indicating whether features in the knot should be splited or not.
#' @usage checkmerged(knot)

checkmerged = function(knot) {
  if (sum((!knot$tail))<4) { return(rep(1, nrow(knot))) }
  
  cc.o = knot$cc
  o = order(knot$mz) 
  
  half = floor(sum(!knot$tail)/2)
  full = sum(!knot$tail)
  
  front = knot[o][tail == F][1:half]
  back = knot[o][tail == F][((half+1):full)]
  
  s1 = mean(diff(head(front$i,3)))/max(head(front$i,3))
  s2 = mean(diff(tail(back$i,3)))/max(tail(back$i,3))
  
  if (s1 < -0.2 & s2 > 0.2) {
    return(c(rep(1, half), rep(2, nrow(knot)))[match(knot$cc, cc.o[o])])
  }
  else return(rep(1, nrow(knot)))
}

#' Calculate the direction of each knot
#' @usage calcdir(ps)
#' @param ps data.table A signle knot with aggregated information.
#' @return numeric A value indicating the trend of the knot. Ascending intensities have positive value, 
#' descending intensities have negeative value. 
#' @seealso \code{\link{findknots}}
#' @import data.table utils

calcdir = function(ps) {
  ps = subset(ps, tail == F)
  
  f = tail(seq_len(nrow(ps)), n=3)
  b = head(seq_len(nrow(ps)), n=3)
  
  ds = c(
    mean(diff(ps$i[order(ps$mz)][f]))/max(ps$i),
    mean(diff(ps$i[order(ps$mz)][b]))/max(ps$i)
  )
  
  val = ds[which.max(abs(ds))]
  if (length(val) < 1) return(as.numeric(NA))
  val
}

#' Re-aggregate information of each knot after splitting
#' @usage knotstats(knot)
#' @param knot data.table A signle knot with aggregated information.
#' @return data.table Aggregated information of new knots.
#' @seealso \code{\link{findknots}}
#' @import data.table

knotstats = function(knot) {
  
  o = order(knot$mz)
  
  stat = data.table(
    meanr = mean(knot$mz %% (13.00335-12)),
    meanmz = mean(knot$mz),
    basemz = min(knot$mz),
    mainmz = knot$mz[which.max(knot$i)],
    rt = mean(knot$rt),
    maxi = max(knot$i),
    n = nrow(knot),
    z = knot$z[1],
    dir = calcdir(knot[,.(mz, i, tail)])
  )
  return(stat)
}