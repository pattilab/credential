#' Step 4 of Credentialing: Find Matched Quipus between Two Credentialing Groups
#'
#' @description This function matches two set of quipus (credentialed knots) by their retention time cutoff, 
#' head&tail m/z, charge state, peak intensities ratio. 
#' @usage credentialquipu(credentialedknots1, credentialedknots2, ppm = 15, rtwin = 1, ratio_ratio = 2, 
#' ratio_ratio_tol = 0.5, tailmatch=T)
#' @param credentialedknots1 list The original output of credentialknots() for the first credentialing group.
#' @param credentialedknots2 list The original output of credentialknots() for the second credentialing group.
#' @param ppm numeric Maximum mass error tolerance when matching head and tail peaks of two quipus.
#' @param rtwin numeric Maximum mean retention time difference (rtmean) when matching two quipus.
#' @param ratio_ratio numeric Ratio of unlabeled/labeled mixing ratios of first and second credentialed samples. The default is 2.
#' @param ratio_ratio_tol numeric Ratio filter. A decimal number (0,1] controling ranges of acceptable deviation of 
#' ratio1/ratio2 from ratio_ratio, The default value is 0.5.
#' @param tailmatch logical Whether or not to include tail peak m/z when matching quipus.
#' @import data.table magrittr
#' @keywords credential credentialquipu
#' @return list A list of three tables "credentialedindex" and credentialedgroups". \code{credentialedindex} data.table
#' The index of quipu-to-quipu assignment between two credentialing conditions. \code{credentialedgroup} data.table
#' The final credentialed groups according to assignment of "credentialedindex", including 'quipu',
#' 'nknot', 'npeak', 'rtmean', 'charge','basemz', 'mainmz1', 'mainmz2', 'int1','int2', 'ratio' and
#' 'ratio_ratio'.
#' @export

credentialquipu <- function(credentialedknots1, credentialedknots2, ppm = 15, rtwin = 1, ratio_ratio = 2, ratio_ratio_tol = 0.5, tailmatch=T){

  quipu1 <- as.data.table(credentialedknots1$quipu)
  quipu2 <- as.data.table(credentialedknots2$quipu)
  
  colnames(quipu1) <- paste0(colnames(quipu1),"1",sep="")
  colnames(quipu2) <- paste0(colnames(quipu2),"2",sep="")
  
  quipu1 <- quipu1[order(quipu1)]
  quipu2 <- quipu2[order(quipu2)]
  
  # retention time match
  rtdiff = outer(quipu1$rtmean1, quipu2$rtmean2, "-")
  rt_ind = which(abs(rtdiff)<=rtwin, arr.ind=T)
  
  # charge state match
  charge = outer(quipu1$charge1, quipu2$charge2, "==")
  cs_ind = which(charge==TRUE, arr.ind = T)
  
  # mainmz match
  mzppm1 = outer(quipu1$mainmz11,quipu2$mainmz12, function(x,y){ifelse(x/(1+ppm*1e-6)<=y/(1-ppm*1e-6) & y/(1+ppm*1e-6)<=x/(1-ppm*1e-6),1,0)})
  mzppm2 = outer(quipu1$mainmz21,quipu2$mainmz22, function(x,y){ifelse(x/(1+ppm*1e-6)<=y/(1-ppm*1e-6) & y/(1+ppm*1e-6)<=x/(1-ppm*1e-6),1,0)})
  
  # head&tail or head only
  if (tailmatch) {
  mz_ind = which(mzppm1 == 1 & mzppm2 ==1, arr.ind=T)
  } else {mz_ind = which(mzppm1 == 1, arr.ind=T)}
  
  match_ind = rbind(mz_ind, rt_ind)[duplicated(rbind(mz_ind,rt_ind)),]
  match_ind = rbind(match_ind, cs_ind)[duplicated(rbind(match_ind,cs_ind)),]
  
  # isolate duplicated match

  if(sum(duplicated(match_ind[,1]))>0) {
    tmp.ind = which(duplicated(match_ind[,1]) & duplicated(match_ind[,1], fromLast = T))
    dup_ind1 = match_ind[tmp.ind,]
    match_ind = match_ind[setdiff(1:nrow(match_ind), tmp.ind),]
  }
  if(sum(duplicated(match_ind[,2]))>0) {
    tmp.ind = which(duplicated(match_ind[,2]) & duplicated(match_ind[,2], fromLast = T))
    dup_ind2 = match_ind[tmp.ind,]
    match_ind = match_ind[setdiff(1:nrow(match_ind), tmp.ind),]
  }
  
  # rank duplicate match, select the closest match
  
  # working on it...
  quipu_match = cbind(quipu1[match_ind[,1]],quipu2[match_ind[,2]])
  # intensity ratio check
  quipu_match[,"ratio1_ratio2" := ratio1/ratio2][,"credentialed":= ratio1_ratio2 <= ratio_ratio / ratio_ratio_tol & ratio1_ratio2 >= ratio_ratio * ratio_ratio_tol]
  
  # clean up data
  quipu_match = quipu_match[order(basemz1)]
  credentialedquipu = list(credentialedindex = quipu_match[credentialed == TRUE,.(quipu1,quipu2)], credentialedgroups = quipu_match[credentialed == TRUE])
  
  return(credentialedquipu)
}

# Under construction
rankmatch = function(dup_ind, quipu_match, ratio_ratio = 2, weight = c(0.6,0.3,0.1,0.1)) {
  
}



  