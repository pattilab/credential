#' Credentialing Untargeted Metabolomics Dataset
#'
#' The dispatching function of credentialing
#' @usage credentialing(peaktable1, peaktable2, ppm=15, rtwin=1, rtcom=3, ratio1= 1/1, ratio2 = 1/2, 
#' ratio_tol=0.1, ratio_ratio_tol = 0.2, cd=13.00335-12, charges= 1:4, mpc=c(12,120), maxnmer=4, 
#' label1 = "1T1", label2 = "1T2",
#' export = T, plot=F, projectName = "Credentialing",...)
#' @description Credentialing is a method to authenticate biologically relevant features in untargeted metabolomics 
#' data. The analysis is based on two peak (feature) tables with index (cc), m/z (mz), retention time (rt) and 
#' intensity (i) values. Each peak table represents a credentialing group where unlabeled and labeled biological 
#' samples are mixed with certain ratio (ratio1, ratio2). First, credentialing searches 'knot' at designated 
#' charge states in each peak table. A knot is a set of features with defined retention time and m/z spacing. 
#' Second, credentialing breaks merged knots with multiple isotopologue patterns. Third, credentialing 
#' searches 'quipu' from knots. A quipu is a set of knots with similar retention time and m/z spacing, charge state, 
#' isotopologue pattern, qualifying mass per carbon and peak intensity ratios. Fourth, credentialing matches 
#' credentialed quipus between two credentialing groups by retention time, charge state, head/tail m/z peaks and 
#' ratio of their intensities.
#' @param peaktable1 data.table Peak table of the first credentialing group (ratio1).
#' @param peaktable2 data.table Peak table of the second credentialing group (ratio2).
#' @param ppm numeric Mass error tolerance (ppm) of credentialing peaks.
#' @param rtwin numeric Retention time tolerance (seconds) for credentialing peaks within group.
#' @param rtcom numeric Retention time tolerance (seconds) for credentialing peaks between groups.
#' @param ratio1 numeric Ratio of 12C/13C (unlabeled/labeled samples) mixing ratio of first credentialing group (peaktable1), default as 1/1.
#' @param ratio2 numeric Ratio of 12C/13C (unlabeled/labeled samples) mixing ratio of second credentialing group (peaktable2), default as 1/2.
#' @param ratio_tol numeric A decimal number (0,1] controling range of acceptable intensity ratios as credentialing peaks within group. The default value is 0.1.
#' @param ratio_ratio_tol numeric A decimal number (0,1] controling range of acceptable intensity ratios as credentialing peaks between groups. The default value is 0.1.
#' @param cd numeric Unit mass difference between unlabeled and labeled atoms. Defalut is 13C-12C = 13.00335 - 12.
#' @param charges numeric Possible charge states of isotopologues. Default value is 1:4.
#' @param mpc numeric A vector of two values setting minimal and maximal theoretical mass per carbon to be considered. The default value is c(12,120).
#' @param maxnmer numeric Maximum possible number of knots accepted in each credentialed groups, The default value is 4.
#' @param label1 character plotting legend of ratio1.
#' @param label2 character plotting legend of ratio2. 
#' @param export logical If TRUE, $credentialedgroups and $credentialedpeaks will be exported in CSV files.
#' @param plot logical If TRUE, a PDF with credentialed peak groups' plots will be generated.
#' @param projectName character Title of the credentialing project.
#' @param ... Further arguments to be passed.
#' @keywords credential credentialing metabolomics
#' @import data.table magrittr utils ggplot2
#' @return list credentialed features. \code{credentialedgroups} data.table A summary of credentialed quipus with
#' following information: quipu#, nknot (number of knots), npeak (number of peaks), rtmean (mean retention time), 
#' basemz (lowest mz), mainmz 1 (credentialed pair - unlabeled mz). mainmz 2 (credentialed pair - labeled mz), 
#' int 1 (credentialed pair - intensity of unlabeled mz). int 2 (credentialed pair - intensity of labeled mz), 
#' ratio (ratio of unlabled/labeled intensity), credentialed (TRUE/FALSE). \code{credentialedpeaks} data.table A summary
#' of all credentialed peaks in two credentialing groups with following information: quipu#, knot#, cc# (feature index),
#' charge (charge state), mz (m/z), rt (retention time), i (intensity), tail (if this peak is tail peak in knot),
#' mainmz 1, mainmz 2, ratio. \code{credentialedknots1} list A list of credentailed knots (quipu) in first 
#' credentialing group. \code{credentialedknots2} list A list of credentailed knots (quipu) in second credentialing group.
#' \code{knots1} list All knots in first credentialing group (peaktable2) \code{knots2} list All knots in second credentialing group 
#' (peaktable2). \code{CredentialParams} list Parameters used in this analysis.
#' @seealso \code{\link{findknots}} \code{\link{credentialknots}} \code{\link{credentialquipu}}
#' @export

credentialing = function(peaktable1, peaktable2, ppm=15, rtwin=1, rtcom=3, 
                         ratio1= 1/1, ratio2 = 1/2, ratio_tol=0.1, ratio_ratio_tol = 0.2, 
                         cd=13.00335-12, charges= 1:4, mpc=c(12,120), maxnmer=4,
                         label1 = "1T1", label2 = "1T2",
                         export = T, plot=F, projectName = "Credentialing",
                         ...){

  #initiation
  peaktable1 = data.table(peaktable1)
  peaktable2 = data.table(peaktable2)
  # 1st round credentialing

  cat("\nCredentialing peaktable1...\n")

  # find possible isotope head and tails at user-defined charge states
  knots1 = findknots(peaktable1, .zs = charges, ppmwid = ppm, rtwid = rtwin, cd = cd)
  # Heuristic search for knots that satisfying credentialing filters
  credentialedknots1 = credentialknots(knots1, ppmwid = ppm, rtwid = rtwin, mpc = mpc, Ratio = ratio1, Ratio.lim = ratio_tol, maxnmer = maxnmer, cd = cd, .zs = charges)

  cat("\nCredentialing peaktable2...\n")
  
  #credentialing with ratio2 combination
  knots2 = findknots(peaktable2, .zs = charges, ppmwid = ppm, rtwid = rtwin, cd = cd)
  credentialedknots2 = credentialknots(knots2, ppmwid = ppm, rtwid = rtwin, mpc = mpc, Ratio = ratio2, Ratio.lim = ratio_tol, maxnmer = maxnmer, cd = cd, .zs = charges)
  
  # 2nd round filtering

  cat("\nCredentialing all quipus...")
  credentialedquipus <- credentialquipu(credentialedknots1, credentialedknots2, ppm = ppm, rtwin = rtcom, ratio_ratio = ratio1/ratio2, ratio_ratio_tol = ratio_ratio_tol, tailmatch = T)
  
  # data clean-up
  
  ft_1 = peaktable1[knots1$cc_knot[credentialedknots1$knot_quipu[!is.na(quipu)],,on="knot"],,on="cc"][credentialedquipus$credentialedgroups[,.(quipu=quipu1,charge1,mainmz11,mainmz21,ncar1,ratio1)],,on="quipu"]
  ft_1 = ft_1[,.(cc1=cc,mz1=mz,rt1=rt,i1=i,knot1=knot,tail1=tail,quipu1=quipu,charge1,mainmz11,mainmz21,ncar1,ratio1)]
  ft_2 = peaktable2[knots2$cc_knot[credentialedknots2$knot_quipu[!is.na(quipu)],,on="knot"],,on="cc"][credentialedquipus$credentialedgroups[,.(quipu=quipu2,charge2,mainmz12,mainmz22,ncar2,ratio2,ratio1_ratio2)],,on="quipu"]
  ft_2 = ft_2[,.(cc2=cc,mz2=mz,rt2=rt,i2=i,knot2=knot,tail2=tail,quipu2=quipu,charge2,mainmz12,mainmz22,ncar2,ratio2,ratio1_ratio2)]
  
  credentialedpeaks = do.call(rbind,apply(credentialedquipus$credentialedindex, MARGIN = 1, function(x){cbind.fill(ft_1[quipu1==x[1]][order(mz1)],ft_2[quipu2==x[2]][order(mz2)])}))
  
  cat("\n Credentialing finished...\nFound", nrow(credentialedquipus$credentialedgroups), "credentialed peaks.")
  
  CredentialParams <- list(ppm=ppm, rtwin=rtwin, rtcom=rtcom, 
                           ratio1 = ratio1, ratio2 = ratio2, ratio_tol = ratio_tol, ratio_ratio_tol = ratio_ratio_tol, 
                           cd = cd, charges = charges, mpc = mpc, maxnmer = maxnmer, label1 = label1, label2 = label2,
                           export = export, plot = plot, projectName = projectName)
  
  credentialed = list(credentialedgroups = credentialedquipus$credentialedgroups,
                      credentialedpeaks = credentialedpeaks,
                      credentialedindex = credentialedquipus$credentialedindex,
                      credentialedknots1 = credentialedknots1,
                      credentialedknots2 = credentialedknots2,
                      knots1 = knots1, knots2 = knots2,
                      CredentialParams = CredentialParams)
  
  if(plot) {
    plotcredpeaks(Credentialedindex = credentialedquipus$credentialedindex, Credentialedpeaks = credentialedpeaks, cred1 = label1, cred2 = label2,
                     filename = paste(paste(projectName,nrow(credentialedquipus$credentialedindex),"credentialed_peak_groups",sep="_"),".pdf",sep=""))
  }
  
  if(export){
    write.csv(credentialed$credentialedgroups,file = paste0(projectName,"_CredentialedGroups.csv"))
    write.csv(credentialed$credentialedpeaks,file = paste0(projectName,"_CredentialedPeaks.csv"))
    cat("\nCredentialing results are exported under:\n", getwd())
  }
  return(credentialed)
}
