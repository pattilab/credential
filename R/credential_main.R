#' #' Credentialing main function
#'
#' The main function of credentialing.
#'
#' @usage credentialing(peaktable1, peaktable2, ppm, rtwin, rtcom, ratio1= 1/1, ratio2 = 1/2, ratio_tol=0.1, ratios_tol = 0.2, 
#' cd=13.00335-12, charges= 1:4, mpc=c(12,120), maxnmer=4)
#' @param peaktable1 data.table Feature table corresponding to the first ratio conditions.
#' @param peaktable2 data.table Feature table corresponding to the second ratio conditions.
#' @param ppm numeric Mass error tolerance for knot searching and credentialing.
#' @param rtwin numeric Retention time window for 1st round credentialing
#' @param rtcom numeric Retention time window for 2nd round credentialing
#' @param ratio1 numeric Ratio of designated 12C/13C ratio in peaktable1, default set to 1/1
#' @param ratio2 numeric Ratio of designated 12C/13C ratio in the first peaktable, by default set to 1/2
#' @param ratio_tol numeric A decimal that controls the ratio range to pass the 1st round credentialing, The default value is 0.1.
#' @param ratios_tol numeric A decimal that controls the ratio range to pass the 2nd round credentialing, The default value is 0.2.
#' @param cd numeric Unit mass difference among isotopologues. Defalut is 13C-12C = 13.00335 - 12.
#' @param charges numeric Possible charge states to be considered when searching isotope pairs. Default value is 1:4.
#' @param mpc numeric A vector of two values setting minimal and maximal theoretical mass per carbon to be considered. The default value is c(12,120)
#' @param maxnmer numeric The maximal number of knots to accept within each credentialed features, The default value is 4.
#' @description Credentialing workflow consists of four major parts. The analysis starts from two feature tables including index, m/z, retention time and intensity values for each isotope conditions.
#' First, credentialing searches for 'knots' at all possible charge states in each feature table. A knot is defined as a set of features with defined retention time window and m/z ppm error.
#' Second, credentialing breaks the knots with multiple isotopologue patterns. Third, credentialing clusteres searches for credentialed knots 
#' based on the similarity, charge state, isotopologue pattern, mass per carbon violation, and peak intensity ratios. Fourth, credentialing matches credentialed 'quipu' from both isotope mixing conditions
#' by checking head and tail peak m/z, retention time and intensity. 1st round credentialing is from step 1 to step 3 and 2nd round credentialing is step 4. 
#' The output of credentialing is a list consisting of credentialed features and credentialed knots at different levels.
#' @keywords credentialing, untargeted metabolomics
#' @import data.table
#' @return A list of credentialed features.
#' @export

credentialing = function(peaktable1, peaktable2, ppm, rtwin, rtcom, ratio1= 1/1, ratio2 = 1/2, ratio_tol=0.1, ratios_tol = 0.2, cd=13.00335-12, charges= 1:4, mpc=c(12,120), maxnmer=4){

  #initiation
  peaktable1 = data.table(peaktable1)
  peaktable2 = data.table(peaktable2)
  # 1st round credentialing

  cat("Start 1st round Credentialing on peaktable1...\n")

  # find possibel isotope head and tails in different charge states
  knots_f1 = findknots(peaktable1, .zs = charges, ppmwid = ppm, rtwid = rtwin, cd = cd)
  # Resolving issues with merged knots
  knots_f1 = fixmergedquipu(knots_f1,peaktable1)
  # heuristic search for knots that satisfying credentialing filters
  credentials_f1 = credentialknots(knots_f1$knot, ppmwid = ppm, rtwid = rtwin, mpc = mpc, Ratio = ratio1, Ratio.lim = ratio_tol, maxnmer = maxnmer, cd = cd)
  # labelled credentialed features at original peaktable
  ft_1 = peaktable1[knots_f1$cc_knot[credentials_f1$knot_quipu[!is.na(quipu)],on="knot"],,on="cc"]

  cat("Start 1st round Credentialing on peaktable2...\n")

  #credentialing with ratio2 combination
  knots_f2 = findknots(peaktable2, .zs = charges, ppmwid = ppm, rtwid = rtwin, cd = cd)
  knots_f2 = fixmergedquipu(knots_f2,peaktable2)
  credentials_f2 = credentialknots(knots_f2$knot, ppmwid = ppm, rtwid = rtwin, mpc = mpc, Ratio = ratio2, Ratio.lim = ratio_tol, maxnmer = maxnmer, cd = cd)

  ft_2 = peaktable2[knots_f2$cc_knot[credentials_f2$knot_quipu[!is.na(quipu)],on="knot"],,on="cc"]

  # rearrangement of credentialed features
  dt1 = ft_1[credentials_f1$quipu[,c("quipu","ratio")][!is.na(quipu)],,on="quipu"]
  dt2= ft_2[credentials_f2$quipu[,c("quipu","ratio")][!is.na(quipu)],,on="quipu"]

  dtR1 = dt1[order(dt1[,"quipu"], dt1[,"mz"]),]
  dtR2 = dt2[order(dt2[,"quipu"], dt2[,"mz"]),]

  # 2nd round filtering
  
  cat("Start 2nd round Credentialing...\n")
  match_cf = matchcredfeature(dtR1,dtR2,ppm=ppm,drt=rtcom,ratio=ratio1/ratio2,ratio_tol = ratios_tol)

   # data output

  credentialing = list(CredentialedFeature1R1=ft_1, CredentialedFeature2R1=ft_2, knots1=knots_f1, knots2=knots_f2, credentialedKnots1=credentials_f1,credentialedKnots2=credentials_f2,
                       CredentialedFeatureGroups = match_cf$Credentialed_FeatureGroups, CredentialedFeatureR2=match_cf$Credentialed_Features, CredentialedFeatureR2F = match_cf$Credentialed_Features_Filtered,
                       CredentialedFeature1N2 = match_cf$NomatchFeatures_Group1, CredentialedFeature2N2 = match_cf$NomatchFeatures_Group2)

  return(credentialing)
}
