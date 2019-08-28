# installation

install.packages("devtools")
devtools::install_github("lingjuewang/credential", ref="release/3.1")

# load required library
library(credential)
library(xcms)
library(data.table)

# extract features from xcmsSet

  # please manually download sample xcmsSet data from followiing url
  # https://github.com/lingjuewang/xcmsdata/blob/master/credentialxcms.Rdata
  # xcms can run on credentialing conditions separately or together ("1T1", "1T2")
  # In our credentialxcms.Rdata, xcms is run together with three technical replicates each.

load(file = "credentialxcms.Rdata")
# if export=T, feature tables and a histogram summary of retention time shift, ppm error, peakwidth will be exported.
features <- credential::getXCMSfeature(xs = credentialxcms, intchoice="into", sampling = 1, sampleclass = NULL, export = T)
feature1t1 <- features$`1T1-credTable`
feature1t2 <- features$`1T2-credTable`

# OR upload features tables in csv format
feature1t1 <- data.table(read.csv(system.file("extdata","features1T1.csv", package = "credential")))
feature1t2 <- data.table(read.csv(system.file("extdata","features1T2.csv", package = "credential")))

# automatic credentialing 
credential_test <- credential::credentialing(peaktable1 = feature1t1, peaktable2 = feature1t2, ppm = 15, rtwin = 1, rtcom =2, ratio1 = 1/1, ratio2 = 1/2, 
                                             ratio_tol = 0.1, ratio_ratio_tol = 0.9, cd = 13.00335-12, charges = 1:4, mpc = c(12,120), maxnmer = 4,
                                             export = T, plot = T, projectName = "credential_demo")
  # function help
  help("credential::credentialing")

# manual step-by-step credentialing
  
  # parameter settings
  ppm = 15
  rtwin = 1
  rtcom =2
  ratio1 = 1/1
  ratio2 = 1/2
  ratio_tol = 0.1
  ratio_ratio_tol = 0.9
  cd = 13.00335-12
  charges = 1:4
  mpc = c(12,120)
  maxnmer = 4
  projectName = "credential_demo"
  
  # step1 find isotope knots of each feature table
  knots1 <- credential::findknots(features = feature1t1, .zs = charges, ppmwid = ppm, rtwid = rtwin, cd = cd)
  knots2 <- credential::findknots(features = feature1t2, .zs = charges, ppmwid = ppm, rtwid = rtwin, cd = cd)
  
  # step2 resolve merged isotope knots
  knots1 = fixmergedknots(knots1, feature1t1)
  knots2 = fixmergedknots(knots2, feature1t2)
  
  # step3 credential knots from each feature table (quipus)
  credentialedknots1 <- credential::credentialknots(Knots = knots1, ppmwid = ppm, rtwid = rtwin, Ratio = ratio1, Ratio.lim = ratio_tol)
  credentialedknots2 <- credential::credentialknots(Knots = knots2, ppmwid = ppm, rtwid = rtwin, Ratio = ratio2, Ratio.lim = ratio_tol)
  
  # step4 match quipus (credentialed knots) to obtain credentialed groups
  credentialedquipus <- credentialquipu(credentialedknots1, credentialedknots2, ppm = ppm, rtwin = rtcom, ratio_ratio = ratio1/ratio2, ratio_ratio_tol = ratio_ratio_tol, tailmatch=T)
  
  # step5 plot credentialed peaks
    
    # merge credentialed peaks
    credpeak1 = feature1t1[knots1$cc_knot[credentialedknots1$knot_quipu[!is.na(quipu)],,on="knot"],,on="cc"][credentialedquipus$credentialedgroups[,.(quipu=quipu1,charge1,mainmz11,mainmz21,ratio1)],,on="quipu"]
    credpeak1 = credpeak1[,.(cc1=cc,mz1=mz,rt1=rt,i1=i,knot1=knot,tail1=tail,quipu1=quipu,charge1,mainmz11,mainmz21,ratio1)]
    credpeak2 = feature1t2[knots2$cc_knot[credentialedknots2$knot_quipu[!is.na(quipu)],,on="knot"],,on="cc"][credentialedquipus$credentialedgroups[,.(quipu=quipu2,charge2,mainmz12,mainmz22,ratio2,ratio1_ratio2)],,on="quipu"]
    credpeak2 = credpeak2[,.(cc2=cc,mz2=mz,rt2=rt,i2=i,knot2=knot,tail2=tail,quipu2=quipu,charge2,mainmz12,mainmz22,ratio2,ratio1_ratio2)]
    credentialedpeaks = do.call(rbind,apply(credentialedquipus$credentialedindex[order(credentialedquipus$credentialedgroups$basemz1)], MARGIN = 1, function(x){credential:::cbind.fill(credpeak1[quipu1==x[1]][order(mz1)],credpeak2[quipu2==x[2]][order(mz2)])}))
    # plotting
    plotcredpeaks(Credentialedindex = credentialedquipus$credentialedindex, Credentialedpeaks = credentialedpeaks,
                  filename = paste(paste(projectName,nrow(credentialedquipus$credentialedindex),"credentialed_peak_groups",sep="_"),".pdf",sep=""))
  
  # step5 export credentialed groups and credentialed peaks as csv files
  write.csv(credentialedquipus$credentialedgroups,file = paste0(projectName,"_CredentialedGroups.csv"))
  write.csv(credentialedpeaks,file = paste0(projectName,"_CredentialedPeaks.csv"))
    
# useful information
  library(magrittr)
  
  # check "Number of peaks per knot"
  credential_test$knots1$cc_knot$knot %>% table %>% table
  
  # check "Charge states of isotopic knots"
  credential_test$knots1$knot$z %>% table

  # check "Ppm error within knots"
  ppm_knot = feature1t1[credential_test$knots1$cc_knot,,on="cc"][,mean(max(mz%%1.003355)-min(mz%%1.003355))/mean(mz)*1E6,by="knot"]$V1
  hist(ppm_knot, breaks = 100)
  
  # check "Ppm error within quipus"
  ppm_quipu = credential_test$knots1$knot[credential_test$credentialedknots1$knot_quipu,,on="knot"][!is.na(quipu),(max(meanr)-min(meanr))/mean(meanmz)*1E6,by="quipu"]$V1
  hist(ppm_quipu, breaks=100)
  
  # check "Number of knots per quipu"
  credential_test$credentialedknots1$knot_quipu$quipu %>% table %>% table

  # check quipus with iso support (more than one isotopic peaks per knot)
  feature1t1[credential_test$knots1$knot[credential_test$knots1$cc_knot[credential_test$credentialedknots1$knot_quipu[credential_test$credentialedknots1$quipu_stat[minsupport > 1],,on="quipu"],,on="knot"],,on="knot"],,on="cc"][,.SD[mz == min(mz)],by="quipu"][,.(cc, mz, rt, i, z, knot, quipu)]
  ß
  # check quipus without iso support (just one isotopic peak per knot)
  feature1t1[credential_test$knots1$knot[credential_test$knots1$cc_knot[credential_test$credentialedknots1$knot_quipu[credential_test$credentialedknots1$quipu_stat[minsupport == 1],,on="quipu"],,on="knot"],,on="knot"],,on="cc"][,.SD[mz == min(mz)],by="quipu"][,.(cc, mz, rt, i, z, knot, quipu)]
  

