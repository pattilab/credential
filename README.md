# credential release/3.1
## Introduction
Credentialing is a method to identify biologically relevant features from untargeted metabolomics data. The analysis is based on two feature tables with index (*cc*), m/z (*mz*), retention time (*rt*) and intensity (*i*) values. Experimentally, credentialing requires two sets of samples by mixing unlabeled and U13C-labeled samples at two user-defined ratios (*ratio1*, *ratio2*). *credential* is a bioinformatics tool (R package) to analyze feature tables from untargeted LC/MS based analysis and infer biologically-relavent features from untargeted metabolomics. 

The current version **credential_3.1.7** was developed as a new algorithm to orginal methods paper:

1. First, the algorithm searches *knot* at defined possible charge states in each feature table. A *knot* is a set of putative isotopologues with close retention time and defined m/z spacing. 
2. Second, the algorithm resolves merged knots with multiple isotopologue patterns. 
3. Third, the algorithm searches *quipu* from *knots*. A *quipu* is a set of *knots* that represent the unlabeled and labeled sides of the isotopologues of one compound. Features within a quipu have shared retention time, m/z spacing, charge state, intensity pattern (ratios), qualifying mass per carbon and peak intensity ratios. *quipus* at each mixing ratios are considered feature groups that pass the first round of credentialing. 
4. Finally, the algorithm matches credentialed quipus between two credentialing groups by retention time, charge state, head/tail m/z values, and thresholds of the ratios of major labeled/unlabled peaks. The resulting *quipus* are considered credentialed features that pass two rounds of credentialing.

## References
- [Mahieu, N. G., Huang, X., Chen, Y. J., & Patti, G. J. (2014). Credentialing features: a platform to benchmark and optimize untargeted metabolomic methods. Analytical Chemistry, 86(19), 9583-9589.](https://doi.org/10.1021/ac503092d)
- [Wang, L., Naser, F. J., Spalding, J. L., & Patti, G. J. (2019). A protocol to compare methods for untargeted metabolomics. In Metabolic Signaling (pp. 1-15). Humana Press, New York, NY.](https://link.springer.com/protocol/10.1007/978-1-4939-8769-6_1)

## Usage

### Installation
```
install.packages("devtools")
devtools::install_github("pattilab/credential", ref="release/3.1", dependencies = T)
```

### Data Input
Credentialing takes a feature table (*data.table* or *data.frame*) that requires at least 4 columns for credentialing algorithm: *cc* - feature index, *mz* - mass-to-charge value, *rt* - retention time, and *i*-intensity. Feature table can be generated from XCMS or other peak detection software. 

- To generate the formatted feature tables from *xcmsSet*:
```
features <- credential::getXCMSfeature(xs = credentialxcms, intchoice="into", sampling = 1, sampleclass = NULL, export = T)
feature1t1 <- features$`1T1-credTable`
feature1t2 <- features$`1T2-credTable`
```
- To import feature tables from external CSV files:
```
feature1t1 <- data.table(read.csv(system.file("extdata","features1T1.csv", package = "credential")))
feature1t2 <- data.table(read.csv(system.file("extdata","features1T2.csv", package = "credential")))
```
Column names must be adjusted accordingly as above using *colnames()*.

```
colnames(features) <- c("cc","mz","rt","i")
```

```
> feature1t1
          cc        mz        rt          i
    1:     1  101.0963 2182.1550  14444.594
    2:     2  101.0963 2227.5400  18668.742
    3:     3  101.0964 2258.0550  17612.537
    4:     4  101.0964 2280.3100  32333.016
    5:     5  102.0455   82.5735   6732.472
   ---                                     
14550: 15448 1597.2860 1551.5900  47496.896
14551: 15449 1597.2964 1578.8200  93539.010
14552: 15450 1598.1958 1605.5200   5061.292
14553: 15451 1598.2902 1551.5900  60039.320
14554: 15452 1598.3012 1577.8100 206334.273
```

### credentialing 3.0 -- integrated processing
```
credential_test <- credential::credentialing(peaktable1 = feature1t1, peaktable2 = feature1t2, ppm = 15, rtwin = 1, rtcom =2, ratio1 = 1/1, ratio2 = 1/2, 
                                             ratio_tol = 0.1, ratio_ratio_tol = 0.9, cd = 13.00335-12, charges = 1:4, mpc = c(12,120), maxnmer = 4,
                                             export = T, plot = T, projectName = "credential_demo")
```
- Documentation of *credentialing()*
```
help("credential::credentialing")
```
### credentialing 3.0 -- step-by-step processing

```
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
  
 # step1 find isotope knots of each feature table (resolve merged isotope knots are performed in this step now)
 knots1 <- credential::findknots(features = feature1t1, .zs = charges, ppmwid = ppm, rtwid = rtwin, cd = cd)
 knots2 <- credential::findknots(features = feature1t2, .zs = charges, ppmwid = ppm, rtwid = rtwin, cd = cd)
  
 # step2 credential knots from each feature table (quipus)
 credentialedknots1 <- credential::credentialknots(Knots = knots1, ppmwid = ppm, rtwid = rtwin, Ratio = ratio1, Ratio.lim = ratio_tol)
 credentialedknots2 <- credential::credentialknots(Knots = knots2, ppmwid = ppm, rtwid = rtwin, Ratio = ratio2, Ratio.lim = ratio_tol)
  
 # step3 match quipus (credentialed knots) to obtain credentialed groups
 credentialedquipus <- credentialquipu(credentialedknots1, credentialedknots2, ppm = ppm, rtwin = rtcom, ratio_ratio = ratio1/ratio2, 
                                       ratio_ratio_tol = ratio_ratio_tol, tailmatch=T)
 
 # step4 plot credentialed peaks
    # merge credentialed peaks
    credpeak1 = feature1t1[knots1$cc_knot[credentialedknots1$knot_quipu[!is.na(quipu)],,on="knot"],,on="cc"][credentialedquipus$credentialedgroups[,. (quipu=quipu1,charge1,mainmz11,mainmz21,ratio1)],,on="quipu"]
    credpeak1 = credpeak1[,.(cc1=cc,mz1=mz,rt1=rt,i1=i,knot1=knot,tail1=tail,quipu1=quipu,charge1,mainmz11,mainmz21,ncar1,ratio1)]
    credpeak2 = feature1t2[knots2$cc_knot[credentialedknots2$knot_quipu[!is.na(quipu)],,on="knot"],,on="cc"][credentialedquipus$credentialedgroups[,.(quipu=quipu2,charge2,mainmz12,mainmz22,ratio2,ratio1_ratio2)],,on="quipu"]
    credpeak2 = credpeak2[,.(cc2=cc,mz2=mz,rt2=rt,i2=i,knot2=knot,tail2=tail,quipu2=quipu,charge2,mainmz12,mainmz22,ncar2,ratio2,ratio1_ratio2)]
    credentialedpeaks = do.call(rbind,apply(credentialedquipus$credentialedindex[order(credentialedquipus$credentialedgroups$basemz1)], MARGIN = 1, function(x){credential:::cbind.fill(credpeak1[quipu1==x[1]][order(mz1)],credpeak2[quipu2==x[2]][order(mz2)])}))
    
    # plotting
    plotcredpeaks(Credentialedindex = credentialedquipus$credentialedindex, Credentialedpeaks = credentialedpeaks,
                  filename = paste(paste(projectName,nrow(credentialedquipus$credentialedindex),"credentialed_peak_groups",sep="_"),".pdf",sep=""))
  
 # step5 export credentialed groups and credentialed peaks as csv files
 write.csv(credentialedquipus$credentialedgroups,file = paste0(projectName,"_CredentialedGroups.csv"))
 write.csv(credentialedpeaks,file = paste0(projectName,"_CredentialedPeaks.csv"))
```

### data output
- *credentialed* -- *list* Including all the processed results of credentialing
- *credentialedgroups* -- *data.table* A summary of credentialed *quipus* with following information: *quipu* (quipu index), *nknot* (number of knots), *npeak* (number of peaks), *rtmean* (mean retention time), *basemz* (lowest mz), *mainmz1* (credentialed pair - unlabeled mz), *mainmz2* (credentialed pair - labeled mz), *int1* (credentialed pair - intensity of unlabeled mz), *int2* (credentialed pair - intensity of labeled mz), *ratio* (ratio of unlabled/labeled intensity), *credentialed* (TRUE/FALSE)
- *credentialedpeaks* *data.table* -- A summary of all credentialed peaks in two credentialing groups with following information: *quipu* (quipu index), *knot* (knot index), *cc* (feature index), *charge* (charge state if determined), *mz* (m/z), *rt* (retention time), *i* (intensity), *tail* (if this feature is tail peak in knot), *mainmz1*, *mainmz2*, *ratio*
- *credentialedknots1* *list* -- A list of credentailed *quipus* (*knots*) in first credentialing group
- *credentialedknots1* *list* -- A list of credentailed *quipus* (*knots*) in second credentialing group
- *knots1* *list* -- All knots in first credentialing group (*feature1t1*) 
- *knots2* *list* -- All knots in second credentialing group (*feature1t2*) 
- *CredentialParams* *list* -- Parameters used in this analysis.

## Author of README (Updated on 04/26/2022)
Lingjue Mike Wang (wang.lingjue@wustl.edu)
