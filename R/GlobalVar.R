# Global variables of all functions
## getXCMSfeature.R 
utils::globalVariables(names = c("rtmax","rtmin","mzmax","mzmin","mzppm","..density..","rtdiff","pwidthdiff"))
## credential_knots.R
utils::globalVariables(names = c("z","direction","mdig","maxi","meanr"))
utils::globalVariables(names = c("g1","z","g2","meanr","rrtg","knot","gtemp",
                                 "meanmz","n","quipu","ratio","minsupport","maxsupport",
                                 "nknot","basemz","mainmz","maxi","ncar","mainmz2","mainmz1","charge"))
## credential_main.R
utils::globalVariables(names = c("quipu1","charge1","mainmz11","mainmz21","cc","i","knot",
                                 "quipu2","charge2","mainmz12","mainmz22","ratio1_ratio2",
                                 "mz1","mz2","ncar1","ncar2"))
## credential_quipu.R
utils::globalVariables(names = c("ratio1","ratio2","ratio1_ratio2","basemz1","credentialed"))
## find_knots.R
utils::globalVariables(names = c("z","i","n","cc","c13r","g4","g3","c13c"))
## fixedmergedknot.R
utils::globalVariables(names = c("cc","knot","i"))
## plotting.R
utils::globalVariables(names = c("knot","i","meanr","z","quipu","quipu1","mz1","rt1","i1",
                                 "ratio1","quipu2","mz2","rt2","i2","ratio2","quipu",
                                 "ncar1","ncar2","charge","label"))