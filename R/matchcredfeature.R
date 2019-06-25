
matchcredfeature = function(dt1T1, dt1T2, ppm, drt, ratio, ratio_tol = 0.1, export) {

match_combined=data.table()
nomatch1_combined=data.table()
nomatch2_combined=data.table()
quipu_match=data.table()

quipu_m1 = numeric()
quipu_m2 = numeric()
rtmean_m1 = numeric()
rtmean_m2 = numeric()
basemz_m1 = numeric()
basemz_m2 = numeric()
quipu_nm1 = numeric()
quipu_nm2 = numeric()
quipu_dm1 = numeric()
quipu_dm2 = numeric()
match_idx2 =numeric()
fratio = numeric()

# Give the accurate mass range:

dt1T1[,"maxTmz":=mz/(1-ppm/1e6)]
dt1T1[,"minTmz":=mz/(1+ppm/1e6)]

dt1T2[,"maxTmz":=mz/(1-ppm/1e6)]
dt1T2[,"minTmz":=mz/(1+ppm/1e6)]


for (tmp.quipu1 in unique(dt1T1[,quipu])){

  idx1 <- which(dt1T1[,"quipu"]==tmp.quipu1)
  cat("\rSearching common credentialed featuregroup in group#2 for #",tmp.quipu1, " featuregroup in group#1 ")

  basemz1 = min(dt1T1[idx1,mz])
  maxmz1 = basemz1/(1-ppm/1e6)
  minmz1 = basemz1/(1+ppm/1e6)
  basert1 = mean(dt1T1[idx1,rt])
  minrt1 = basert1-drt
  maxrt1 = basert1+drt

  idx2 <- c(which(dt1T2[,"maxTmz"]>=minmz1 & dt1T2[,"minTmz"]<=minmz1 & dt1T2[,"rt"]>=minrt1 & dt1T2[,"rt"]<=maxrt1), which(dt1T2[,"maxTmz"]>=maxmz1 & dt1T2[,"minTmz"]<=maxmz1 & dt1T2[,"rt"]>=minrt1 & dt1T2[,"rt"]<=maxrt1))

  if(length(idx2)<1){
    #cat("\nNo match is found for quipu#",tmp.quipu1, "feature")
    nomatch1_combined = rbind(nomatch1_combined,dt1T1[idx1])
  }
  if(length(idx2)>1){
    tmp.quipu2 = unique(dt1T2[idx2,quipu])
    #cat("\nDuplicate match between quipu#:",tmp.quipu1,"in group #1 and quipu#:",tmp.quipu2, "in group #2")
    quipu_dm1 = c(quipu_dm1,tmp.quipu1)
  }
  if(length(idx2)==1){

    tmp.quipu2 <- dt1T2[idx2,quipu]
    if(sum(tmp.quipu2==quipu_m2)>=1){
    quipu_dm2 = c(quipu_dm2,tmp.quipu2)
    #cat("\nDuplicate match between quipu#:",tmp.quipu1,"in group #1 and quipu#:",tmp.quipu2,"in group #2")
    }else{
    basemz2 = dt1T2[idx2,mz]
    idx2 <- which(dt1T2[,"quipu"]==tmp.quipu2)  #extract all quipu in set 2
    merged = cbind.fill(dt1T1[idx1,],dt1T2[idx2,])
    basert2 = mean(dt1T2[idx2,rt])
    quipu_m1 = c(quipu_m1,tmp.quipu1)
    quipu_m2 = c(quipu_m2,tmp.quipu2)
    #cat("\nUnique match is found between quipu#:",tmp.quipu1,"in group #1 and quipu#:",tmp.quipu2,"in group #2.")

    basemz_m1 = c(basemz_m1,basemz1)
    basemz_m2 = c(basemz_m2,basemz2)
    rtmean_m1 = c(rtmean_m1,basert1)
    rtmean_m2 = c(rtmean_m2,basert2)
    match_combined = rbind(match_combined,merged)
    match_idx2 = c(match_idx2,idx2)
    }
   }
}


quipu_match = data.table(cbind(quipu_m1,rtmean_m1,basemz_m1,quipu_m2,rtmean_m2,basemz_m2))
nomatch2_combined = dt1T2[!match_idx2]
nomatch1 = length(unique(dt1T1[,quipu]))-length(quipu_match[,quipu_m1])
nomatch2 = length(unique(dt1T2[,quipu]))-length(quipu_match[,quipu_m2])

#rename the matched feature table
match_combined[,c("maxTmz","minTmz"):= NULL]
match_combined[,c("maxTmz","minTmz"):= NULL]
colnames(match_combined) <- c("cc_1","mz_1", "rt_1", "int_1", "knot_1", "tail_1","quipu_1","ratio_1","cc_2","mz_2", "rt_2", "int_2", "knot_2", "tail_2","quipu_2","ratio_2")
match_combined[,"combined_ratio":=ratio_1/ratio_2]
match_filtered= match_combined[which(combined_ratio>=ratio*ratio_tol & combined_ratio<=ratio/ratio_tol),]

cat("\n",length(unique(match_filtered[,quipu_1])), "common feature groups are finally credentialed.\n")

#cat("\n",length(quipu_match[,quipu_m1]), " unique matches are found between two groups.\n")

#cat(nomatch1, "unmatched features in group #1 and",nomatch2,"unmatched features in group #2.\n")

#cat(length(quipu_dm1),"features in group #1 and",length(quipu_dm2),"features in group #2 have duplicated matches.\n")

list(Credentialed_FeatureGroups=quipu_match,Credentialed_Features=match_combined,Credentialed_Features_Filtered = match_filtered, NomatchFeatures_Group1=nomatch1_combined,NomatchFeatures_Group2=nomatch2_combined)
}
