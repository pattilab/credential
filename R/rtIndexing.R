roughlyCoelutingPeakIndices = function(
  rt, 
  the_rt_index,
  isotope_rt_delta_s
) {  #only preform search on coeluting peaks - optimized with rough exclusion
  
  middle_index = round(as.numeric(rt)/(isotope_rt_delta_s/2)) +1
  
  start = if ((middle_index-3) >= 1) middle_index-3 else 1
  end = if ((middle_index+3) <= nrow(the_rt_index)) middle_index+3 else nrow(the_rt_index)
  
  return(the_rt_index[start,2]:the_rt_index[end,2])
}

indexRt = function(
  peak_analysis, 
  isotope_rt_delta_s
) {
  if (any(order(peak_analysis[,"rt"]) != 1:nrow(peak_analysis))) { stop("Peaks must be ordered by retention time before building the retention time index.") } 
  
  trunc_rts = trunc(peak_analysis[,"rt"])
  unique_rts = unique(trunc_rts)
  indices = unlist(lapply(unique_rts, function(x){
    which(trunc_rts == x)[1]
  }))
  rt_index = cbind(unique_rts, indices)
  desired_breaks = seq(0, max(peak_analysis[,"rt"]), isotope_rt_delta_s/2)
  
  the_rt_index = matrix(unlist(lapply(desired_breaks, function(x) {
    foo = which(abs(x - rt_index[,1]) == min(abs(x - rt_index[,1])))[1]
    return(c(x, rt_index[foo,2]))
  })),ncol=2, byrow=T)
  
  last_val= the_rt_index[nrow(the_rt_index),1]+isotope_rt_delta_s/2
  the_rt_index = rbind(the_rt_index, c(last_val,nrow(peak_analysis)))
  
  return(the_rt_index)
}