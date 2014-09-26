matrixLapplyUniqueRows = function(m, unique_column, the_fun) {
  cat("\nStarting rows:", nrow(m), "- Unique", paste(sep="",unique_column,":"), length(unique(m[,unique_column])))
  
  unique_rows = unique(m[,unique_column])
  cat("\nWorking on unique row of", length(unique_rows), "total:\n")
  
  tmp = lapply(1:length(unique_rows), function(i) {
    cat("\r", i); flush.console();
    x = unique_rows[i]
    
    grouped_rows = m[m[,unique_column] == x,,drop=F]
    return(the_fun(grouped_rows))
  })
  
  m2 = do.call("rbind", tmp)
  cat("\nEnding rows:", nrow(m2), "- Unique", paste(sep="",unique_column,":"),length(unique(m2[,unique_column])))
  
  m2
}