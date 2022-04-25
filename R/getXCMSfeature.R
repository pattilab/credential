#' Extract peak-table from an xcmsSet object
#'
#'The adaptive function to automatically extract features (after grouping) from xcmsSet object
#'@usage getXCMSfeature(xs,intchoice=c("into","intb","maxo"),sampling=1, sampleclass = NULL, export = F)
#'@description The function extracts features table from xcmsSet object and return two types of data.table, one for reference, one for credentialing input
#'@param xs xcmsSet object
#'@param intchoice character the intensity options provided by xcms. Options: \code{into}(default) -raw peak area, \code{intb}-baseline corrected peak area, \code{maxo}- max MS peak height.
#'@param sampling numeric The method to extract peak intensity in each feature group. If 1, sampling mean intensities of all peaks. If 2, sampling the maximum intensity of all peaks.
#'@param sampleclass character Specify sample classes to be extracted. The default is NULL, where all sample classes in "phenoData" table of xcmsSet.
#'@param export TRUE/FALSE Specify if the output features table should be exported to csv files and statistics exported to pdf file.
#'@import data.table xcms gridExtra ggplot2 utils
#'@keywords credential getXCMSfeature getxcmsfeature xcms
#'@return a list The output includes two sets of data.table objects. \code{sampleclass-featuretable} data.table The features table from xcmsSet.
#'\code{sampleclass-credTable} data.table The structured features table for credentialing analysis.
#'@author Lingjue "Mike" Wang <wang.lingjue@wustl.edu>
#'@export

getXCMSfeature = function(xs, intchoice=c("into","intb","maxo"),sampling=1, sampleclass = NULL, export = F){

  if(class(xs)[1]!="xcmsSet") {
    stop("The first input is not a xcmsSet object.")
  }
  intchoice = match.arg(intchoice)
  # test samplename existance
  if(is.null(sampleclass)) {
    sampleclass = unique(as.character(xs@phenoData$class))
    cat("sampleclass is not assigned. Using sample classes in phenoData:", sampleclass)
  } else if (any(sampleclass == xs@phenoData$class)) {
    sampleclass = unique(intersect(sampleclass,xs@phenoData$class))
    cat("Following sampleclass(es) match xs@phenoData:", sampleclass)
  } else {
    stop("The sampleclass doesn't match any classes in the xcmsSet object.")
  }

  if(sampling==1) {
    intfun<- function(x) mean(x)
    samp_method = "mean intensty"
  } else if (sampling == 2) {
    intfun<- function(x) max(x)
    samp_method = "maximum intensity"
  } else {
    stop("Sampling parameter is not valid.")
  }
  cat("\nExtracting feature information with",samp_method, "by sampleclass:",unique(sampleclass),"...\n" )

  #intchoice = "into" # for testing
  #sampling = 1 # for testing
  #xs = xs_1 # for testing

# Extract the intensity based on intchoice and sampling
  groups = data.table(copy(xs@groups))
  peaks = data.table(copy(xs@peaks))
  sampleid = list()
  featuretable = list()
  credTable = list()

  for (class in sampleclass) {
    #class = sampleclass[1] # for testing

    tmp_feature = cbind(index=seq.int(nrow(groups)),groups[,c("mzmed","rtmed")])
    sampleid[class] = list(which(xs@phenoData==class))
    class_ind = which(groups[,class, with=F]>=1)

    for (ind in class_ind) {
    #ind =  class_ind[1] #testing

    id = xs@groupidx[[ind]]
    tmp_peak = peaks[id,][sample %in% sampleid[[class]]][,':='(peakwidth = rtmax-rtmin, mzppm = (mzmax-mzmin)*1E6/mz)]
    tmp_feature[ind,':='(int = intfun(unlist(tmp_peak[,intchoice, with=F])), npeaks = nrow(tmp_peak), mzppm = mean(tmp_peak$mzppm), rtdiff = max(tmp_peak$rt)-min(tmp_peak$rt),
                         peakwidth = mean(tmp_peak[,peakwidth]), pwidthdiff = max(tmp_peak$peakwidth)-min(tmp_peak$peakwidth))]
    }

    featuretable[class] <- list(tmp_feature[class_ind,])
    credtable = tmp_feature[class_ind,][,c("index","mzmed","rtmed","int")]
    colnames(credtable) <- c("cc","mz","rt","i")
    credTable[class] <- list(credtable)
    cat(nrow(featuretable[[class]])," features are extracted from class:", class, "\n")
    cat("\tAverage mass error is",round(mean(featuretable[[class]]$mzppm),1), "ppm", paste0("(min: ", round(min(featuretable[[class]]$mzppm),1)," ppm, max: ", round(max(featuretable[[class]]$mzppm),1)," ppm)."),"\n")
    cat("\tAverage retention time shift is",round(mean(featuretable[[class]]$rtdiff),1), "second", paste0("(min: ", round(min(featuretable[[class]]$rtdiff),1)," second, max: ", round(max(featuretable[[class]]$rtdiff),1)," second)."),"\n")
    cat("\tAverage peakwidth is",round(mean(featuretable[[class]]$peakwidth),1), "second", paste0("(min: ", round(min(featuretable[[class]]$peakwidth),1)," second, max: ", round(max(featuretable[[class]]$peakwidth),1)," second)."),"\n")
    cat("\tAverage peakwidth variance is",round(mean(featuretable[[class]]$pwidthdiff),1), "second", paste0("(min: ", round(min(featuretable[[class]]$pwidthdiff),1)," second, max: ", round(max(featuretable[[class]]$pwidthdiff),1)," second)."),"\n")
  }

  # dataset clean-up
  if(export){
    for (class in sampleclass){
    filename <- paste("xcms",class,paste0(nrow(featuretable[[class]]),"feature",sep=""),Sys.Date(), sep="_")
    write.csv(featuretable[[class]], paste0(filename,".csv",sep=""))

    pdf(file=paste(filename,".pdf",sep=""), width = 8.5, height = 11)
    temp_grid <- list()
    temp_grid[[1]] <- ggplot2::ggplot(featuretable[[class]], aes(x=mzppm)) +
                                geom_histogram(aes(y=..density..), bins = 20, colour="black", fill="white")+ geom_density(alpha=.2, fill="#FF6666")  +
                                theme_classic() + ggtitle(paste(class,"Distribution of ppm in",nrow(featuretable[[class]]),"features"))

    temp_grid[[2]] <- ggplot2::ggplot(featuretable[[class]], aes(x=rtdiff)) +
                                      geom_histogram(aes(y=..density..), bins = 20, colour="black", fill="white")+ geom_density(alpha=.2, fill="#FF6666")  +
                                      theme_classic() + ggtitle(paste(class,"Distribution of retention time shift in",nrow(featuretable[[class]]),"features"))

    temp_grid[[3]] <- ggplot2::ggplot(featuretable[[class]], aes(x=peakwidth)) +
                                      geom_histogram(aes(y=..density..), bins = 20, colour="black", fill="white")+ geom_density(alpha=.2, fill="#FF6666")  +
                                      theme_classic() + ggtitle(paste(class,"Distribution of peakwidth in",nrow(featuretable[[class]]),"features"))

    temp_grid[[4]] <- ggplot2::ggplot(featuretable[[class]], aes(x=pwidthdiff)) +
                                      geom_histogram(aes(y=..density..), bins = 20, colour="black", fill="white")+ geom_density(alpha=.2, fill="#FF6666")  +
                                      theme_classic() + ggtitle(paste(class,"Distribution of peakwidth varince in",nrow(featuretable[[class]]),"features"))

    do.call(gridExtra::grid.arrange,c(temp_grid,nrow=4,ncol=1))
    dev.off()
    }
    cat("\n Features table and summary is exported under:\n",getwd())
  }

  extractFeatures = c(featuretable,credTable)
  names(extractFeatures) <- c(paste0(sampleclass,"-featuretable"),paste0(sampleclass,"-credTable"))
  cat("\nFeature extraction are done.")
  return(extractFeatures)

}
