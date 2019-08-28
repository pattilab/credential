#' Diagnostic Plotting Function for Knot
#'
#' @description This function generate plots of knots with m/z   
#' @usage plot.knots(.knot, knots, features, scales = c(1,1), ...)
#' @param .knot numeric Knot index.
#' @param knots list Knots. The output of \code{\link{findknots}}
#' @param features numeric The maximum retention time difference in seconds.
#' @param scales numeric Vector containing the minimum and maximum allowable mass per carbon.
#' @param ... Further arguments to be passed.
#' @import data.table magrittr ggplot2 gridExtra
#' @return a ggplot object
#' 

plot.knots = function(.knot, knots, features, scales = c(1,1), ...) {
    
    df = features[knots$cc_knot,,on="cc"][knots$knot[knot%in%.knot],,on="knot", nomatch=0] %>% data.table
    xlim = c(min(df$mz)-15, max(df$mz)+15)
    ylimmeanr = c(mean(df$meanr)-scales[1]*1, mean(df$meanr) + scales[1]*1)
    ylimrt = c(mean(df$rt)-scales[2]*1, mean(df$rt) + scales[2]*1)
    
    spec = ggplot(df) + geom_segment(aes(x = mz, xend = mz, y = i, yend = 0, colour=factor(knot))) + ggtitle(paste(.knot)) + theme_classic() + theme(legend.position="none") + geom_segment(aes(y = 0, yend = 0, x = mean(range(df$mz)) - 5, xend = mean(range(df$mz)) + 5)) + coord_cartesian(xlim = xlim)
    mrspec = ggplot(df) + geom_point(aes(x = mz, y = meanr, colour=factor(knot))) + theme_classic() + theme(legend.position="none") + coord_cartesian(xlim = xlim, ylim=ylimmeanr)
    rtspec = ggplot(df) + geom_point(aes(x = mz, y = rt, colour=factor(knot))) + theme_classic() + theme(legend.position="none") + coord_cartesian(xlim = xlim, ylim = ylimrt)
    
    psdf = knots$knot[knot%in%.knot] %>% as.data.frame
    distspec = ggplot(psdf) + geom_point(aes(x = rt, y = meanr, colour = factor(knot), shape = factor(z))) + theme_classic() + theme(legend.position="none") + coord_cartesian(xlim = ylimrt, ylim = ylimmeanr)
    gridExtra::grid.arrange(mrspec, spec, rtspec, distspec, nrow = 4, newpage = TRUE)
    
}

#' Diagnostic Plotting Function for Quipu
#'
#' @description This function generate plots of quipu (credentialed knots).
#' @usage plot.quipu(.quipu=NULL, quipus, knots, features, export = TRUE, file = "quipu.pdf", ...) {
#' @param .quipu numeric quipu index.
#' @param knots list Knots. The output of \code{\link{findknots}}
#' @param quipus list Quipus aka credentialed knots. The output of \code{\link{credentialknots}}
#' @param features data.table The featue table.
#' @param export logical Whether to export or not.
#' @param file file name of the export pdf.
#' @param ... Further arguments to be passed.
#' @import data.table magrittr ggplot2 gridExtra grDevices
#' @seealso \code{\link{plot.knots}}
#' @return This function is a plotting function with no output but an exported pdf file.
#'
plot.quipu = function(.quipu=NULL, quipus, knots, features, export = TRUE, file = "quipu.pdf", ...) {
  
  if(is.null(.quipu))
    .quipu = unique(quipus$knot_quipu$quipu)
  if(export){
    pdf(file = file, width = 8.5, height = 11)
    for(q in .quipu) {
    cat("\rGenerating diagnostic plot of credentialed feature #", q, "...")
    plot.knots(quipus$quipu[quipus$knot_quipu[quipu %in% q],,on="quipu"]$knot, knots, features)
    }
    dev.off()
  }
}

#' Plotting all credentialed peaks
#'
#' @description This function generates pseudo MS spectrum of all credentialed peak
#' groups. Each plot includes two pseudo-MS spectrum with each centroid peak represent a
#' credentialed isotopologue (the i refers to the intensity from feature table). The red spectra (up) is the
#' credentialed peak group in first ratio condition, the blue spectra (down) is the credentialed
#' peak group in second ratio condition. The legend represent the quipu number corresponding to
#' the matched credentialed groups in first and second conditions.
#' @usage plotcredpeaks(Credentialedindex, Credentialedpeaks, filename = "credentialedpeaks.pdf", ...) {
#' @param Credentialedindex data.table The index table of quipu-to-quipu assignment. A output table from \code{\link{credentialquipu}} or \code{\link{credentialing}}
#' @param Credentialedpeaks data.table The peak table includes the aligned peaks of all credentialed groups. A output table from \code{\link{credential}}
#' @param filename character the file name of the export pdf file Include .
#' @param ... Further arguments to be passed.
#' @import data.table ggplot2 gridExtra grDevices
#' @seealso \code{\link{findknots}}
#' @keywords credential plotcredpeaks
#' @return This plotting function have no functional output but an exported pdf file under working directory.
#' @export
#'

plotcredpeaks = function(Credentialedindex, Credentialedpeaks, filename = "credentialedpeaks.pdf", ...) {

  cat("\nPlotting",nrow(Credentialedindex),"credentialed peak groups...")
  index = Credentialedindex
  peaks = data.table(Credentialedpeaks)

  pdf(file= filename, width = 8.5, height = 11)
  temp_grid <- list()
  j=1
  for (i in seq.int(nrow(index))) {
  last = nrow(index)
  q1 = index$quipu1[i]
  q2 = index$quipu2[i]
  peaks1 = peaks[quipu1==q1,.(quipu=quipu1,mz=mz1,rt=rt1, i=i1,ratio=ratio1)]
  peaks2 = peaks[quipu2==q2,.(quipu=quipu2,mz=mz2,rt=rt2, i=-i2,ratio=ratio2)]
  ratio_ratio = peaks1$ratio[1]/peaks2$ratio[1]
  peak = rbind(peaks1,peaks2)[,quipu:=as.factor(quipu)][]
  
  title = paste("#",i,": basemz ",round(min(peak$mz),4),", rt ", round(mean(peak$rt),2),", ratio ", round(ratio_ratio,2),sep="")

  plot <- ggplot(peak, aes(x=mz,y=i,color=quipu, label=mz)) + geom_bar(width = 0.04, stat="identity",fill="white")
  grob <- ggplotGrob(plot)
  ax <- grob[["grobs"]][grob$layout$name == "axis-b"][[1]]
  plot <- plot + annotation_custom(grid::grobTree(ax, vp = grid::viewport(y=1, height=sum(ax$height))),
                         ymax=0, ymin=0) +
          geom_hline(aes(yintercept=0)) +
          theme_classic() + labs(title=title) +
          theme(axis.text.x = element_blank(),
                plot.title = element_text(hjust=0.5, size=10, face="bold"),
                axis.line.x = element_blank(),
                axis.ticks.x=element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=1))
  temp_grid[[j]] <- plot
  j = j+1
    if (j==7 | i==last) {
    do.call(grid.arrange,c(temp_grid,nrow=3,ncol=2,newpage=TRUE))
    temp_grid <- list()
    j=1
    }
  }
  dev.off()
  cat("\nPlotting finished. The file is exported under:\n",getwd())
}