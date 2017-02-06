plot.knots = function(.knot, Nmacha, scales = c(1,1)) {
  df = Nmacha$cc[Nmacha$cc_knot,,on="cc"][Nmacha$knot[knot%in%.knot],,on="knot", nomatch=0] %>% as.data.frame
  xlim = c(min(df$mz)-15, max(df$mz)+15)
  
  ylimmeanr = c(mean(df$meanr)-scales[1]*1, mean(df$meanr) + scales[1]*1)
  ylimrt = c(mean(df$rt)-scales[2]*1, mean(df$rt) + scales[2]*1)
  
  spec = ggplot(df) + geom_segment(aes(x = mz, xend = mz, y = i, yend = 0, colour=factor(knot))) + ggtitle(paste(.knot)) + theme_nate() + theme(legend.position="none") + geom_segment(aes(y = 0, yend = 0, x = mean(range(df$mz)) - 5, xend = mean(range(df$mz)) + 5)) + coord_cartesian(xlim = xlim)
  mrspec = ggplot(df) + geom_point(aes(x = mz, y = meanr, colour=factor(knot))) + theme_nate() + theme(legend.position="none") + coord_cartesian(xlim = xlim, ylim=ylimmeanr)
  rtspec = ggplot(df) + geom_point(aes(x = mz, y = rt, colour=factor(knot))) + theme_nate() + theme(legend.position="none") + coord_cartesian(xlim = xlim, ylim = ylimrt)
  
  psdf = Nmacha$knot[knot%in%.knot] %>% as.data.frame
  distspec = ggplot(psdf) + geom_point(aes(x = rt, y = meanr, colour = factor(knot), shape = factor(z))) + theme_nate() + theme(legend.position="none") + coord_cartesian(xlim = ylimrt, ylim = ylimmeanr)
  
  gridExtra::grid.arrange(mrspec, spec, rtspec, distspec, nrow = 4)
}

plot.quipu = function(.quipu, Nmacha) {
  plot.knots(Nmacha$quipu[Nmacha$knot_quipu[quipu %in% .quipu],,on="quipu"]$knot, Nmacha)
}