loadXs = function(file) {
  the_xs = eval(parse(text=load(file)))
  if (class(the_xs) != "xcmsSet") { stop(paste(sep="","The xcmsSet (", file,") was not loaded correctly.")) }
  return(the_xs)
}