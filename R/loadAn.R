loadAn = function(file) {
  the_xs = eval(parse(text=load(file)))
  if (class(the_xs) != "xsAnnotate") { stop(paste(sep="","The xsAnnotate (", file,") was not loaded correctly.")) }
  return(the_xs)
}