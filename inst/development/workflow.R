library(data.table)
library(magrittr)
library(igraph)
library(matrixStats)
library(ggplot2)

source("R/credential_knots.R")
source("R/find_knots.R")
source("R/helpers.R")
source("R/plot.R")

dt = readRDS("../macha/inst/development_nmacha/cache/2NM133A_200t1900_13C.nmacha.wg.rds")

features = dt$cc

knots = findknots(features, .zs=1:2, ppmwid=4, rtwid = 1, cd = 13.00335-12)

credentials = credentialknots(knots$knot, ppmwid = 9, rtwid = 1, mpc = c(12, 120), ratio = 1/1, ratio.lim = 0.1, maxnmer = 4, cd = 13.00335-12)

#For toy data

#Find a good credentialed peak

credentials$quipu[minsupport > 4]

.quipu = 1

df = features[knots$cc_knot[credentials$knot_quipu[quipu == 1],,on="knot"],,on="cc"]

ggplot(df) + geom_segment(aes(x = mz, xend = mz, yend = 0, y = i))

### Redo with our data
features = dt$cc[cc %in% .cc][,.(cc, mz, rt, i)]
saveRDS(features, "features.rds")

features = readRDS(system.file("data", "features.rds", package="credential"))
features


knots = findknots(features, .zs=1:2, ppmwid=4, rtwid = 1, cd = 13.00335-12)

credentials = credentialknots(knots$knot, ppmwid = 9, rtwid = 1, mpc = c(12, 120), ratio = 1/1, ratio.lim = 0.1, maxnmer = 4, cd = 13.00335-12)



### Summarize
features = dt$cc

knots = findknots(features, .zs=1:2, ppmwid=4, rtwid = 1, cd = 13.00335-12)

credentials = credentialknots(knots$knot, ppmwid = 9, rtwid = 1, mpc = c(12, 120), ratio = 1/1, ratio.lim = 0.1, maxnmer = 4, cd = 13.00335-12)



summarize = function(feautres, knots, credentials) {
  
  
  cat("Number of peaks per knot")
  knots$cc_knot$knot %>% table %>% table
  
  cat("Charge states of isotopic knots")
  knots$knot$z %>% table
  
  cat("Ppm error within knots")
  hist(features[knots$cc_knot,,on="cc"][,mean(abs(diff(mz %% 1.003355))/mean(mz)*1E6),by="knot"]$V1, breaks = 100)
  
  
  cat("Ppm error between knots")
  hist(knots$knot[credentials$knot_quipu,,on="knot"][!is.na(quipu),mean(abs(diff(meanr))/mean(meanmz)*1E6),by="quipu"]$V1, breaks = 100)
  
  
  cat("Number of knots per credential")
  credentials$knot_quipu$quipu %>% table %>% table
  
  
  cat("Credentialed peaks with iso support")
  features[knots$knot[knots$cc_knot[credentials$knot_quipu[credentials$quipu[minsupport > 1],,on="quipu"],,on="knot"],,on="knot"],,on="cc"][,.SD[mz == min(mz)],by="quipu"][,.(cc, mz, rt, i, z, knot, quipu)]
  
  cat("Credentialed peaks without iso support")
  features[knots$knot[knots$cc_knot[credentials$knot_quipu[credentials$quipu[minsupport == 1],,on="quipu"],,on="knot"],,on="knot"],,on="cc"][,.SD[mz == min(mz)],by="quipu"][,.(cc, mz, rt, i, z, knot, quipu)]
  
  
  
  }

