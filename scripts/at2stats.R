## Multivariate stats on AT2 Clustering methods originally from cuns macro methods, 
##    see multivariatethesis and multivariate2 in R/

## load necessary packages
# library("mvpart") not available for R 3.2.3
library("rioja")

## read in old strati objects originally from R scripts in R/
opts <- read.csv("data/private/at2macrowaalt.csv") # optima along strati plot
macrowa <- read.csv("data/private/at2macrowaset.csv") # taxa for strati plots
macrowa <- macrowa[,-which(names(macrowa) == 'X')]
terrsub <- read.csv("data/private/at2macroterrsubset.csv") # subset of terrestrials
terrsub <- terrsub[,-which(names(terrsub) == 'X')]

## read in basic ones put on Dropbox
ani <- read.csv("data/private/at2macroallanicorr.csv")
terr <- read.csv("data/private/at2macroallterrcorr.csv")
plant <- read.csv("data/private/at2macroallplantcorr.csv")


## =================================================================================
## Stratigraphic zones for AT2 macrofossils: aquatics and terrestrials separately
## =================================================================================
## Aquatics: standardise data but no centring since using Bray and want to preserve
##    0 abundances
greppedp <- grep('rundepth', names(plant))
greppeda <- grep('rundepth', names(ani))
aqua <- cbind(ani[,-greppeda], plant[,-greppedp])
aquascale <- scale(aqua, center = FALSE, scale = TRUE) 
aquascale <- as.data.frame(aquascale)
aquadist <- vegdist(aquascale, method = "bray", binary = FALSE)
aquazones <- chclust(aquadist)
plot(aquazones)
aquastick <- bstick(aquazones, 8) 
## see stratigraphic-prac.pdf; I think that the most
##    parsimonious point is where the red line first goes above the black line going
##    from left to right on the plot.... maybe...
aquastick ## 5 groups for good fit?? ...1-15, 18-37, 38-55, 56-75, (76-83, 84-104)

## clustering for terrestrials; used terrsub originally, to reduce redundancy of
##    different parts, same plant, when showing exactly same pattern: ALSO:
##    In terrsub, some occasionals have been lumped together taxonomically, e.g.
##    Draba + Cochlearie --> Brassicaceae
##    Sagina + Minuartia --> Caryophyllaceae
terrscale <- scale(terrsub, center = FALSE, scale = TRUE) 
terrscale <- as.data.frame(terrscale)
rowSums(terrscale) ## 91 and 92 yield 0 which is why vegdist no likey
terrdist <- vegdist(terrscale[-c(92,91),], method = "bray", binary = FALSE)
terrzones <- chclust(terrdist)
plot(terrzones)
terrstick <- bstick(terrzones, 6) ## see stratigraphic-prac.pdf
terrstick ## even 3 groups for good fit, but taking more will yield ..1-4, 5-17, 18-33, 
##  34-57, 58-78, 79-104 (need to convert to strat depths)

## let's try clustering for all other terrestrials too.... all aquatics in the start
##    of macrowa, so will remove aqua cols --> 30 vs 21
terrsub2 <- macrowa[,-c(1:15)]
terrscale2 <- scale(terrsub2, center = FALSE, scale = TRUE) 
terrscale2 <- as.data.frame(terrscale2)
remove <- which(rowSums(terrscale2) == 0)
terrdist2 <- vegdist(terrscale2[-remove,], method = "bray", binary = FALSE)
terrzones2 <- chclust(terrdist2)
plot(terrzones2)
terrstick2 <- bstick(terrzones2, 6) ## see stratigraphic-prac.pdf
terrstick2 ## 4 groups for good fit
## depends whether or not I believe that redundancy reveals more than just redundancy, e.g.
##    in form of catchment stability etc. 
## FIXME: WHY did I not include only 3 or 4 in the summary diagram that I sent to everyone?
## WHY are there 5 (with a and b)


## NJA wants PCA scores for the data, separately for terrestrial and aquatic.... 
##  though the one dominant gradient thing argues against this.... but will do.
## centring and standardising here should be ok, as I will not do non-symmetric distance measures 
##  (e.g. Bray Curtis) which require absences to be 0s. 
##"Who advised on the way to treat macrofossil data prior to analysis? It should be standard to 
##  centre and standardize Mfossil data, or maybe normalize it as you did the pigment data –
##  I would say this must be done as the ## units of the different remains are not comparable – 
##  centering and standardizing is a good way for getting around that. If you want to use Bray Curtis 
##  it is important to return former absences (they will now be a minus ## values as the mean has been 
##  subtracted prior to dividing by the standard deviation) to zeros prior to the nMDS though – as bc 
##  ignores zeros. What you did might be ok, but the nMDS will tend to more affected by things 
## that leave more numerous remains – i.e. chara and nympheae trichosclerids..  Log (x+1) might be enough  - 
##  but I think log10 is a stronger transformation that ln so might be better in this case"
## however this from multivariate2: "Gav recommends log transformation for abundance data for PCA"

## now for pcas
terrsubstand <- scale(terrsub, scale = TRUE, center = TRUE) ## should've logged?? 
terrsub.pca <- prcomp(terrsub, tol = sqrt(.Machine$double.eps), scale = FALSE, center = FALSE)
## so this could have centered and standardised for me...
scores(terrsub.pca)[,c(1:2)] ## PCs 1 and 2
ordiplot(terrsub.pca, type = "t", display = "sites")
points(as.numeric(scores(terrsub.pca)[,1]), as.numeric(scores(terrsub.pca)[,2]), pch = 16, 
       col = rep(c("red","blue","green","black"),times = c(25,25,25,29)), cex = 1) 
## horseshoe 
lines(as.numeric(scores(terrsubset.pca)[,1]), as.numeric(scores(terrsubset.pca)[,2])) 
biplot(terrsubset.pca) ## plots points on according to top and right scales, but alters scores 
##  to make plotting nicer, without letting me access those points... aaaghhhh how to plot onto this? 
##  hence used ordiplot
terrsubset.pca$sdev

## and for the logged alternative..:
terrsubstandlog <- scale(log(terrsub + 1), scale = TRUE, center = TRUE) 
terrsublog.pca <- prcomp(terrsubstandlog, tol = sqrt(.Machine$double.eps), scale = FALSE, center = FALSE)
scores(terrsublog.pca)[,c(1:2)] ## PCs 1 and 2
terrscores <- data.frame(scores(terrsublog.pca)[,c(1:2)])
terrscores$depth <- terr[,1]

ordiplot(terrsublog.pca, type = "t", display = "sites")
terrsub.pca$sdev ## vs other. not massively different
terrsubeigen <- terrsub.pca$sdev^2 
terrsubvarexp <- vector()
for (i in 1:21) {terrsubvarexp[i] <- terrsubeigen[i]/sum(terrsubeigen) } 

##let's check level of community turnover:
## All row sums must be >0 in the community matrix: remove empty sites.
remove <- which(rowSums(terrsub) == 0)
at2terrdca <- decorana(log(terrsub[-remove,] + 1))
#DCA 1 and 2 > 4sd.
decorana(log(aqua + 1))
#under 3 here.

##so let's do ca on terrestrials
at2terrca <- cca(log(terrsub[-remove,] + 1))
#Eigenvalues for unconstrained axes:
#   CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8 
#0.5369 0.4380 0.3432 0.3344 0.2265 0.2092 0.1726 0.1557
plot(at2terrca) # shows how wildly different some are... this may be an artefact of the rare species, 
##  as most of my data here are quite low in number i.e. occur in the data set rarely (low abundances, 
##  but frequent occurrences, are ok). I could downweight!
at2terrcadown <- cca(downweight((terrsub[-remove,] + 1)))
plot(at2terrcadown)
ordiplot(at2terrcadown, type = "t", display = "sites") # looks a lot tighter. and doesn't seem to have an 
##  arch effect in that they are not arranged in an arch dictated by time, and it is not hollow...
scores(at2terrcadown) # but still confused, if e.g. one axis is a distortion of the other, how can the 
##  scores be valid for any interpretation with axis 1 or 2 vs stratigraphic depth?


## now for aquatics....
aquastand <- scale(aqua, scale = TRUE, center = TRUE) 
## took out redundant sessoblasts, and chydorids which show nothing
aqua.pca <- prcomp(aquastand, tol = sqrt(.Machine$double.eps), scale = FALSE, center = FALSE) 
scores(aqua.pca)[,c(1:2)] ## PCs 1 and 2
ordiplot(aqua.pca, type = "t", display = "sites")
points(as.numeric(scores(aqua.pca)[,1]), as.numeric(scores(aqua.pca)[,2]), pch = 16, 
       col = rep(c("red","blue","green","black"),times = c(25,25,25,29)), cex = 1) 
## top and bottom again close to each other...
lines(as.numeric(scores(aqua.pca)[,1]), as.numeric(scores(aqua.pca)[,2])) 
biplot(aqua.pca) ## plots points on according to top and right scales, but alters scores to make 
##  plotting nicer, without letting me access those points... aaaghhhh how to plot onto this?
##  hence used ordiplot

#and logged alternative:
aqualogstand <- scale(log(aqua + 1), scale = TRUE, center = TRUE) 
aqualog.pca <- prcomp(aqualogstand, tol = sqrt(.Machine$double.eps), scale = FALSE, center = FALSE) 
scores(aqualog.pca)[,c(1:2)] ## PCs 1 and 2
aquascores <- data.frame(scores(aqualog.pca)[,c(1:2)])
aquascores$depth <- terr[,1]
ordiplot(aqualog.pca, type = "t", display = "sites")
points(as.numeric(scores(aqualog.pca)[,1]), as.numeric(scores(aqualog.pca)[,2]), pch = 16, 
       col = rep(c("red","blue","green","black"),times = c(25,25,25,29)), cex = 1) 
## top and bottom again close to each other...
lines(as.numeric(scores(aqualog.pca)[,1]), as.numeric(scores(aqualog.pca)[,2])) 
biplot(aqualog.pca) ## plots points on according to top and right scales, but alters scores to make 
##  plotting nicer, without letting me access those points... aaaghhhh how to plot onto this? 
##  hence used ordiplot
# there is still a horseshoe effect so not too confident about this, I think CA would be better!

