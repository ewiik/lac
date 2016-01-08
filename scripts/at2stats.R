## Multivariate stats on AT2 Clustering methods originally from cuns macro methods, 
##    see multivariatethesis and multivariate2 in R/

## load necessary packages
library("mvpart")
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


## standardise data
macrowascale <- scale(macrowa, center = FALSE, scale = TRUE) # didn't center, as for Bray, 
#   I want to keep 0 abundances as 0.
macrowascale <- as.data.frame(macrowascale)
at2macrodist <- vegdist(macrowascale, method = "bray", binary = FALSE)
at2zones <- chclust(at2macrodist)
plot(at2zones) 
at2stick <- bstick(at2zones, 6) ## see stratigraphic-prac.pdf

## clustering for aquatics 
aqua <- cbind(ani[,-c(1,5,6)], plant)
aquascale <- scale(aqua, center = FALSE, scale = TRUE) ## didn't center, as for Bray,
##  I want to keep 0 abundances as 0.
aquascale <- as.data.frame(aquascale)
aquadist <- vegdist(aquascale, method = "bray", binary = FALSE)
aquazones <- chclust(aquadist)
plot(aquazones)
aquastick <- bstick(aquazones, 6) ## see stratigraphic-prac.pdf
aquastick ## 5 groups for good fit...1-15, 18-38, 39-55, 56-75, 76-104, OR 76-86, 87-104, 
##  emerge (need to convert to strat depths)
## FIXME: check if with new rerun, I get real depths not rownumbers. i.e. first column shit

## ================================================================GOT HERE
## clustering for terrestrials
terrscale <- scale(terrsub, center = FALSE, scale = TRUE) 
## didn't center, as for Bray, I want to keep 0 abundances as 0.
terrscale <- as.data.frame(terrscale)
rowSums(terrscale) ## 91 and 92 yield 0 which is why vegdist no likey
terrdist <- vegdist(terrscale[-c(92,91),], method = "bray", binary = FALSE)
terrzones <- chclust(terrdist)
plot(terrzones)
terrstick <- bstick(terrzones, 6) ## see stratigraphic-prac.pdf
terrstick ## even 3 groups for good fit, but taking more will yield ..1-4, 5-17, 18-33, 
##  34-57, 58-78, 79-104 (need to convert to strat depths)

at2depths <- list(c(at2macroalltot$rundepth), c(1:104))
at2depths[[1]][c(1,4,17,33,57,78,104)]


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
at2macroterrsubset <- subset(at2macroallterr, select = c(1,3,4,5,6,8,10,11,12,15,17,18,19,20,21,23,
                                                         24,26,28,29,30,31,32)) 
## not Vaccinium in yet, awaiting ID
at2macroterrsubset[,3] <- at2macroterrsubset[,3] + at2macroterrsubset[,5] 
## draba E Brassicaceae, so put together with Cochlearia
names(at2macroterrsubset)[3] <- "Brassicaceae \n (seed)"
at2macroterrsubset <- at2macroterrsubset[,-5]
at2macroterrsubset[,11] <- at2macroterrsubset[,11] + at2macroterrsubset[,12]  
## put sagina and minuartia together into Caryophyllaceae
names(at2macroterrsubset)[11] <- "Caryophyllaceae \n (seed)"
at2macroterrsubset <- at2macroterrsubset[,-12]
write.csv(at2macroterrsubset, "at2macroterrsubset.csv")
at2macroterrsubset <- read.csv("at2macroterrsubset.csv", row.names = 1)
terrsubsetnames <- names(at2macroterrsubset)

## now for pcas
terrsubsetstand <- scale(at2macroterrsubset, scale = TRUE, center = TRUE) ## should've logged?? 
##  see above for clustering separately
terrsubset.pca <- prcomp(terrsubsetstand, tol = sqrt(.Machine$double.eps), scale = FALSE, center = FALSE)
## so this could have centered and standardised for me...
scores(terrsubset.pca)[,c(1:2)] ## PCs 1 and 2
ordiplot(terrsubset.pca, type = "t", display = "sites")
points(as.numeric(scores(terrsubset.pca)[,1]), as.numeric(scores(terrsubset.pca)[,2]), pch = 16, 
       col = rep(c("red","blue","green","black"),times = c(25,25,25,29)), cex = 1) 
## shows how lovely and complete the horseshoe is....Complete arse!!
lines(as.numeric(scores(terrsubset.pca)[,1]), as.numeric(scores(terrsubset.pca)[,2])) 
biplot(terrsubset.pca) ## plots points on according to top and right scales, but alters scores 
##  to make plotting nicer, without letting me access those points... aaaghhhh how to plot onto this? 
##  hence used ordiplot
terrsubset.pca$sdev

## and for the logged alternative..:
terrsubsetstandlog <- scale(log(at2macroterrsubset + 1), scale = TRUE, center = TRUE) 
terrsubsetlog.pca <- prcomp(terrsubsetstandlog, tol = sqrt(.Machine$double.eps), scale = FALSE, center = FALSE)
scores(terrsubsetlog.pca)[,c(1:2)] ## PCs 1 and 2
terrscores <- data.frame(scores(terrsubsetlog.pca)[,c(1:2)])
terrscores$depth <- at2macroalltot$rundepth
write.csv(terrscores, "at2_terrestrials_pcascores.csv")

ordiplot(terrsubsetlog.pca, type = "t", display = "sites")
terrsubset.pca$sdev ## vs other. not massively different
terrsubseteigen <- terrsubset.pca$sdev^2 
terrsubsetvarexp <- vector()
for (i in 1:21) {terrsubsetvarexp[i] <- terrsubseteigen[i]/sum(terrsubseteigen) } 

##let's check level of community turnover:
decorana(log(at2macroterrsubset + 1)) ## All row sums must be >0 in the community matrix: remove empty sites.
at2terrdca <- decorana(log(at2macroterrsubset[-c(91,92),] + 1))
#DCA 1 and 2 > 4sd.
aquastandsubset <- cbind(at2macroallani[,-c(1,5,6)], at2macroallplant)
decorana(log(aquastandsubset + 1))
#under 3 here.

##so let's do ca on terrestrials
at2terrca <- cca(log(at2macroterrsubset[-c(91,92),] + 1))
#Eigenvalues for unconstrained axes:
#   CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8 
#0.5369 0.4380 0.3432 0.3344 0.2265 0.2092 0.1726 0.1557
plot(at2terrca) # shows how wildly different some are... this may be an artefact of the rare species, 
##  as most of my data here are quite low in number i.e. occur in the data set rarely (low abundances, 
##  but frequent occurrences, are ok). I could downweight!
at2terrcadown <- cca(downweight((at2macroterrsubset[-c(91,92),] + 1)))
plot(at2terrcadown)
ordiplot(at2terrcadown, type = "t", display = "sites") # looks a lot tighter. and doesn't seem to have an 
##  arch effect in that they are not arranged in an arch dictated by time, and it is not hollow...
scores(at2terrcadown) # but still confused, if e.g. one axis is a distortion of the other, how can the 
##  scores be valid for any interpretation with axis 1 or 2 vs stratigraphic depth?


## now for aquatics....
plantstand <- scale(at2macroallplant, scale = TRUE, center = TRUE)
plant.pca <- prcomp(plantstand, tol = sqrt(.Machine$double.eps), scale = FALSE, center = FALSE) 
## so this could have centered and standardised for me...
scores(plant.pca)[,c(1:2)] ## PCs 1 and 2
ordiplot(plant.pca, type = "t", display = "sites")
points(as.numeric(scores(plant.pca)[,1]), as.numeric(scores(plant.pca)[,2]), pch = 16, 
       col = rep(c("red","blue","green","black"),times = c(25,25,25,29)), cex = 1) 
## shows how lovely and complete the horseshoe is....Complete arse!!
lines(as.numeric(scores(plant.pca)[,1]), as.numeric(scores(plant.pca)[,2])) 
biplot(plant.pca) ## plots points on according to top and right scales, but alters scores to make 
##  plotting nicer, without letting me access those points... aaaghhhh how to plot onto this? 
##  hence used ordiplot

anistand <- scale(at2macroallani[,-c(1,5,6)], scale = TRUE, center = TRUE) 
## took out redundant sessoblasts,
##  and chydorids which show nothing
ani.pca <- prcomp(anistand, tol = sqrt(.Machine$double.eps), scale = FALSE, center = FALSE) 
## so this could have centered and standardised for me...
scores(ani.pca)[,c(1:2)] ## PCs 1 and 2
ordiplot(ani.pca, type = "t", display = "sites")
points(as.numeric(scores(ani.pca)[,1]), as.numeric(scores(ani.pca)[,2]), pch = 16, 
       col = rep(c("red","blue","green","black"),times = c(25,25,25,29)), cex = 1) 
## top and bottom again close to each other...
lines(as.numeric(scores(ani.pca)[,1]), as.numeric(scores(ani.pca)[,2])) 
biplot(ani.pca) ## plots points on according to top and right scales, but alters scores 
##  to make plotting nicer, without letting me access those points... aaaghhhh how to plot onto this? 
##  hence used ordiplot

combining ani and plant together:
  aquastand <- scale(cbind(at2macroallani[,-c(1,5,6)], at2macroallplant), scale = TRUE, center = TRUE) 
## took out redundant sessoblasts, and chydorids which show nothing
aqua.pca <- prcomp(aquastand, tol = sqrt(.Machine$double.eps), scale = FALSE, center = FALSE) 
## so this could have centered and standardised for me...
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
aqualogstand <- scale(log(cbind(at2macroallani[,-c(1,5,6)], at2macroallplant) +1), scale = TRUE, center = TRUE) 
## took out redundant sessoblasts, and chydorids which show nothing
aqualog.pca <- prcomp(aqualogstand, tol = sqrt(.Machine$double.eps), scale = FALSE, center = FALSE) 
## so this could have centered and standardised for me...
scores(aqualog.pca)[,c(1:2)] ## PCs 1 and 2
aquascores <- data.frame(scores(aqualog.pca)[,c(1:2)])
aquascores$depth <- at2macroalltot$rundepth
write.csv(aquascores, "at2_aquatics_pcascores.csv")
ordiplot(aqualog.pca, type = "t", display = "sites")
points(as.numeric(scores(aqualog.pca)[,1]), as.numeric(scores(aqualog.pca)[,2]), pch = 16, 
       col = rep(c("red","blue","green","black"),times = c(25,25,25,29)), cex = 1) 
## top and bottom again close to each other...
lines(as.numeric(scores(aqualog.pca)[,1]), as.numeric(scores(aqualog.pca)[,2])) 
biplot(aqualog.pca) ## plots points on according to top and right scales, but alters scores to make 
##  plotting nicer, without letting me access those points... aaaghhhh how to plot onto this? 
##  hence used ordiplot
# there is still a horseshoe effect so not too confident about this, I think CA would be better!

## LOI data to be included:
at2loi <- read.csv("at2loirundepths.csv", header = TRUE)
at2pcastack <- stack(as.data.frame(scores(at2terrcadown)[[2]]))
at2pcastack$depth <- rep(as.numeric(at2macroalltot$rundepth[-c(91,92)]))
at2pcastack$ind <- as.character(at2pcastack$ind)
at2pcastack$ind[at2pcastack$ind == "CA1"] <- "terrCA1"
at2pcastack$ind[at2pcastack$ind == "CA2"] <- "terrCA2"
at2pcastack <- rbind(at2pcastack, data.frame(values = at2loi$loi[1:166], ind = rep("loi"), 
                                             depth = at2loi$rundepth[1:166]))
at2pcastackaqua <- stack(data.frame(aquaPC1 = scores(aqualog.pca)[,1], aquaPC2 = scores(aqualog.pca)[,2]))
at2pcastackaqua$depth <- rep(as.numeric(at2macroalltot$rundepth))
at2pcastack <- rbind(at2pcastack, at2pcastackaqua)
