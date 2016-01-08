## initial plotting and looking at AT2 data

## load necessary packages
library("vegan")
library("analogue")

## read in run depth-corrected data files
if(!file.exists("data/private/at2macroallanicorr.csv")) {
  source("scripts/at2rundepths.R")
} # all csvs are created in same script so this should be enough of a check for the files
ani <- read.csv("data/private/at2macroallanicorr.csv")
terr <- read.csv("data/private/at2macroallterrcorr.csv")
plant <- read.csv("data/private/at2macroallplantcorr.csv")

## read in old strati objects from R scripts in R/
opts <- read.csv("data/private/at2macrowaalt.csv") # optima along strati plot
macrowa <- read.csv("data/private/at2macrowaset.csv") # taxa for strati plots
macrowa <- macrowa[,-which(names(macrowa) == 'X')]

## plot stratiplots
terrplot <- Stratiplot(terr[,-1], ylab = list(label = "Depth (cm)", cex = 1.2), 
                       xlab = list(label = expression("Counts as n(individuals)" ~ 100 ~ cm^{-3} ~ 
                                                        "wet sediment OR abundance score"), cex = 1.2), 
                       terr[,1], type = "h", pages = 1, rev = TRUE, sort = "wa",rev.sort = TRUE, 
                       strip = FALSE, topPad =6, lwd = 3, varTypes = "absolute",col = "black", 
                       absoluteSize = 0.5, drawLegend = TRUE)
plantplot <- Stratiplot(plant[,-1], ylab = list(label = "Depth (cm)", cex = 1.2), 
                        xlab = list(label = expression("Counts as n(individuals)" ~ 100 ~ cm^{-3} ~ 
                                                         "wet sediment OR abundance score"), cex = 1.2), 
                        plant[,1], type = "h", pages = 1, rev = TRUE, sort = "wa",rev.sort = TRUE, 
                        strip = FALSE, topPad =6, lwd = 3, varTypes = "absolute",col = "black", 
                        absoluteSize = 0.5, drawLegend = TRUE)
aniplot <- Stratiplot(ani[,-1], ylab = list(label = "Depth (cm)", cex = 1.2), 
                      xlab = list(label = expression("Counts as n(individuals)" ~ 100 ~ cm^{-3} ~ 
                                                       "wet sediment OR abundance score"), cex = 1.2), 
                      ani[,1], type = "h", pages = 1, rev = TRUE, sort = "wa",rev.sort = TRUE, 
                      strip = FALSE, topPad =6, lwd = 3, varTypes = "absolute",col = "black", 
                      absoluteSize = 0.5, drawLegend = TRUE)

## FIXME: check if numeric selections correct!
at2sumtaxa <- cbind(subset(at2macroallterr, select = c(7,8,9,11,12,15,19,22,24,25,27,31)), 
                    subset(at2macroallani, select = c(3,6,7,9,10,11)), at2macroallplant)
at2sumplot <- Stratiplot(at2sumtaxa, ylab = list(label = "Depth (cm)", cex = 1.2), 
                         xlab = list(label = expression("Counts as n(individuals)" ~ 100 ~ cm^{-3} ~ 
                                                          "wet sediment OR abundance score"), cex = 1.2), 
                         at2sumtaxa[,1], type = "h", pages = 1, rev = TRUE, sort = "wa",rev.sort = TRUE, 
                         strip = FALSE, topPad =6, lwd = 3, varTypes = "absolute",col = "black", 
                         absoluteSize = 0.5, drawLegend = TRUE)

at2terrstratiplotwa <- Stratiplot(at2macrowaset[,c(18:42,45)], at2sumtaxa[,1], type = "h", ylab = "Depth (cm)", 
                                  xlab = expression("Counts as n(individuals)" ~ 100 ~ cm^{-3} ~ 
                                                      "wet sediment"), pages = 1, rev = TRUE, sort = "var", 
                                  svar = at2macrowaalt$opt[c(18:42,45)], rev.sort = TRUE, strip = FALSE, 
                                  topPad= 11, varTypes = "absolute", absoluteSize = 0.5, drawLegend = TRUE, 
                                  col = "black", zones = c(4.5,30.5,62.5,111.5,153.5), 
                                  zoneNames = c(1,2,3,4,"5a","5b"), lwd = 4) # zones from terr cluster analysis
at2aquastratiplotwa <- Stratiplot(at2macrowaset[,-c(18:45)], at2sumtaxa[,1], type = "h", ylab = "Depth (cm)", 
                                  xlab = expression("Counts as n(individuals)" ~ 100 ~ cm^{-3} ~ 
                                                      "wet sediment"), pages = 1, rev = TRUE, sort = "var", 
                                  svar = at2macrowaalt$opt[-c(18:45)], rev.sort = TRUE, strip = FALSE, 
                                  topPad= 11, varTypes = "absolute", absoluteSize = 0.5, drawLegend = TRUE, 
                                  col = "black", zones = c(26.5,72.5,107.5,147.5), zoneNames = c(1,2,3,4,5), 
                                  lwd = 4) # zones from terr cluster analysis

## create pdfs
pdf("data/private/at2macros.pdf",width = 15)
at2macroallterrplot
at2macroallaniplot
at2macroallplantplot
at2sumplot
dev.off()

pdf("at2macrostratiplotwa.pdf", width = 15)
at2macrostratiplotwa
dev.off()
png("at2macrostratiplotwa.png", width = 1500) # appears fuzzy =(
at2macrostratiplotwa
dev.off()

setEPS()
postscript("at2terrstratiplotwa.eps", width = 10, height = 7)
at2terrstratiplotwa
dev.off()

setEPS()
postscript("at2aquastratiplotwa.eps", width = 8, height = 7)
at2aquastratiplotwa[-1]
dev.off()



at2ord <- xyplot(-depth ~ values | ind, data = at2pcastack, type = "o", layout = c(5,1), scales = list(x = list(relation = "free")), col = "black")
## !!! to be made, csv files with the scores for NJA
setEPS()
postscript("at2ord.eps", width = 7, height = 5)
at2ord
dev.off()

## loi diagram with AT1 and AT2 LOIs for radiocarbon justification...
gloistack <- at2loi
gloistack$core <- rep("AT2")
at1loi <- read.csv("at1loirundepths.csv", header = TRUE)
at1loi$core <- rep("AT1")
gloistack <- rbind(gloistack, at1loi)
at1antloi <- read.csv("at1antloirundepths.csv", header = TRUE)
gloistack <- rbind(gloistack, at1antloi)
prettydepths <- seq(0, 280, by = 10)
gloiplot <- xyplot(-rundepth ~ loi | core, data = gloistack, main = "", type = "o", layout = c(3,1), xlab = "% organic", ylab = "Depth (cm)", scales = list(relation = "same",y = list(at = -prettydepths,labels = prettydepths)), col = "black")
trellis.focus(name = "panel", column = 2, row = 1)
panel.text(x = c(8, rep(5, times = 5)), y = c(0,-15,-107, -155, -212, -263), labels = c("Age (ky BP)",".24", "3.6", "5.8", "7.5", "9.5"), cex = 1)
trellis.unfocus()
trellis.focus(name = "panel", column = 3, row = 1)
panel.points(x = c(rep(5, times = 7), 13), y = c(-46,-67.5,-81,-107,-134,-177,-202,-209), pch = 8, col = "black")
panel.text(x = 10, y = -280, labels = "* = Submitted")
trellis.unfocus()
setEPS()
postscript("greenlanddates.eps")
#pdf("greenlanddates.pdf", onefile = FALSE)
gloiplot
trellis.focus(name = "panel", column = 2, row = 1)
panel.text(x = c(8, rep(5, times = 5)), y = c(0,-15,-107, -155, -212, -263), labels = c("Age (ky BP)",".24", "3.6", "5.8", "7.5", "9.5"), cex = 1)
trellis.unfocus()
trellis.focus(name = "panel", column = 3, row = 1)
panel.points(x = c(rep(5, times = 7), 13, 2), y = c(-46,-67.5,-81,-107,-134,-177,-202,-209,-280), pch = 8, col = "black")
panel.text(x = 15, y = -280, labels = " = Submitted")
trellis.unfocus()
trellis.focus(name = "panel", column = 1, row = 1)
panel.points(x = 3, y = -c(35, 75, 130, 170, 220, 280), pch = 3, col = "black")
panel.text(x = 15, y = -280, labels = " = Requested")
trellis.unfocus()
dev.off()
##prelim dates for at2: 
prelimdepths <- c(81,134,202,209)
prelimdates <- c(6.3, 7.8, 8.8, 9.4)
trellis.focus(name = "panel", column = 3, row = 1)
panel.text(x = 35, y = -c(prelimdepths), labels = c(prelimdates))
trellis.unfocus()
# made a prelim eps

##xyplot(-depth~values|ind, data = isostack, main = "Core data from pilot study", ylab = "depth",xlab = "Isotope values as ppm ",scales = list(relation = "sliced",y = list(at = pretty(-prettydepths),labels = rev(format(pretty(loihaw2$depth))))), type = "o", strip=strip.custom(strip.levels = c(TRUE, TRUE),strip.names=c(FALSE, FALSE), factor.levels= c(expression(paste("carbonate"~delta^{13}~"C")),expression(paste("carbonate"~delta^{18}~"O")),expression(paste("bulk organic"~delta^{13}~"C")),"C/N"), par.strip.text = list(cex = 0.7)), as.table = TRUE,layout = c(4,1))
expression(paste("total phosphorus, "~mu*gL^{-1}))