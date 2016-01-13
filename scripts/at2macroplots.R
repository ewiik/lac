## initial plotting and looking at AT2 data
## only the first depth (0) is a 1cm slice. All others are 2cm slices. Top depth given.

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
## these already have taxa lumped e.g. into Saxifragaceae, and unnecessaries *mostly*
##    removed. FIXME: "rabbit", "black ball", "bubblecushion"
## list of columns in existing summary plot:
## Sphagnum pres/abs, cf Nostoc, Oribatid, potobtleaf, simeph, tolyoo, plumrfloat, bconfseed,
##    hvulgseed, sphaersheath, bryoabun, tubicocoon, rhabdococoon, fredpipt, dapheph,
##    odigyseed, chameriumseed, sedumseed, poaceaeseed, saxiseed, gnivaseed, juncusseed, 
##    carexseed, luzulaseed, equispine, mycohyphae, coeno, ledumseed, vaccseed, harrileaf, 
##    salixseed cassiopeseed, phylloseed, vulilead, sherbleaf, salixbud, enigrseed, callunaseed,
##    enigrleaf, dintleaf, betulafruit i.e. 41... not leaf abun, rabbit, bubblecushion or blac ball
opts <- read.csv("data/private/at2macrowaalt.csv") # optima along strati plot
macrowa <- read.csv("data/private/at2macrowaset.csv") # taxa for strati plots
macrowa <- macrowa[,-which(names(macrowa) == 'X')] # remove row numbers

## plot stratiplots
terrplot <- Stratiplot(terr[,-c(1,3)], ylab = list(label = "Depth (cm)", cex = 1.2), 
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

## combine plots for terrestrials + aquatics, order by ecological group. Drop juniperinum since there
##    throughout and not showing anything
terrstratiplotwa <- Stratiplot(macrowa[,c(18:42,45)], terr[,1], type = "h", ylab = "Depth (cm)", 
                                  xlab = expression("Counts as n(individuals)" ~ 100 ~ cm^{-3} ~ 
                                                      "wet sediment"), pages = 1, rev = TRUE, sort = "var", 
                                  svar = opts$opt[c(18:42,45)], rev.sort = TRUE, strip = FALSE, 
                                  topPad= 11, varTypes = "absolute", absoluteSize = 0.5, drawLegend = TRUE, 
                                  col = "black", zones = c(4.5,30.5,62.5,111.5,153.5), 
                                  zoneNames = c(1,2,3,4,"5a","5b"), lwd = 4) # zones from terr cluster analysis
aquastratiplotwa <- Stratiplot(macrowa[,-c(16,18:45)], terr[,1], type = "h", ylab = "Depth (cm)", 
                                  xlab = expression("Counts as n(individuals)" ~ 100 ~ cm^{-3} ~ 
                                                      "wet sediment"), pages = 1, rev = TRUE, sort = "var", 
                                  svar = opts$opt[-c(16,18:45)], rev.sort = TRUE, strip = FALSE, 
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



at2ord <- xyplot(-depth ~ values | ind, data = at2pcastack, type = "o", layout = c(5,1), 
                 scales = list(x = list(relation = "free")), col = "black")
## !!! to be made, csv files with the scores for NJA
setEPS()
postscript("at2ord.eps", width = 7, height = 5)
at2ord
dev.off()
