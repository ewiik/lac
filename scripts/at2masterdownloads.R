## read in all supporting data for AT2 and get it organised
## using master file for pigs and C/N stuff....
## FIXME: no actual diatom counts in Dropbox????

## read in files
master <- read.csv("data/private/AT2_MasterSpreadsheet_15-12-15.csv") # rundepth is topdepth
cladorel <- read.csv("data/private/AT2-Cladocera-counts.csv") # this is file last modified 8th Feb 2015
cladoraw <- read.csv("data/private/AT2-Cladocera-counts-raw.csv") # this is file last modified 8th Feb 2015
cladoraw[is.na(cladoraw)] <- 0 # replace NA with 0, since these are true 0s

chiroraw <- read.csv("data/private/AT2-chiro-counts-raw.csv") # this is the file last modified
## FIXME: check with Maarten that this one (sheet "cleaned") is actually raw data
chiroraw[is.na(chiroraw)] <- 0 # replace NA with 0, since these are true 0s

plant <- read.csv("data/private/at2macroallplantcorr.csv")

## create pigments
pigseq <- grep("Phaeo|ytin.a", names(master))
pigs <- cbind(master$Running.Depth, master[,pigseq[1]:pigseq[2]])
names(pigs)[1] <- "rundepthtop"

## create chemistry; X. denotes %.
geoseq <- grep("BioS|C.N", names(master))
geos <- cbind(master$Running.Depth, master$LOI, master[,geoseq[1]:geoseq[2]])
names(geos)[1:2] <- c("rundepthtop", "LOI")

## correct rundepth for clados; use rundepth for macros since same sample material used
## --> know that last sample same since clados also terminate at 200something with larger gap
##    in the last two samples
take <- nrow(cladorel)
max <- nrow(plant)
taken <- plant[(max-take + 1):max,1]

cladorel$Depth <- taken
names(cladorel)[names(cladorel) == "Depth"] <- "rundepthtop"

cladoraw <- cbind(taken, cladoraw)
names(cladoraw)[1] <- "rundepthtop"

chiroraw <- cbind(taken, chiroraw)
names(chiroraw)[1] <- "rundepthtop"

## create initial Stratiplots for initial discussion
pdf("data/private/allplots.pdf", width = 15, onefile = TRUE)
Stratiplot(geos[,-1], geos[,1], type = "h", varTypes = "absolute", col = "black")
Stratiplot(pigs[,-1], pigs[,1], type = "h", varTypes = "absolute", col = "black")
Stratiplot(cladoraw[,-1], cladoraw[,1], type = "h", varTypes = "absolute", col = "black")
Stratiplot(chiroraw[,-1], chiroraw[,1], type = "h", varTypes = "absolute", col = "black")
dev.off()