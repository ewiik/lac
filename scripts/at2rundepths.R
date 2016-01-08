## correcting and updating running depths for at2 macrofossils
## source run depth file: dropbox: AT2 running depths_version2.xlsx
##  AT2 running depths_2014_version1.xlsx

## read in files with run depths
oldrun <- read.csv("data/oldrun.csv")
newrun <- read.csv("data/newrun.csv")
prevrun <- read.csv("data/prevrun.csv") # the depths that previously were on Dropbox

## read in macrofossil files
if(!file.exists("data/private/at2macroallani.csv")) {
  stop("get macrofossil files from Dropbox") }
ani <- read.csv("data/private/at2macroallani.csv")
ani <- ani[,-which(names(ani) == 'X')]

if(!file.exists("data/private/at2macroallterr.csv")) {
  stop("get macrofossil files from Dropbox") }
terr <- read.csv("data/private/at2macroallterr.csv")
terr <- terr[,-which(names(terr) == 'X')]
## FIXME: correct name of 'rabbit'

if(!file.exists("data/private/at2macroallplant.csv")) {
  stop("get macrofossil files from Dropbox") }
plant <- read.csv("data/private/at2macroallplant.csv")
plant <- plant[,-which(names(plant) == 'X')]

if(!file.exists("data/private/at2_aquatics_pcascores.csv")) {
  stop("get macrofossil files from Dropbox") 
}
aquapca <- read.csv("data/private/at2_aquatics_pcascores.csv")
aquapca <- aquapca[,-which(names(aquapca) == 'X')]

if(!file.exists("data/private/at2_terrestrials_pcascores.csv")) {
  stop("get macrofossil files from Dropbox") 
}
terrpca <- read.csv("data/private/at2_terrestrials_pcascores.csv")
terrpca <- terrpca[,-which(names(terrpca) == 'X')]

## reassign drive and depth denotations in prevrun to allow merging
prevrun$drive[prevrun$drive == 1] <- 'AT2-R1-1'
prevrun$drive[prevrun$drive == 2] <- 'AT2-R1-2'
prevrun$drive[prevrun$drive == 4] <- 'AT2-R1-4'

names(prevrun)[which(names(prevrun) == 'depth')] <- 'prevdepth'
prevrun$prevtop <- prevrun$prevdepth # create topdepth (these middepths)
prevrun$prevtop[prevrun$prevtop == 0.5] <- 1
prevrun$prevtop <- prevrun$prevtop - 1

## merge previous with newest
combo <- merge(newrun[,-which(names(newrun) == 'newrunbottom' | names(newrun) == 'newbottom')],
               prevrun[,-which(names(prevrun) == 'prevdepth')], 
               by.x = c("newdrive", "newtop"), by.y = c("drive","prevtop"), all.x = TRUE,
               all.y = TRUE)

## which depths have I analysed that are not in the newrun excel sheet?
reassign <- which(is.na(combo$newruntop))

## assign correct rundepth to them
combo$newruntop[reassign] <- c(125,127)

## remove rows that were not in macrofossil analysis
macros <- na.omit(combo)

## put updated depth into macro tables
anicorr <- cbind(macros$newruntop, ani)
names(anicorr)[1] <- "rundepthtop"

terrcorr <- cbind(macros$newruntop, terr)
names(terrcorr)[1] <- "rundepthtop"

plantcorr <- cbind(macros$newruntop, plant)
names(plantcorr)[1] <- "rundepthtop"

aquapcacorr <- cbind(macros$newruntop, aquapca[,-which(names(aquapca) == "depth")])

terrpcacorr <- cbind(macros$newruntop, terrpca[,-which(names(terrpca) == "depth")])

## write csvs
write.csv(anicorr, "data/private/at2macroallanicorr.csv", row.names = FALSE)
write.csv(terrcorr, "data/private/at2macroallterrcorr.csv", row.names = FALSE)
write.csv(plantcorr, "data/private/at2macroallplantcorr.csv", row.names = FALSE)

write.csv(aquapcacorr, "data/private/at2_aquatics_pcascores_corr.csv", row.names = FALSE)
write.csv(terrpcacorr, "data/private/at2_terrestrials_pcascores_corr.csv", row.names = FALSE)

