## create ages for AT2 based on Maarten's CLAM modeling
## files received from Erika 14.01.16 (Regina email)
## in case uncertainties check with Maarten what everything means

## read files in
ages <- read.csv("data/private/AT2_interpolated_ages.txt", sep = "")

## put in ages in AD, assuming that current form refers to distance from 1950
ages <- transform(ages, yearAD = 1950 - best)
ages <- transform(ages, plus = max95. - best)
ages <- transform(ages, minus = min95. - best)

## rename depth so that can match with other files
names(ages)[names(ages) == "depth"] <- "runtopdepth"

## save for further use
write.csv(ages, "data/private/AT2yearAD.csv")
