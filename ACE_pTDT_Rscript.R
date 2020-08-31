## PTDT script
### Author: Varun Warrier


### Paper: https://www.nature.com/articles/ng.3863#methods

###Formulas
###PGSmidparent = (PGSfather + PGSmother)/2 ; 
###pTDT deviation is (PGSchild - PGSmidparent)/SD(PGSmidparent);
###tPDT = mean (pTDT deviation)/ (SD(pTDT deviation)/sqrt(N))



## STEP 1: Read files and merge

library(data.table)

#First, read the fam files
Mother = fread("~/SFARI/mothersdata.txt", header= T)
Father = read.table("~/SFARI/fathersdata.txt", header= T)
Cases = fread("~/SFARI/cases.txt", header= T)
siblings = fread("~/SFARI/siblingsdata.txt", header= T)

#Next, read the prs scores
PRS = fread("~/ALSPAC/PRSice2results/SFARI_autosomes_retro_prospcildtraumasum.all.score", header = TRUE)

#Create merged files
fatherpgs = merge(Father, PRS, by = "IID")
motherpgs = merge(Mother, PRS, by = "IID")
casespgs = merge(Cases, PRS, by = "IID")

##Step 2: Calculate mid-parent PGS and SD of midparent PGS

parentpgs = merge(motherpgs, fatherpgs, by = "FID.x")
parentpgs$midparent8 = (parentpgs$`1.000000.x` + parentpgs$`1.000000.y`)/2
parentpgs$midparent7 = (parentpgs$`0.750000.x` + parentpgs$`0.750000.y`)/2
parentpgs$midparent6 = (parentpgs$`0.500000.x` + parentpgs$`0.500000.y`)/2
parentpgs$midparent5 = (parentpgs$`0.250000.x` + parentpgs$`0.250000.y`)/2
parentpgs$midparent4 = (parentpgs$`0.100000.x` + parentpgs$`0.100000.y`)/2
parentpgs$midparent3 = (parentpgs$`0.010000.x` + parentpgs$`0.010000.y`)/2
parentpgs$midparent2 = (parentpgs$`0.001000.x` + parentpgs$`0.001000.y`)/2
parentpgs$midparent1 = (parentpgs$`0.000100.x` + parentpgs$`0.000100.y`)/2

triopgs = merge(parentpgs, casespgs, by = "FID.x")

Sd1 = sd(triopgs$midparent1)
Sd2 = sd(triopgs$midparent2)
Sd3 = sd(triopgs$midparent3)
Sd4 = sd(triopgs$midparent4)
Sd5 = sd(triopgs$midparent5)
Sd6 = sd(triopgs$midparent6)
Sd7 = sd(triopgs$midparent7)
Sd8 = sd(triopgs$midparent8)

## Step 3: Calculate pTDT deviation
triopgs$diff8 = (triopgs$`1.000000` - triopgs$midparent8)/Sd8
triopgs$diff7 = (triopgs$`0.750000` - triopgs$midparent7)/Sd7
triopgs$diff6 = (triopgs$`0.500000` - triopgs$midparent6)/Sd6
triopgs$diff5 = (triopgs$`0.250000` - triopgs$midparent5)/Sd5
triopgs$diff4 = (triopgs$`0.100000` - triopgs$midparent4)/Sd4
triopgs$diff3 = (triopgs$`0.010000` - triopgs$midparent3)/Sd3
triopgs$diff2 = (triopgs$`0.001000` - triopgs$midparent2)/Sd2
triopgs$diff1 = (triopgs$`0.000100` - triopgs$midparent1)/Sd1


## Step 4: Calculate the T score for pTDT deviation
N = sqrt(nrow(triopgs))

One = mean(triopgs$diff1)/(sd(triopgs$diff1)/N)
Two = mean(triopgs$diff2)/(sd(triopgs$diff2)/N)
Three = mean(triopgs$diff3)/(sd(triopgs$diff3)/N)
Four = mean(triopgs$diff4)/(sd(triopgs$diff4)/N)
Five = mean(triopgs$diff5)/(sd(triopgs$diff5)/N)
Six = mean(triopgs$diff6)/(sd(triopgs$diff6)/N)
Seven = mean(triopgs$diff7)/(sd(triopgs$diff7)/N)
Eight = mean(triopgs$diff8)/(sd(triopgs$diff8)/N)

One
Two
Three
Four
Five
Six
Seven
Eight



### In siblings

siblingpgs = merge(siblings, PRS, by = "IID")


controltriopgs = merge(parentpgs, siblingpgs, by = "FID.x")
controltriopgs = controltriopgs[!duplicated(controltriopgs[,c("FID.x")]),]

Sd1 = sd(controltriopgs$midparent1)
Sd2 = sd(controltriopgs$midparent2)
Sd3 = sd(controltriopgs$midparent3)
Sd4 = sd(controltriopgs$midparent4)
Sd5 = sd(controltriopgs$midparent5)
Sd6 = sd(controltriopgs$midparent6)
Sd7 = sd(controltriopgs$midparent7)
Sd8 = sd(controltriopgs$midparent8)

controltriopgs$diff8 = (controltriopgs$`1.000000` - controltriopgs$midparent8)/Sd8
controltriopgs$diff7 = (controltriopgs$`0.750000` - controltriopgs$midparent7)/Sd7
controltriopgs$diff6 = (controltriopgs$`0.500000` - controltriopgs$midparent6)/Sd6
controltriopgs$diff5 = (controltriopgs$`0.250000` - controltriopgs$midparent5)/Sd5
controltriopgs$diff4 = (controltriopgs$`0.100000` - controltriopgs$midparent4)/Sd4
controltriopgs$diff3 = (controltriopgs$`0.010000` - controltriopgs$midparent3)/Sd3
controltriopgs$diff2 = (controltriopgs$`0.001000` - controltriopgs$midparent2)/Sd2
controltriopgs$diff1 = (controltriopgs$`0.000100` - controltriopgs$midparent1)/Sd1


N = sqrt(nrow(controltriopgs))

One_Control = mean(controltriopgs$diff1)/(sd(controltriopgs$diff1)/N)
Two_Control = mean(controltriopgs$diff2)/(sd(controltriopgs$diff2)/N)
Three_Control = mean(controltriopgs$diff3)/(sd(controltriopgs$diff3)/N)
Four_Control = mean(controltriopgs$diff4)/(sd(controltriopgs$diff4)/N)
Five_Control = mean(controltriopgs$diff5)/(sd(controltriopgs$diff5)/N)
Six_Control = mean(controltriopgs$diff6)/(sd(controltriopgs$diff6)/N)
Seven_Control = mean(controltriopgs$diff7)/(sd(controltriopgs$diff7)/N)
Eight_Control = mean(controltriopgs$diff8)/(sd(controltriopgs$diff8)/N)

One_Control
Two_Control
Three_Control
Four_Control
Five_Control
Six_Control
Seven_Control
Eight_Control


rm(list = ls())


### Validation in the SPARK cohort


Mother = fread("~/SPARK/Phenotypes/SPARK_mother.txt", header= T)
setnames(Mother, c("subject_sp_id", "family_id"), c("IID", "FID"))

Father = read.table("~/SPARK/Phenotypes/SPARK_father.txt", header= T)
setnames(Father, c("subject_sp_id", "family_id"), c("IID", "FID"))

Cases = fread("~/SPARK/Phenotypes/Autistic_proband.txt", header= T)
setnames(Cases, c("subject_sp_id", "family_id"), c("IID", "FID"))

siblings = fread("~/SPARK/Phenotypes/Nonautistic_sibling.txt", header= T)
setnames(siblings, c("subject_sp_id", "family_id"), c("IID", "FID"))

#Next, read the prs scores
PRS = fread("~/ALSPAC/PRSice2results/SPARK_autosomes_retro_prospcildtraumasum.all.score", header = TRUE)

library(tidyr)

#PRS = PRS %>% separate(IID,c("FID1", "IID"))
#PRS = PRS[,-c("FID")]
#setnames(PRS, "FID1", "FID")

fatherpgs = merge(Father, PRS, by = "IID")
motherpgs = merge(Mother, PRS, by = "IID")
casespgs = merge(Cases, PRS, by = "IID")


parentpgs = merge(motherpgs, fatherpgs, by = "FID.x")
parentpgs$midparent8 = (parentpgs$`1.000000.x` + parentpgs$`1.000000.y`)/2
parentpgs$midparent7 = (parentpgs$`0.750000.x` + parentpgs$`0.750000.y`)/2
parentpgs$midparent6 = (parentpgs$`0.500000.x` + parentpgs$`0.500000.y`)/2
parentpgs$midparent5 = (parentpgs$`0.250000.x` + parentpgs$`0.250000.y`)/2
parentpgs$midparent4 = (parentpgs$`0.100000.x` + parentpgs$`0.100000.y`)/2
parentpgs$midparent3 = (parentpgs$`0.010000.x` + parentpgs$`0.010000.y`)/2
parentpgs$midparent2 = (parentpgs$`0.001000.x` + parentpgs$`0.001000.y`)/2
parentpgs$midparent1 = (parentpgs$`0.000100.x` + parentpgs$`0.000100.y`)/2

triopgs = merge(parentpgs, casespgs, by = "FID.x")

triopgs = triopgs[!duplicated(triopgs[,c("FID.x")]),]

Sd1 = sd(triopgs$midparent1)
Sd2 = sd(triopgs$midparent2)
Sd3 = sd(triopgs$midparent3)
Sd4 = sd(triopgs$midparent4)
Sd5 = sd(triopgs$midparent5)
Sd6 = sd(triopgs$midparent6)
Sd7 = sd(triopgs$midparent7)
Sd8 = sd(triopgs$midparent8)

## Step 3: Calculate pTDT deviation
triopgs$diff8 = (triopgs$`1.000000` - triopgs$midparent8)/Sd8
triopgs$diff7 = (triopgs$`0.750000` - triopgs$midparent7)/Sd7
triopgs$diff6 = (triopgs$`0.500000` - triopgs$midparent6)/Sd6
triopgs$diff5 = (triopgs$`0.250000` - triopgs$midparent5)/Sd5
triopgs$diff4 = (triopgs$`0.100000` - triopgs$midparent4)/Sd4
triopgs$diff3 = (triopgs$`0.010000` - triopgs$midparent3)/Sd3
triopgs$diff2 = (triopgs$`0.001000` - triopgs$midparent2)/Sd2
triopgs$diff1 = (triopgs$`0.000100` - triopgs$midparent1)/Sd1


## Step 4: Calculate the T score for pTDT deviation
N = sqrt(nrow(triopgs))

One = mean(triopgs$diff1)/(sd(triopgs$diff1)/N)
Two = mean(triopgs$diff2)/(sd(triopgs$diff2)/N)
Three = mean(triopgs$diff3)/(sd(triopgs$diff3)/N)
Four = mean(triopgs$diff4)/(sd(triopgs$diff4)/N)
Five = mean(triopgs$diff5)/(sd(triopgs$diff5)/N)
Six = mean(triopgs$diff6)/(sd(triopgs$diff6)/N)
Seven = mean(triopgs$diff7)/(sd(triopgs$diff7)/N)
Eight = mean(triopgs$diff8)/(sd(triopgs$diff8)/N)

One
Two
Three
Four
Five
Six
Seven
Eight


### In siblings

siblingpgs = merge(siblings, PRS, by = "IID")


controltriopgs = merge(parentpgs, siblingpgs, by = "FID.x")
controltriopgs = controltriopgs[!duplicated(controltriopgs[,c("FID.x")]),]

Sd1 = sd(controltriopgs$midparent1)
Sd2 = sd(controltriopgs$midparent2)
Sd3 = sd(controltriopgs$midparent3)
Sd4 = sd(controltriopgs$midparent4)
Sd5 = sd(controltriopgs$midparent5)
Sd6 = sd(controltriopgs$midparent6)
Sd7 = sd(controltriopgs$midparent7)
Sd8 = sd(controltriopgs$midparent8)

controltriopgs$diff8 = (controltriopgs$`1.000000` - controltriopgs$midparent8)/Sd8
controltriopgs$diff7 = (controltriopgs$`0.750000` - controltriopgs$midparent7)/Sd7
controltriopgs$diff6 = (controltriopgs$`0.500000` - controltriopgs$midparent6)/Sd6
controltriopgs$diff5 = (controltriopgs$`0.250000` - controltriopgs$midparent5)/Sd5
controltriopgs$diff4 = (controltriopgs$`0.100000` - controltriopgs$midparent4)/Sd4
controltriopgs$diff3 = (controltriopgs$`0.010000` - controltriopgs$midparent3)/Sd3
controltriopgs$diff2 = (controltriopgs$`0.001000` - controltriopgs$midparent2)/Sd2
controltriopgs$diff1 = (controltriopgs$`0.000100` - controltriopgs$midparent1)/Sd1


N = sqrt(nrow(controltriopgs))

One_Control = mean(controltriopgs$diff1)/(sd(controltriopgs$diff1)/N)
Two_Control = mean(controltriopgs$diff2)/(sd(controltriopgs$diff2)/N)
Three_Control = mean(controltriopgs$diff3)/(sd(controltriopgs$diff3)/N)
Four_Control = mean(controltriopgs$diff4)/(sd(controltriopgs$diff4)/N)
Five_Control = mean(controltriopgs$diff5)/(sd(controltriopgs$diff5)/N)
Six_Control = mean(controltriopgs$diff6)/(sd(controltriopgs$diff6)/N)
Seven_Control = mean(controltriopgs$diff7)/(sd(controltriopgs$diff7)/N)
Eight_Control = mean(controltriopgs$diff8)/(sd(controltriopgs$diff8)/N)

One_Control
Two_Control
Three_Control
Four_Control
Five_Control
Six_Control
Seven_Control
Eight_Control


rm(list = ls())
