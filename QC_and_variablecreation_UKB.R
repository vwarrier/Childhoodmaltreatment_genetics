#Files for ACE GWAS

# R program ukb22918.tab created 2018-09-17 by ukb2r.cpp Mar 14 2018 14:22:05
library(data.table)
library(plyr)

bd <- fread("~/UKB_v2/ukb24196.tab", header=TRUE, sep="\t")

## subset it to self-identified white europeans
file1 = subset(bd, f.21000.0.0 == "1"| f.21000.0.0 == "1001" | f.21000.0.0 == "1002" | f.21000.0.0 == "1003"| f.21000.0.0 == "1004") 

file1$checksex = ifelse(file1$f.31.0.0 == file1$f.22001.0.0, "correct", "incorrect")
file1 = subset(file1, checksex == "correct") # remove sex mismatches
pheno2 = file1[is.na(file1$f.22027.0.0),] #remove excessive heterozygosity

rm(bd)

PC1_mean = mean(na.omit(pheno2$f.22009.0.1)) 
PC1_sd = sd(na.omit(pheno2$f.22009.0.1))

PC2_mean = mean(na.omit(pheno2$f.22009.0.2))
PC2_sd = sd(na.omit(pheno2$f.22009.0.2))

pheno3 = subset(pheno2, f.22009.0.1 < (PC1_mean + 5*PC1_sd) & f.22009.0.1 > (PC1_mean - 5*PC1_sd)) #Remove outliers based on first PC
pheno2 = subset(pheno3, f.22009.0.2 < (PC2_mean + 5*PC2_sd) & f.22009.0.2 > (PC2_mean - 5*PC2_sd)) #Remove outliers based on second PC


pheno2$f.20487.0.0[pheno2$f.20487.0.0<0]<- NA #Felt hated by family member as a child
pheno2$f.20488.0.0[pheno2$f.20488.0.0<0]<- NA #Physically abused by family as a child
pheno2$f.20491.0.0[pheno2$f.20491.0.0<0]<- NA #Someone to take to doctor when needed as a child
pheno2$f.20490.0.0[pheno2$f.20490.0.0<0]<- NA #Sexually molested as a child
pheno2$f.20489.0.0[pheno2$f.20489.0.0<0]<- NA # Felt loved as a child 

pheno2$f.20489.0.0 <- 4 - pheno2$f.20489.0.0   # Felt loved as a child, inverse score to make it maltreatment
pheno2$f.20491.0.0 <- 4 - pheno2$f.20491.0.0  # Someone to take me to the doctor as a child, inverse score to make it maltreatment

pheno2$childtraumasum = pheno2$f.20489.0.0 + pheno2$f.20490.0.0 + pheno2$f.20488.0.0 + pheno2$f.20487.0.0 + pheno2$f.20491.0.0

pheno2 = pheno2[!is.na(pheno2$childtraumasum),]

pheno2$log_childtraumasum = log(1 + pheno2$childtraumasum)


####FILES FOR BOLT-LMM###

phenotypefile = pheno2[,c("log_childtraumasum", "f.eid", "f.eid", "childtraumasum", "f.20489.0.0", "f.20491.0.0", "f.20487.0.0", "f.20488.0.0", "f.20490.0.0")]
setnames(phenotypefile, 1, "FID")
setnames(phenotypefile, 2, "IID")


write.table(phenotypefile, file = "~/UKB_v2/ACE_Psychiatry/ACEphenobolt.txt", row.names = F, col.names = T, quote = F)

covar = pheno2[,c("f.eid", "f.eid","f.22000.0.0", "f.22001.0.0", "f.34.0.0","f.22009.0.1",
                  "f.22009.0.2", "f.22009.0.3", "f.22009.0.4", "f.22009.0.5",
                  "f.22009.0.6", "f.22009.0.7", "f.22009.0.8", "f.22009.0.9",
                  "f.22009.0.10", "f.22009.0.11", "f.22009.0.12", "f.22009.0.13",
                  "f.22009.0.14", "f.22009.0.15", "f.22009.0.16", "f.22009.0.17",
                  "f.22009.0.18", "f.22009.0.19", "f.22009.0.20", "f.22009.0.21",
                  "f.22009.0.22", "f.22009.0.23", "f.22009.0.24", "f.22009.0.25",
                  "f.22009.0.26", "f.22009.0.27", "f.22009.0.28", "f.22009.0.29",
                  "f.22009.0.30", "f.22009.0.31", "f.22009.0.32", "f.22009.0.33",
                  "f.22009.0.34", "f.22009.0.35", "f.22009.0.36", "f.22009.0.37",
                  "f.22009.0.38", "f.22009.0.39", "f.22009.0.40")]

setnames(covar, 1, "FID")
setnames(covar, 2, "IID")

write.table(covar, file = "~/UKB_v2/ACE_Psychiatry/covarbolt.txt", row.names = F, col.names = T, quote = F)


remove_file = bd[!(bd$f.eid %in% phenotypefile$FID),]

remove_file = remove_file[,c("f.eid", "f.eid")]
setnames(remove_file, 1, "FID")
setnames(remove_file, 2, "IID")

fam = fread("~/UKB_v2/Plink_files/ukbchr21.fam")
remove3 = remove_file[remove_file$FID %in% fam$V1,]

write.table(remove3, file = "~/UKB_v2/ACE_Psychiatry/removebolt.txt", row.names = F, col.names = T, quote = F)