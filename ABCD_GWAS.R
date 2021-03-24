##ABCD GWAS

#Step 1: Read the files
setwd("~/ABCD/ABCDgenotype/GWAS")
library(data.table)
library(RNOmni)

keep_file = fread("~/ABCD/ABCDgenotype/Genotype_preimputation/QC4_european.fam")
pheno = fread("~/ABCD/ABCDgenotype/GWAS/Childhood_trauma/ABCD_phenocovar.txt", fill = TRUE)
PCs = fread("~/ABCD/ABCDgenotype/Genotype_preimputation/ABCD_PCsforGWAS_European.txt")

# Step 2: Keep only individuals who are not ancestry outliers
PC1_mean = mean(na.omit(PCs$X1))
PC1_sd = sd(na.omit(PCs$X1))
PC2_mean = mean(na.omit(PCs$X2))
PC2_sd = sd(na.omit(PCs$X2))

PCs = subset(PCs, X1 < (PC1_mean + 5*PC1_sd) & X1 > (PC1_mean - 5*PC1_sd))
PCs = subset(PCs, X2 < (PC2_mean + 5*PC2_sd) & X2 > (PC2_mean - 5*PC2_sd))

keep_file = keep_file[,c("V1", "V2")]
setnames(keep_file, "V2", "IID")
keep_file = keep_file[(keep_file$IID %in% PCs$Sample_name),] #5591

pheno_merged = merge(keep_file, pheno, by = "IID")
setnames(pheno_merged, 2, "FID")

write.table(pheno_merged[,c("FID", "IID", "trauma_full")], file = "./Childhood_trauma/traum_pheno.txt", row.names = F, col.names = T, quote = F)
setnames(pheno_merged, 3, "Age")
setnames(pheno_merged, 4, "Sex")
covariates = pheno_merged[,c("FID", "IID", "Age", "Sex")]

setnames(PCs, 1, "IID")
covar_merged = merge(covariates, PCs, by = "IID")

qcovar = covar_merged[,c("FID", "IID", "Age", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", 
                         "X12", "X13", "X14", "X15", "X16", "X17", "X18", "X19", "X20")]

covar = covar_merged[,c("FID", "IID", "Sex")]


write.table(qcovar, file = "~/ABCD/ABCDgenotype/GWAS/Childhood_trauma/qcovar.txt", row.names = F, col.names = T, quote = F)
write.table(covar, file = "~/ABCD/ABCDgenotype/GWAS/Childhood_trauma/covar.txt", row.names = F, col.names = T, quote = F)



###Rank_transformed model

library(RNOmni)

phenotype = fread("~/ABCD/ABCDgenotype/GWAS/Childhood_trauma/traum_pheno.txt")
qcovar = fread("~/ABCD/ABCDgenotype/GWAS/Childhood_trauma/qcovar.txt")
covar = fread("~/ABCD/ABCDgenotype/GWAS/Childhood_trauma/covar.txt")

merged = merge(phenotype, covar, by = "IID")
merged = merged[,c("FID.x", "IID", "trauma_full", "Sex")]
merged = merge(merged, qcovar, by = "IID")

model = lm(trauma_full ~ Sex + Age, data = merged)
merged$residuals = model$residuals

merged$res_ranktrans = rankNorm(merged$residuals)

pheno2 = merged[,c("FID.x", "IID", "res_ranktrans")]
setnames(pheno2, 1, "FID")

write.table(pheno2, file = "~/ABCD/ABCDgenotype/GWAS/Childhood_trauma/trauma_pheno_rnktransformed.txt", row.names = F, col.names = T, quote = F)
write.table(qcovar[,-c("Age")], file = "~/ABCD/ABCDgenotype/GWAS/Childhood_trauma/qcovar_withoutage.txt", row.names = F, col.names = T, quote = F )


###Run model directly on bash - including the script below for convenience (don't run it in R!)
./gcta64 --fastGWA-mlm --mbfile Fast_GWAS_pathfile.txt --grm-sparse sp_grm_european --pheno ./Childhood_trauma/ABCD_pheno_forGCTA.txt --mpheno 21  --qcovar ./Childhood_trauma/qcovar.txt --covar ./Childhood_trauma/covar.txt --threads 20 --out ./Childhood_trauma/trauma_all__rnknorm --maf 0.001 

