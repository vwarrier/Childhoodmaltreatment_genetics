###PGS

./PRSice2_linux --base ./pgssumstats/log_childtraumasum_wosiblings_prsice.txt --target ~/UKB_v2/Plink_files/UKB2_autosomes --thread 20 --stat BETA --binary-target T --out ./PRSice2results/UKB2_autosomes_logchildtraumasum_wosiblings_v2 --bar-levels 0.00000005,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.25,0.5,0.75,1 --no-regress --fastscore



###Polygenic scores for childhood maltreatment

pheno = fread("~/UKB_v2/ACE_Psychiatry/ACEphenobolt.txt")
covariates = fread("~/UKB_v2/ACE_Psychiatry/covarbolt.txt")
pgs = fread("/mnt/b2/home4/arc/vw260/ALSPAC/PRSice2results/UKB2_autosomes_logchildtraumasum_wosiblings_v2.all.score")

pheno$log_childtraumasum = log(1 + pheno$childtraumasum)

merged = merge(pheno, covariates, by = "FID")
merged = merge(merged, pgs, by = "FID")

Unrelated = fread("~/GCTA/ace_siblingsonly_grm_unrelated_oldversion.grm.id")

merged_unrelated = merged[merged$FID %in% Unrelated$V1,]


summary(lm(log_childtraumasum ~ scale(`1.000000`) + f.34.0.0 + f.22000.0.0 + f.22001.0.0 + f.22009.0.1 +
             f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + 
             f.22009.0.9 + f.22009.0.10 + f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 +
             f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20 + f.22009.0.21 + f.22009.0.22 +
             f.22009.0.23 + f.22009.0.24 + f.22009.0.25 + f.22009.0.26 + f.22009.0.27 + f.22009.0.28 + f.22009.0.29 + 
             f.22009.0.30 + f.22009.0.31 + f.22009.0.32 + f.22009.0.33 + f.22009.0.34 + f.22009.0.35 + f.22009.0.36 + 
             f.22009.0.37 + f.22009.0.38 + f.22009.0.39 + f.22009.0.40, data = merged_unrelated))

summary(lm(log_childtraumasum ~  f.34.0.0 + f.22000.0.0 + f.22001.0.0 + f.22009.0.1 +
             f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + 
             f.22009.0.9 + f.22009.0.10 + f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 +
             f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20 + f.22009.0.21 + f.22009.0.22 +
             f.22009.0.23 + f.22009.0.24 + f.22009.0.25 + f.22009.0.26 + f.22009.0.27 + f.22009.0.28 + f.22009.0.29 + 
             f.22009.0.30 + f.22009.0.31 + f.22009.0.32 + f.22009.0.33 + f.22009.0.34 + f.22009.0.35 + f.22009.0.36 + 
             f.22009.0.37 + f.22009.0.38 + f.22009.0.39 + f.22009.0.40, data = merged_unrelated))


##Now with SES
covar_ses= fread("~/GCTA/traum_qcovar.txt")
covar_ses = covar_ses[,c("f.eid", "f.189.0.0")]
setnames(covar_ses, "f.eid", "FID")
merged_unrelated = merge(merged_unrelated, covar_ses, by = "FID")

summary(lm(log_childtraumasum ~ scale(`1.000000`) + f.189.0.0 + f.34.0.0 + f.22000.0.0 + f.22001.0.0 + f.22009.0.1 +
             f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + 
             f.22009.0.9 + f.22009.0.10 + f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 +
             f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20 + f.22009.0.21 + f.22009.0.22 +
             f.22009.0.23 + f.22009.0.24 + f.22009.0.25 + f.22009.0.26 + f.22009.0.27 + f.22009.0.28 + f.22009.0.29 + 
             f.22009.0.30 + f.22009.0.31 + f.22009.0.32 + f.22009.0.33 + f.22009.0.34 + f.22009.0.35 + f.22009.0.36 + 
             f.22009.0.37 + f.22009.0.38 + f.22009.0.39 + f.22009.0.40, data = merged_unrelated))

summary(lm(log_childtraumasum ~  f.189.0.0 + f.34.0.0 + f.22000.0.0 + f.22001.0.0 + f.22009.0.1 +
             f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + 
             f.22009.0.9 + f.22009.0.10 + f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 +
             f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20 + f.22009.0.21 + f.22009.0.22 +
             f.22009.0.23 + f.22009.0.24 + f.22009.0.25 + f.22009.0.26 + f.22009.0.27 + f.22009.0.28 + f.22009.0.29 + 
             f.22009.0.30 + f.22009.0.31 + f.22009.0.32 + f.22009.0.33 + f.22009.0.34 + f.22009.0.35 + f.22009.0.36 + 
             f.22009.0.37 + f.22009.0.38 + f.22009.0.39 + f.22009.0.40, data = merged_unrelated))



