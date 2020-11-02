###PGS in the UK Biobank#####################

./PRSice2_linux --base ./pgssumstats/log_childtraumasum_wosiblings_prsice.txt --target ~/UKB_v2/Plink_files/UKB2_autosomes --thread 20 --stat BETA --binary-target T --out ./PRSice2results/UKB2_autosomes_logchildtraumasum_wosiblings_v2 --bar-levels 0.00000005,0.000001,0.00001,0.0001,0.001,0.01,0.1,0.25,0.5,0.75,1 --no-regress --fastscore


###Polygenic scores for childhood maltreatment in the UKB

pheno = fread("~/UKB_v2/ACE_Psychiatry/ACEphenobolt.txt")
covariates = fread("~/UKB_v2/ACE_Psychiatry/covarbolt.txt")
pgs = fread("~/ALSPAC/PRSice2results/UKB2_autosomes_logchildtraumasum_wosiblings_v2.all.score")

merged = merge(pheno, covariates, by = "FID")
merged = merge(merged, pgs, by = "FID")

Unrelated = fread("~/GCTA/ace_siblingsonly_grm_unrelated_oldversion.grm.id")

merged_unrelated = merged[merged$FID %in% Unrelated$V1,]

summary(lm(scale(childtraumasum) ~ scale(`1.000000`) + f.34.0.0 + f.22000.0.0 + f.22001.0.0 + f.22009.0.1 +
             f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + 
             f.22009.0.9 + f.22009.0.10 + f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 +
             f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20 + f.22009.0.21 + f.22009.0.22 +
             f.22009.0.23 + f.22009.0.24 + f.22009.0.25 + f.22009.0.26 + f.22009.0.27 + f.22009.0.28 + f.22009.0.29 + 
             f.22009.0.30 + f.22009.0.31 + f.22009.0.32 + f.22009.0.33 + f.22009.0.34 + f.22009.0.35 + f.22009.0.36 + 
             f.22009.0.37 + f.22009.0.38 + f.22009.0.39 + f.22009.0.40, data = merged_unrelated))

summary(lm(scale(childtraumasum) ~  f.34.0.0 + f.22000.0.0 + f.22001.0.0 + f.22009.0.1 +
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

summary(lm(scale(childtraumasum) ~ scale(`0.750000`) + f.189.0.0 + f.34.0.0 + f.22000.0.0 + f.22001.0.0 + f.22009.0.1 +
             f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + 
             f.22009.0.9 + f.22009.0.10 + f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 +
             f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20 + f.22009.0.21 + f.22009.0.22 +
             f.22009.0.23 + f.22009.0.24 + f.22009.0.25 + f.22009.0.26 + f.22009.0.27 + f.22009.0.28 + f.22009.0.29 + 
             f.22009.0.30 + f.22009.0.31 + f.22009.0.32 + f.22009.0.33 + f.22009.0.34 + f.22009.0.35 + f.22009.0.36 + 
             f.22009.0.37 + f.22009.0.38 + f.22009.0.39 + f.22009.0.40, data = merged_unrelated))

summary(lm(scale(childtraumasum) ~  f.189.0.0 + f.34.0.0 + f.22000.0.0 + f.22001.0.0 + f.22009.0.1 +
             f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + 
             f.22009.0.9 + f.22009.0.10 + f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 +
             f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20 + f.22009.0.21 + f.22009.0.22 +
             f.22009.0.23 + f.22009.0.24 + f.22009.0.25 + f.22009.0.26 + f.22009.0.27 + f.22009.0.28 + f.22009.0.29 + 
             f.22009.0.30 + f.22009.0.31 + f.22009.0.32 + f.22009.0.33 + f.22009.0.34 + f.22009.0.35 + f.22009.0.36 + 
             f.22009.0.37 + f.22009.0.38 + f.22009.0.39 + f.22009.0.40, data = merged_unrelated))





################################
### Between-sibling analyses####
################################
#This script has been modified from Saskia Selzam's between siblings and Rosa Cheeman and Perline Demange's genetic nurture
###Perline's github: https://github.com/PerlineDemange/GeneticNurtureNonCog/blob/master/UKB/Siblings/siblings_UKB.R
###############################


sibling_data <- fread("~/UKB_v2/ACE_Psychiatry/sibling_problem/rels/sibs_famid3",h=F, data.table=F, verbose=T) #41498
setnames(sibling_data, old = c(1,2), new = c("IID","FID"))

#colnames(sibling_data)[1:2]<-c()


###Read polygenic scores, phenotype, and covariates
pheno = fread("~/UKB_v2/ACE_Psychiatry/ACEphenobolt.txt")
covariates = fread("~/UKB_v2/ACE_Psychiatry/covarbolt.txt")
pgs = fread("/mnt/b2/home4/arc/vw260/ALSPAC/PRSice2results/UKB2_autosomes_logchildtraumasum_wosiblings_v2.all.score")
pheno = pheno[,c("IID", "childtraumasum")]
pgs = pgs[,c("IID", "1.000000")]
setnames(pgs, 2, "pgs_trauma")

dat2<-merge(pheno, b ,by = "IID")

finalsib1 <- merge(dat2, pgs, by='IID')
covariates = covariates[-c("FID")]
finalsib <- merge(finalsib1, covariates, by = "IID")

finalsib$trauma_sc = scale(finalsib$childtraumasum)
finalsib$log_childtraumasum_sc = scale(finalsib$log_childtraumasum)
finalsib$trauma_pgs_scaled = scale(finalsib$pgs_trauma)



### ICC

ICCest <- function(model) {
  icc <- sqrt(diag(getVarCov(model)))^2 / (sqrt(diag(getVarCov(model)))^2 + model$sigma^2 )
  as.vector(icc)
}

# intercept model
m0 <- lme(log_childtraumasum_sc~1, random=~1|FID, method="ML", na.action=na.omit,data=finalsib)
ICCest(m0) #get ICC #0.2807447

m0 <- lme(trauma_sc~1, random=~1|FID, method="ML", na.action=na.omit,data=finalsib)
ICCest(m0) #get ICC #0.321328

m0 <- lme(trauma_sc~1, random=~1|FID, method="ML", na.action=na.omit,data=finalsib)
ICCest(m0) #get ICC #0.321328

m0 <- lme(trauma_pgs_scaled~1, random=~1|FID, method="ML", na.action=na.omit,data=finalsib)
ICCest(m0) #get ICC #0.4906724


# Create between family estimates of PGS: average across family member
######################################################################

library(dplyr)
mean_trauma_pgs_scaled<-group_by(finalsib,FID) %>% summarize(m=mean(trauma_pgs_scaled))
colnames(mean_trauma_pgs_scaled) <- c("FID", "GPS_B_trauma")

finalsib<-merge(finalsib,mean_trauma_pgs_scaled,by="FID")

# Create within-family variable
################################

finalsib$GPS_W_childtrauma <- finalsib$trauma_pgs_scaled  - finalsib$GPS_B_trauma
 

cor.test(finalsib$GPS_W_childtrauma, finalsib$GPS_B_trauma) #-2.394342e-18  p=1


# Mixed effects model
######################

final <- lme((trauma_sc) ~ GPS_B_trauma + GPS_W_childtrauma +
               f.34.0.0 + f.22000.0.0 + f.22001.0.0 + f.22009.0.1 +
               f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + 
               f.22009.0.9 + f.22009.0.10 + f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 +
               f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20 + f.22009.0.21 + f.22009.0.22 +
               f.22009.0.23 + f.22009.0.24 + f.22009.0.25 + f.22009.0.26 + f.22009.0.27 + f.22009.0.28 + f.22009.0.29 + 
               f.22009.0.30 + f.22009.0.31 + f.22009.0.32 + f.22009.0.33 + f.22009.0.34 + f.22009.0.35 + f.22009.0.36 + 
               f.22009.0.37 + f.22009.0.38 + f.22009.0.39 + f.22009.0.40, random=~1|FID,  method="ML", na.action=na.omit,data=finalsib)

summary(final)

names(summary(final))
summary(final)$tTable

direct_trauma <- summary(final)$tTable[3,1]#  Value    Std.Error    DF     t-value      p-value
                                          # 0.0543163650 0.0219276323  2745  2.47707387 1.330625e-02
total_trauma <- summary(final)$tTable[2,1] # 0.0940901282 0.0097579061 10064  9.64245071 6.579784e-22
indirect_trauma <- total_trauma - direct_trauma #0.0397737
ratio_trauma <- indirect_trauma/direct_trauma #0.7322


resultsib <- summary(final)$tTable

###
nboot <- 10000
bootcoef<-function(data,index){
  datx<-data[index,]
  mod<-lme((trauma_sc) ~ GPS_B_trauma + GPS_W_childtrauma +
             f.34.0.0 + f.22000.0.0 + f.22001.0.0 + f.22009.0.1 +
             f.22009.0.2 + f.22009.0.3 + f.22009.0.4 + f.22009.0.5 + f.22009.0.6 + f.22009.0.7 + f.22009.0.8 + 
             f.22009.0.9 + f.22009.0.10 + f.22009.0.11 + f.22009.0.12 + f.22009.0.13 + f.22009.0.14 + f.22009.0.15 +
             f.22009.0.16 + f.22009.0.17 + f.22009.0.18 + f.22009.0.19 + f.22009.0.20 + f.22009.0.21 + f.22009.0.22 +
             f.22009.0.23 + f.22009.0.24 + f.22009.0.25 + f.22009.0.26 + f.22009.0.27 + f.22009.0.28 + f.22009.0.29 + 
             f.22009.0.30 + f.22009.0.31 + f.22009.0.32 + f.22009.0.33 + f.22009.0.34 + f.22009.0.35 + f.22009.0.36 + 
             f.22009.0.37 + f.22009.0.38 + f.22009.0.39 + f.22009.0.40, random=~1|FID, method="ML", na.action=na.omit, data=datx,control=lmeControl(opt = "optim"))
  fixef(mod) #get fixed effects
}

# carry out bootstrap
boot.out<-boot(finalsib,bootcoef,nboot, parallel = "multicore", ncpus=20) 

saveRDS(boot.out, "bootstrapped_output_sib_UKB_EA_20200529.Rda")

#plot to check bootstrapping
options(bitmapType='cairo')
png("UKB.sib.bootstrap_20200529.png",
    width = 10,
    height = 6,
    units = 'in',
    res = 600)
plot(boot.out)
dev.off()

# look at  CI
boot.ci(boot.out, type = c("norm", "basic"))

# save results for each bootstrap in data frame
bootoutput <- as.data.frame(boot.out$t)
colnames(bootoutput) <- rownames(as.data.frame(boot.out$t0))
#write.table(bootoutput, "Data_scores_siblings_UKB_bootstrapped_2020529.csv", row.names=F, quote=F)


# Get values out of boot.out for all estimates + create indirect and ratio estimates
original <- as.data.frame(t(boot.out$t0)) # estimates of the original sample #best estimates of the effects

original$direct_trauma <- original$GPS_W_childtrauma
original$total_trauma <- original$GPS_B_trauma
original$indirect_trauma <- original$total_trauma - original$direct_trauma
original$ratio_trauma <- original$indirect_trauma / original$direct_trauma
original$ratio_tot_trauma <- original$indirect_trauma / original$total_trauma


bootoutput$direct_trauma <- bootoutput$GPS_W_childtrauma
bootoutput$total_trauma <- bootoutput$GPS_B_trauma
bootoutput$indirect_trauma <- bootoutput$total_trauma - bootoutput$direct_trauma
bootoutput$ratio_trauma <- bootoutput$indirect_trauma / bootoutput$direct_trauma
bootoutput$ratio_tot_trauma <- bootoutput$indirect_trauma / bootoutput$total_trauma

mean <- apply(bootoutput, 2, mean) # mean of the estimates of the bootstrap resamples
bias <- mean - original
se <- apply(bootoutput, 2, sd) #the standard deviation of the bootstrap estimates is the standard error of the sample estimates

error <- qnorm(0.975)*se
leftCI <- original - bias - error # normal Ci from boot.ci 
rightCI <- original - bias + error
# Other kind of CI given by boot.ci, not saved
# leftCI3 <- quantile(bootoutput$X, 0.025) # percentile CI from boot.ci 
# rightCI3 <- quantile(bootoutput$X, 0.975)
# leftCI4 <- 2*original$SCORE.Nontrans.Cog_sc - quantile(bootoutput$SCORE.Nontrans.Cog_sc, 0.975) #basic ci from boot.ci
# rightCI4 <- 2*original$SCORE.Nontrans.Cog_sc - quantile(bootoutput$SCORE.Nontrans.Cog_sc, 0.025)

statsoutput <- rbind(original, mean, bias, se, error, leftCI, rightCI)
statsoutput$Estimates <- c('original', 'mean', 'bias', 'se', 'error', 'leftCI', 'rightCI')
statsoutput
tot <- statsoutput[,c(ncol(statsoutput), 1:(ncol(statsoutput)-1))]
tot

write.table(tot, "summary_mean_CI_siblings_UKB_20200529.csv", row.names=T, quote=F)

### Compare estimates 
######################

difftrauma <- original$direct_trauma - original$indirect_trauma

SD_sampling_difftrauma <- sd(bootoutput$direct_trauma - bootoutput$indirect_trauma)


Z_difftrauma <- difftrauma/SD_sampling_difftrauma


P_difftrauma <- 2*pnorm(-abs(Z_difftrauma))
P_diffnoncog <- 2*pnorm(-abs(Z_diffnoncog))
P_diffratio <- 2*pnorm(-abs(Z_diffratio))

compare <- cbind(Z_diffcog, P_diffcog, Z_diffnoncog, P_diffnoncog, Z_diffratio, P_diffratio)

write.table(compare, "Ztests_sib_UKB_20200529.csv", row.names=T, quote=F)


