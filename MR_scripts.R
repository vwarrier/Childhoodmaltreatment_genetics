
###MR with twosample MR and MRPResso

library(TwoSampleMR)
library(data.table)
library(MRPRESSO)

######
#Step 1: First convert everything to BETAs and SEs and standardize format - see example script below


outcome = fread("~/ALSPAC/pgssumstats/daner_PGC_BIP32b_mds7a_0416a")

outcome$BETA = log(outcome$OR)

outcome$Z = qnorm(1 - outcome$P/2) 
outcome$SE = abs(outcome$BETA/outcome$Z)

write.table(outcome, file = "daner_BIP32b_fortwosampleMR.txt", row.names = F, col.names = T, quote = F)


#### Run two sample MR - forward direction


SNP = as.data.table(c("rs12031035", "rs61818983", "rs3851357", "rs611531", "rs4895718", "rs7763390", "rs1859100", "rs4305836", "rs3896224", "rs35560901", "rs77987546", "rs4702", "rs5928362", "rs6633421"))
setnames(SNP, "V1", "SNP")

childhood_mal = fread("~/UKB_v2/ACE_Psychiatry/GWAS_summarystats/Retro_prosp_meta_forPRSice.txt")
Exposure = merge(SNP, childhood_mal, by = "SNP")
#Exposure$Chr_pos = paste(Exposure$CHR, Exposure$BP, sep = ":")
#setnames(Exposure, old = c('SNP','BETA','SE', 'A1', 'A2', 'Chr_pos'), new = c('rsid','beta','se', 'effect_allele', 'other_allele', 'SNP'))
setnames(Exposure, old = c('BETA','SE', 'A1', 'A2'), new = c('beta','se', 'effect_allele', 'other_allele'))


Exposure <- format_data(Exposure, type="exposure")


outcome_dat <- read_outcome_data(
  snps = Exposure$SNP,
  filename = "daner_pgc_mdd_meta_w2_no23andMe_rmUKBB",
  sep = " ",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P"
)

#Change the outcome data as needed

dat = harmonise_data(Exposure, outcome_dat, action = 1)
res = mr(dat, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
res

res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)


dat$samplesize.exposure = 185000
dat$samplesize.outcome =  142646  #185000 - CAD #149821 - diabetes #28520 -CRP  #142646 -depression   #51000 - bipolar #45000 - autism   #55000 - ADHD #135236 - scz #change as needed
out <- directionality_test(dat)
table(out)



###plots

mr_scatter_plot(res, dat) #mrplot
mr_forest_plot(res_single) #forest plot
mr_leaveoneout_plot(res_loo) # leave one out plot



mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.8)

##### Reverse direction


#SNP_data = fread("daner_pgc_mdd_meta_w2_no23andMe_rmUKBB", header = T)



#SNP_data = subset(SNP_data, P < 0.00000005)
SNP_data = as.data.table(c("rs1950829", "rs76025409", "rs6689226", "rs1936365", "rs17499892", "rs4811079"))
setnames(SNP_data, "V1", "SNP")

Exposure = fread("daner_pgc_mdd_meta_w2_no23andMe_rmUKBB")
setnames(Exposure, old = c('SE', "A1", "A2"), new = c('se', "effect_allele", "other_allele"))
#setnames(Exposure, old = c('markername','se_dgc', 'p_dgc', 'noneffect_allele'), new = c('SNP','se', 'P', 'other_allele'))
Exposure = Exposure[Exposure$SNP %in% SNP_data$SNP,]


Exposure <- format_data(Exposure, type="exposure")



outcome_dat <- read_outcome_data(
  snps = Exposure$SNP,
  filename = "~/UKB_v2/ACE_Psychiatry/GWAS_summarystats/Retro_prosp_meta_forPRSice.txt",
  sep = " ",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P"
)

#Change the outcome data as needed

dat = harmonise_data(Exposure, outcome_dat, action = 1)
res = mr(dat, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
res

res_single <- mr_singlesnp(dat)
res_loo <- mr_leaveoneout(dat)

mr_heterogeneity(dat)
mr_pleiotropy_test(dat)


dat$samplesize.exposure = 142646 #28520 - CRP #185000 -cardiogram  #138884 - depression    #51000 - bipolar #45000 - autism   #55000 - ADHD #135236 - scz #185000 - cad #change as needed 
dat$samplesize.outcome =  185000 
out <- directionality_test(dat)
table(out)



###plots

mr_scatter_plot(res, dat) #mrplot
mr_forest_plot(res_single) #forest plot
mr_leaveoneout_plot(res_loo) # leave one out plot

###MRPRESSO
mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.outcome", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  SignifThreshold = 0.05)
