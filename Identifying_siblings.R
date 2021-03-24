# The sibling problem ##
  

##Identify siblings
rel = fread("~/UKB_v2/ukb20904_rel_s488302.dat")
siblings = subset(rel, IBS0 > 0.0012)
siblings = subset(siblings, Kinship  < 0.35 & Kinship > 0.176)

##22660 pairs
##22666 in the Bycroft et al paper

##Understand how many unique individuals are there
siblings_1= siblings[,1]
siblings_2= siblings[,2]
siblings_unique_ind = unique(rbind(siblings_1, siblings_2, use.names = F))
setnames(siblings_unique_ind, 1, "FID")

##41504 individuals who are unique

##Read_phenotype
pheno =fread("~/UKB_v2/ACE_Psychiatry/ACEphenobolt.txt")
merged = merge(siblings_unique_ind, pheno, by = "FID") 

###12,855 unique individuals

###Now, keep it in pairs
siblings_withpheno_pair = siblings[(siblings$ID1 %in% pheno$FID),]
siblings_withpheno_pair2 = siblings[(siblings$ID2 %in% pheno$FID),]

siblings_alpha= siblings_withpheno_pair[,1]
siblings_beta= siblings_withpheno_pair2[,2]

siblings_unique_ind_alphabet = unique(rbind(siblings_alpha, siblings_beta, use.names = F))

#12,855 unique individuals


##OK now the issue with the previous two steps are that you're selecting individuals where at least
##one of the siblings have a phenotype, and not where both of them have a phenotype. 
##To do this, lets do the steps below


##Keep it in pairs #2
siblings_withpheno_pair_v2 = siblings[(siblings$ID1 %in% pheno$FID),]
siblings_withpheno_pair_v3 = siblings_withpheno_pair_v2[(siblings_withpheno_pair_v2$ID2 %in% pheno$FID),]

write.table(siblings_withpheno_pair_v3, file = "~/UKB_v2/ACE_Psychiatry/sibling_problem/list_of_siblingpairs_withpheno_ACE.txt", row.names = F, col.names = T, quote = F )


siblings_alpha= siblings_withpheno_pair_v3[,1]
siblings_beta= siblings_withpheno_pair_v3[,2]

siblings_unique_ind_alphabet = unique(rbind(siblings_alpha, siblings_beta, use.names = F)) #5515
setnames(siblings_unique_ind_alphabet, 1, "FID")
siblings_unique_ind_alphabet$IID = siblings_unique_ind_alphabet$FID

write.table(siblings_unique_ind_alphabet, file = "~/UKB_v2/ACE_Psychiatry/sibling_problem/list_of_siblings_ACE.txt", row.names = F, col.names = T, quote = F )

## 2849 sibling pairs with phenotype representing 5515 individuals. Run the command below directly in bash, not in R. 
./gcta64  --grm ace_siblingsonly_grm  --keep ~/UKB_v2/ACE_Psychiatry/sibling_problem/list_of_siblings_ACE.txt  --grm-cutoff  0.05  --make-grm  --out ace_siblingsonly_grm_unrelated --thread-num 15
