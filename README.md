# Childhoodmaltreatment_genetics
This repository provides the scripts used in the genetic study of childhood maltreatment. Please adapt the scripts for your own use. 

1. QC_and_variablecreation_UKB.R - R script for running some QC (need to do other QC at an earlier level), and creating variables.

2. GWAS_UKB.sh - Script for running BOLT-LMM (see: https://www.nature.com/articles/ng.3190)

3. GCTA_GREML_commands.sh - Scripts for running GCTA-GREML heritability and bivariate genetic correlation. See GCTA-GREML for how to create GRM and doing QC.

4. Identifying_siblings.R - Scripts for identifying siblings in the UKB. I'd recommend this excellent script for further, detailed info: https://github.com/PerlineDemange/GeneticNurtureNonCog/blob/master/UKB/Siblings/siblings_UKB.R

5. Polygenic_scoring_scripts.R - Scripts for both between individuals and between-sibling based PGS. Between-siblings based on work done by Saskia Selzam, Perline Demange, and Rosa Cheesman. Please see the github link above. 

6. ABCD_GWAS.R - An R script for data QC etc in ABCD. At the end of the script is a line for running GWAS using FastGWA in ABCD, to be run on the terminal.

7. ACE_pTDT.R - pTDT script for SSC and SFARI

8. MR_scripts.R - TwosampleMR package based MR analyses

9. Code_for_figs.R - R code to create some of the figures
