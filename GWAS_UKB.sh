###Run bolt

./bolt \
--bed=ukb_cal_chr{1:22}_v2.bed \
--bim=ukb_snp_chr{1:22}_v2.bim \
--fam=ukb20904_cal_chr21_v2_s488292.fam \
--phenoFile=./ACE_Psychiatry/ACEphenobolt.txt \
--phenoCol=log_childtraumasum \
--covarFile=./ACE_Psychiatry/covarbolt.txt \
--covarCol=f.22000.0.0 \
--covarCol=f.22001.0.0 \
--qCovarCol=f.34.0.0 \
--qCovarCol=f.22009.0.{1:20} \
--remove=./ACE_Psychiatry/removebolt.txt \
--covarMaxLevels=200 \
--lmm \
--LDscoresFile=./Plink_files/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=./Plink_files/genetic_map_hg19_withX.txt.gz \
--lmmForceNonInf \
--numThreads=15 \
--modelSnps=modelsnps2.txt \
--bgenMinMAF=0.001 \
--bgenMinINFO=0.40 \
--statsFileBgenSnps=./ACE_Psychiatry/log_childtraumasum_try2_withX_bgenfile.txt \
--statsFile=./ACE_Psychiatry/log_childtraumasum_try2_withX_statfile.txt \
--bgenSampleFileList=UKB_bgenfiles.txt \
--remove=bolt.in_plink_but_not_imputed.FID_IID.93.txt \
--remove=bolt.in_plink_but_not_imputed.FID_IID.978.txt
