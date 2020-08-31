###GCTA_GREML commands
##Subtypes and childhood_maltreatment

./gcta64  --reml-bivar 17 7 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/logchildraumasum_sexabuse_binary --thread-num 15
./gcta64  --reml-bivar 17 8 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/logchildraumasum_emotionalabuse_binary --thread-num 15
./gcta64  --reml-bivar 17 9 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/logchildraumasum_physicalabuse_binary --thread-num 15
./gcta64  --reml-bivar 17 10 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/logchildraumasum_physicalneglect_binary --thread-num 15
./gcta64  --reml-bivar 17 11 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/logchildraumasum_emotionalneglect_binary --thread-num 15
./gcta64  --reml-bivar 7 8 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/emotionalabusebinary_sexabuse_binary --thread-num 15
./gcta64  --reml-bivar 7 9 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/physicalabusebinary_sexabuse_binary --thread-num 15
./gcta64  --reml-bivar 7 10 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/physicalneglectbinary_sexabuse_binary --thread-num 15
./gcta64  --reml-bivar 7 11 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/emotionalneglectbinary_sexabuse_binary --thread-num 15
./gcta64  --reml-bivar 8 9 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/emotionalabusebinary_physicalabusebinary --thread-num 15
./gcta64  --reml-bivar 8 10 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/emotionalabusebinary_physicalneglectbinary --thread-num 15
./gcta64  --reml-bivar 8 11 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/emotionalabusebinary_emotionalneglectbinary --thread-num 15
./gcta64  --reml-bivar 9 10 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/physicalneglectbinary_physicalabusebinary --thread-num 15
./gcta64  --reml-bivar 9 11 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/emotionalneglectbinary_physicalabusebinary --thread-num 15
./gcta64  --reml-bivar 10 11 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/emotionalneglectbinary_physicalneglectbinary --thread-num 15



##Operationalizations of childhood_maltreatment
./gcta64  --reml-bivar 17 18 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/logchildraumasum_childtraumasumbinary --thread-num 15
./gcta64  --reml-bivar 17 39 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/logchildtraumasum_childtraumasumsevere --thread-num 15
./gcta64  --reml-bivar 17 40 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/logchildtraumasum_childtraumasumsevereabuseonly --thread-num 15
./gcta64  --reml-bivar 18 39 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/childtraumasumbinary_childtraumasumsevere --thread-num 15
./gcta64  --reml-bivar 18 40 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/childtraumasumbinary_childtraumasumsevereabuseonly --thread-num 15
./gcta64  --reml-bivar 39 40 --grm trauma_grm_unrelated  --pheno trauma_phenotypes.txt --covar traum_covar.txt --qcovar traum_qcovar_withoutses.txt  --out ./trauma_results/childtraumasumsevere_childtraumasumsevereabuseonly --thread-num 15


