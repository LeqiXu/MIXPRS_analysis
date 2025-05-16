# SDPRX
# Step1.1 subsample_full [dsq job: 18575464]
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/SDPRX/subsample_full_GLGC_PRS.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_full

pop1=EUR

for pop2 in EAS AFR SAS AMR; do
for trait in HDL LDL TC logTG; do
for rpt in {1..4}; do

# sample size
if [[ ${pop2} == "EAS" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=87303; rho=0.99
elif  [[ ${pop2} == "AFR" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=68103; rho=0.93
elif  [[ ${pop2} == "SAS" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=25465; rho=0.99
elif  [[ ${pop2} == "AMR" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=35457; rho=0.95
elif [[ ${pop2} == "EAS" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=59770; rho=0.90
elif [[ ${pop2} == "AFR" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=65819; rho=0.70
elif [[ ${pop2} == "SAS" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=25244; rho=0.99
elif [[ ${pop2} == "AMR" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=25492; rho=0.92
elif [[ ${pop2} == "EAS" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=108434; rho=0.99
elif [[ ${pop2} == "AFR" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=69416; rho=0.74
elif [[ ${pop2} == "SAS" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=25601; rho=0.99
elif [[ ${pop2} == "AMR" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=36041; rho=0.92
elif [[ ${pop2} == "EAS" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=60803; rho=0.99
elif [[ ${pop2} == "AFR" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=67100; rho=0.93
elif [[ ${pop2} == "SAS" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=25517; rho=0.91
elif [[ ${pop2} == "AMR" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=27955; rho=0.96
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_repeat${rpt}_chr${chr}_2.txt"
if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/SDPRX/SDPRX.py --load_ld /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/SDPRX/EUR_${pop2} --valid /gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal.bim --ss1 data/summary_data/SDPRX/${trait}_${pop1}_inter_SDPRX_full_snplist.txt --ss2 data/summary_data/subsample/SDPRX/${trait}_full_snplist_${pop2}_train_SDPRX_approxFALSE_ratio3.00_repeat${rpt}.txt --N1 ${sample_size1} --N2 ${sample_size2} --mcmc_samples 2000 --burn 1000 --force_shared True --chr ${chr} --rho ${rho} --out result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_repeat${rpt}_chr${chr}" >> $job_file
fi

done

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/SDPRX/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/SDPRX/subsample_full_GLGC_PRS.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_full_GLGC_PRS-$(date +%Y-%m-%d).sh



# Step1.2 subsample_prune [dsq job: 18575532; 18575533]
for approx in TRUE FALSE; do

job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/SDPRX/subsample_prune_GLGC_PRS_approx${approx}.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune

pop1=EUR

for pop2 in EAS AFR SAS AMR; do
for trait in HDL LDL TC logTG; do
for rpt in {1..4}; do
for i in {1..4}; do

# sample size
if [[ ${pop2} == "EAS" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=87303; rho=0.99
elif  [[ ${pop2} == "AFR" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=68103; rho=0.93
elif  [[ ${pop2} == "SAS" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=25465; rho=0.99
elif  [[ ${pop2} == "AMR" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=35457; rho=0.95
elif [[ ${pop2} == "EAS" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=59770; rho=0.90
elif [[ ${pop2} == "AFR" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=65819; rho=0.70
elif [[ ${pop2} == "SAS" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=25244; rho=0.99
elif [[ ${pop2} == "AMR" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=25492; rho=0.92
elif [[ ${pop2} == "EAS" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=108434; rho=0.99
elif [[ ${pop2} == "AFR" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=69416; rho=0.74
elif [[ ${pop2} == "SAS" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=25601; rho=0.99
elif [[ ${pop2} == "AMR" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=36041; rho=0.92
elif [[ ${pop2} == "EAS" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=60803; rho=0.99
elif [[ ${pop2} == "AFR" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=67100; rho=0.93
elif [[ ${pop2} == "SAS" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=25517; rho=0.91
elif [[ ${pop2} == "AMR" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=27955; rho=0.96
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_chr${chr}_2.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/SDPRX/SDPRX.py --load_ld /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/SDPRX/EUR_${pop2} --valid /gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal.bim --ss1 data/summary_data/SDPRX/${trait}_${pop1}_inter_SDPRX_${pop2}_prune_snplist_${i}.txt --ss2 data/summary_data/subsample/SDPRX/${trait}_prune_snplist_${i}_${pop2}_train_SDPRX_approx${approx}_ratio3.00_repeat${rpt}.txt --N1 ${sample_size1} --N2 ${sample_size2} --mcmc_samples 2000 --burn 1000 --force_shared True --chr ${chr} --rho ${rho} --out result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_chr${chr}" >> $job_file
fi

done

done
done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/SDPRX/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/SDPRX/subsample_prune_GLGC_PRS_approx${approx}.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_GLGC_PRS_approx${approx}-$(date +%Y-%m-%d).sh

done

# Step2: Organize beta by chr pop for each param in each trait
library(data.table)

pop1="EUR"

for (pop2 in c("EAS","AFR","SAS","AMR")){
for (trait in c("HDL","LDL","TC","logTG")){
for (rpt in c(1:4)){

## GWAS_type = subsample_full
GWAS_type = "subsample_full"

SDPRX_all <- data.table()
for (chr in 1:22){
    SDPRX_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_repeat",rpt,"_chr",chr,"_2.txt"))

    names(SDPRX_pop_chr) = c("SNP","A1","BETA")

    SDPRX_all = rbind(SDPRX_all,SDPRX_pop_chr)
    
}

write.table(SDPRX_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_repeat",rpt,"_beta_",pop2,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

SDPRX_all <- data.table()
for (chr in 1:22){
    SDPRX_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_repeat",rpt,"_chr",chr,"_1.txt"))

    names(SDPRX_pop_chr) = c("SNP","A1","BETA")

    SDPRX_all = rbind(SDPRX_all,SDPRX_pop_chr)
    
}

write.table(SDPRX_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_repeat",rpt,"_beta_",pop1,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

## GWAS_type = subsample_prune
GWAS_type = "subsample_prune"

for (i in c(1:4)){
for (approx in c("TRUE","FALSE")){

SDPRX_all <- data.table()
for (chr in 1:22){
    SDPRX_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_chr",chr,"_2.txt"))

    names(SDPRX_pop_chr) = c("SNP","A1","BETA")

    SDPRX_all = rbind(SDPRX_all,SDPRX_pop_chr)
    
}

write.table(SDPRX_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop2,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

SDPRX_all <- data.table()
for (chr in 1:22){
    SDPRX_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_chr",chr,"_1.txt"))

    names(SDPRX_pop_chr) = c("SNP","A1","BETA")

    SDPRX_all = rbind(SDPRX_all,SDPRX_pop_chr)
    
}

write.table(SDPRX_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop1,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}

}
}
}

# Step3: Clean the chr file
for trait in HDL LDL TC logTG; do

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/
rm -rf ${trait}_*_chr*.txt

done