# SDPRX
# Step1 Obtain PRS beta
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/SDPRX/subsample_prune_PAGE_PRS.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune

pop1=EUR
i=1
approx=TRUE

for pop2 in AFR; do
for trait in Height BMI SBP DBP PLT; do
for rpt in {1..4}; do

# sample size
if [[ ${pop2} == "AFR" && ${trait} == "Height" ]]; then
sample_size1=252357; sample_size2=42146; rho=0.99
elif [[ ${pop2} == "AFR" && ${trait} == "BMI" ]]; then
sample_size1=233787; sample_size2=41802; rho=0.99
elif [[ ${pop2} == "AFR" && ${trait} == "SBP" ]]; then
sample_size1=728893; sample_size2=31276; rho=0.99
elif [[ ${pop2} == "AFR" && ${trait} == "DBP" ]]; then
sample_size1=746038; sample_size2=31276; rho=0.99
elif [[ ${pop2} == "AFR" && ${trait} == "PLT" ]]; then
sample_size1=539667; sample_size2=26605; rho=0.99
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/${trait}_UKB_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_chr${chr}_2.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/SDPRX/SDPRX.py --load_ld /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/SDPRX/EUR_${pop2} --valid /gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal.bim --ss1 data/summary_data/SDPRX/${trait}_${pop1}_inter_SDPRX_${pop2}_prune_snplist_${i}.txt --ss2 data/summary_data/subsample/SDPRX/${trait}_UKB_prune_snplist_${i}_${pop2}_train_SDPRX_approx${approx}_ratio3.00_repeat${rpt}.txt --N1 ${sample_size1} --N2 ${sample_size2} --mcmc_samples 2000 --burn 1000 --force_shared True --chr ${chr} --rho ${rho} --out result/summary_result/SDPRX/${trait}_UKB_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_chr${chr}" >> $job_file
fi

done

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/SDPRX/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/SDPRX/subsample_prune_PAGE_PRS.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_PAGE_PRS-$(date +%Y-%m-%d).sh


# Step2: Organize beta by chr pop for each param in each trait
library(data.table)

pop1="EUR"

for (pop2 in c("AFR")){
for (trait in c("Height","BMI","SBP","DBP","PLT")){
for (rpt in c(1:4)){

## GWAS_type = subsample_prune
GWAS_type = "subsample_prune"
i = 1
approx = "TRUE"

SDPRX_all <- data.table()
for (chr in 1:22){
    SDPRX_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_UKB_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_chr",chr,"_2.txt"))

    names(SDPRX_pop_chr) = c("SNP","A1","BETA")

    SDPRX_all = rbind(SDPRX_all,SDPRX_pop_chr)
    
}

write.table(SDPRX_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_UKB_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop2,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

SDPRX_all <- data.table()
for (chr in 1:22){
    SDPRX_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_UKB_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_chr",chr,"_1.txt"))

    names(SDPRX_pop_chr) = c("SNP","A1","BETA")

    SDPRX_all = rbind(SDPRX_all,SDPRX_pop_chr)
    
}

write.table(SDPRX_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/",trait,"_UKB_",pop1,"_",pop2,"_SDPRX_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop1,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}

}


# Step3: Clean the chr file
for trait in Height BMI SBP DBP PLT; do

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/SDPRX/
rm -rf ${trait}_*_chr*.txt

done