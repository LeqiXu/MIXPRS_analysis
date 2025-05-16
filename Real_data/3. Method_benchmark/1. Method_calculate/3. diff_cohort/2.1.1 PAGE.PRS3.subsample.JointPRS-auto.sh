# JointPRS
# Step1 Obtain PRS beta
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_prune_PAGE_PRS.txt"
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
sample_size1=252357; sample_size2=159095; sample_size3=42146
elif [[ ${pop2} == "AFR" && ${trait} == "BMI" ]]; then
sample_size1=233787; sample_size2=158284; sample_size3=41802
elif [[ ${pop2} == "AFR" && ${trait} == "SBP" ]]; then
sample_size1=728893; sample_size2=179000; sample_size3=31276
elif [[ ${pop2} == "AFR" && ${trait} == "DBP" ]]; then
sample_size1=746038; sample_size2=179000; sample_size3=31276
elif [[ ${pop2} == "AFR" && ${trait} == "PLT" ]]; then
sample_size1=539667; sample_size2=179000; sample_size3=26605
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

if [[ ${pop2} == "AFR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_UKB_EUR_EAS_subAFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/subsample/PRScsx/${trait}_UKB_prune_snplist_${i}_AFR_train_PRScsx_approx${approx}_ratio3.00_repeat${rpt}.txt --rho_cons=1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3} --chrom=${chr} --pop=EUR,EAS,AFR --out_dir=result/summary_result/JointPRS --out_name=${trait}_UKB_EUR_EAS_subAFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}" >> $job_file
fi

fi


done

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_prune_PAGE_PRS.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_PAGE_PRS-$(date +%Y-%m-%d).sh


# Step2: Organize beta by chr pop for each param in each trait
library(data.table)

for (pop2 in c("AFR")){
for (pop in c("EUR","EAS","AFR")){
for (trait in c("Height","BMI","SBP","DBP","PLT")){
for (rpt in c(1:4)){

## GWAS_type = subsample_prune
GWAS_type = "subsample_prune"
i = 1
approx = "TRUE"

JointPRS_all <- data.table()
for (chr in 1:22){

    if (pop2 == "AFR"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_UKB_EUR_EAS_subAFR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    }

    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("SNP","A1","BETA")

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

if (pop2 == "AFR"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_UKB_EUR_EAS_subAFR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}

}
}

}
}


# Step3: Clean the chr file
for trait in Height BMI SBP DBP PLT; do

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/
rm -rf ${trait}_*_chr*.txt

done