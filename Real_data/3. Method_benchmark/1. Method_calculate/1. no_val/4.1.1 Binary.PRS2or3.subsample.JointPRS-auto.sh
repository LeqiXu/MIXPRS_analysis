# JointPRS
# Step1 Obtain PRS beta
# Step1.1 subsample_full [dsq job: 48464472; 48464579]
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_full_Binary_PRS3.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_full

pop1=EUR

for pop2 in EAS AFR; do
for trait in T2D BrC; do
for rpt in {1..4}; do

# sample size
if [[ ${pop2} == "EAS" && ${trait} == "T2D" ]]; then
sample_size1=88825; sample_size2=88870; sample_size3=38919
elif [[ ${pop2} == "EAS" && ${trait} == "BrC" ]]; then
sample_size1=245620; sample_size2=20354; sample_size3=7434
elif [[ ${pop2} == "AFR" && ${trait} == "T2D" ]]; then
sample_size1=88825; sample_size2=118493; sample_size3=29189
elif [[ ${pop2} == "AFR" && ${trait} == "BrC" ]]; then
sample_size1=245620; sample_size2=27138; sample_size3=5576
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

if [[ ${pop2} == "EAS" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_full_snplist.txt,data/summary_data/subsample/PRScsx/${trait}_full_snplist_EAS_train_PRScsx_approxFALSE_ratio3.00_repeat${rpt}.txt,data/summary_data/PRScsx/${trait}_AFR_inter_PRScsx_full_snplist.txt --rho_cons=1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3} --chrom=${chr} --pop=EUR,EAS,AFR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "AFR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_full_snplist.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_full_snplist.txt,data/summary_data/subsample/PRScsx/${trait}_full_snplist_AFR_train_PRScsx_approxFALSE_ratio3.00_repeat${rpt}.txt --rho_cons=1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3} --chrom=${chr} --pop=EUR,EAS,AFR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_repeat${rpt}" >> $job_file
fi

fi

done

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_full_Binary_PRS3.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_full_Binary_PRS3-$(date +%Y-%m-%d).sh


job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_full_Binary_PRS2.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_full

pop1=EUR
pop2=EAS

for trait in CAD LuC; do
for rpt in {1..4}; do

# sample size
if [[ ${pop2} == "EAS" && ${trait} == "CAD" ]]; then
sample_size1=61333; sample_size2=75819
elif [[ ${pop2} == "EAS" && ${trait} == "LuC" ]]; then
sample_size1=77095; sample_size2=11918
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"
if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_${pop1}_inter_PRScsx_full_snplist.txt,data/summary_data/subsample/PRScsx/${trait}_full_snplist_${pop2}_train_PRScsx_approxFALSE_ratio3.00_repeat${rpt}.txt --rho_cons=1,1 --n_gwas=${sample_size1},${sample_size2} --chrom=${chr} --pop=${pop1},${pop2} --out_dir=result/summary_result/JointPRS --out_name=${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_repeat${rpt}" >> $job_file
fi


done

done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_full_Binary_PRS2.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_full_Binary_PRS2-$(date +%Y-%m-%d).sh


# Step1.2 subsample_prune [dsq job: 48464582; 48464584]
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_prune_Binary_PRS3.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune

pop1=EUR

for pop2 in EAS AFR; do
for trait in T2D BrC; do
for rpt in {1..4}; do
for i in {1..4}; do
for approx in TRUE FALSE; do

# sample size
if [[ ${pop2} == "EAS" && ${trait} == "T2D" ]]; then
sample_size1=88825; sample_size2=88870; sample_size3=38919
elif [[ ${pop2} == "EAS" && ${trait} == "BrC" ]]; then
sample_size1=245620; sample_size2=20354; sample_size3=7434
elif [[ ${pop2} == "AFR" && ${trait} == "T2D" ]]; then
sample_size1=88825; sample_size2=118493; sample_size3=29189
elif [[ ${pop2} == "AFR" && ${trait} == "BrC" ]]; then
sample_size1=245620; sample_size2=27138; sample_size3=5576
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

if [[ ${pop2} == "EAS" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/subsample/PRScsx/${trait}_prune_snplist_${i}_EAS_train_PRScsx_approx${approx}_ratio3.00_repeat${rpt}.txt,data/summary_data/PRScsx/${trait}_AFR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt --rho_cons=1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3} --chrom=${chr} --pop=EUR,EAS,AFR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "AFR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/subsample/PRScsx/${trait}_prune_snplist_${i}_AFR_train_PRScsx_approx${approx}_ratio3.00_repeat${rpt}.txt --rho_cons=1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3} --chrom=${chr} --pop=EUR,EAS,AFR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}" >> $job_file
fi

fi


done

done
done
done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_prune_Binary_PRS3.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_Binary_PRS3-$(date +%Y-%m-%d).sh


job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_prune_Binary_PRS2.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune

pop1=EUR
pop2=EAS

for trait in CAD LuC; do
for rpt in {1..4}; do
for i in {1..4}; do
for approx in TRUE FALSE; do

# sample size
if [[ ${pop2} == "EAS" && ${trait} == "CAD" ]]; then
sample_size1=61333; sample_size2=75819
elif [[ ${pop2} == "EAS" && ${trait} == "LuC" ]]; then
sample_size1=77095; sample_size2=11918
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"
if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_${pop1}_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/subsample/PRScsx/${trait}_prune_snplist_${i}_${pop2}_train_PRScsx_approx${approx}_ratio3.00_repeat${rpt}.txt --rho_cons=1,1 --n_gwas=${sample_size1},${sample_size2} --chrom=${chr} --pop=${pop1},${pop2} --out_dir=result/summary_result/JointPRS --out_name=${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}" >> $job_file
fi


done

done
done
done
done


module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_prune_Binary_PRS2.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_Binary_PRS2-$(date +%Y-%m-%d).sh

# Step2: Organize beta by chr pop for each param in each trait
library(data.table)

for (pop2 in c("EAS","AFR")){
for (pop in c("EUR","EAS","AFR")){
for (trait in c("T2D","BrC")){
for (rpt in c(1:4)){

## GWAS_type = subsample_full
GWAS_type = "subsample_full"

JointPRS_all <- data.table()
for (chr in 1:22){
    
    if (pop2 == "EAS"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_subEAS_AFR_JointPRS_",GWAS_type,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    } 
    if (pop2 == "AFR"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_subAFR_JointPRS_",GWAS_type,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    }

    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("SNP","A1","BETA")

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

if (pop2 == "EAS"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_subEAS_AFR_JointPRS_",GWAS_type,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}
if (pop2 == "AFR"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_subAFR_JointPRS_",GWAS_type,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}


## GWAS_type = subsample_prune
GWAS_type = "subsample_prune"

for (i in c(1:4)){
for (approx in c("TRUE","FALSE")){

JointPRS_all <- data.table()
for (chr in 1:22){

    if (pop2 == "EAS"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_subEAS_AFR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    } 
    if (pop2 == "AFR"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_subAFR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    }

    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("SNP","A1","BETA")

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

if (pop2 == "EAS"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_subEAS_AFR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}
if (pop2 == "AFR"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_subAFR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}

}
}

}
}
}
}


library(data.table)

pop1="EUR"
pop2="EAS"
for (trait in c("CAD","LuC")){
for (rpt in c(1:4)){

## GWAS_type = subsample_full
GWAS_type = "subsample_full"

JointPRS_all <- data.table()
for (chr in 1:22){
    JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_",pop1,"_",pop2,"_JointPRS_",GWAS_type,"_repeat",rpt,"_",pop2,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))

    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("SNP","A1","BETA")

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_",pop1,"_",pop2,"_JointPRS_",GWAS_type,"_repeat",rpt,"_beta_",pop2,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

JointPRS_all <- data.table()
for (chr in 1:22){
    JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_",pop1,"_",pop2,"_JointPRS_",GWAS_type,"_repeat",rpt,"_",pop1,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))

    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("SNP","A1","BETA")

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_",pop1,"_",pop2,"_JointPRS_",GWAS_type,"_repeat",rpt,"_beta_",pop1,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

## GWAS_type = subsample_prune
GWAS_type = "subsample_prune"

for (i in c(1:4)){
for (approx in c("TRUE","FALSE")){

JointPRS_all <- data.table()
for (chr in 1:22){
    JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_",pop1,"_",pop2,"_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop2,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))

    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("SNP","A1","BETA")

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_",pop1,"_",pop2,"_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop2,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

JointPRS_all <- data.table()
for (chr in 1:22){
    JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_",pop1,"_",pop2,"_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop1,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))

    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("SNP","A1","BETA")

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_",pop1,"_",pop2,"_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop1,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}

}
}


# Step3: Clean the chr file
for trait in T2D BrC CAD LuC; do

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/
rm -rf ${trait}_*_chr*.txt

done