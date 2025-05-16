# JointPRS
# Step1 Obtain PRS beta
approx=TRUE

job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_prune_GLGC_PRS_approx${approx}.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune

pop1=EUR
i=1

for pop2 in AFR AMR; do
for trait in HDL LDL TC logTG; do
for rpt in {1..4}; do

# sample size
if [[ ${pop2} == "AFR" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=116404; sample_size3=72877; sample_size4=40172; sample_size5=47276
elif [[ ${pop2} == "AFR" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=79693; sample_size3=70967; sample_size4=40472; sample_size5=33989
elif [[ ${pop2} == "AFR" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=144579; sample_size3=74573; sample_size4=40962; sample_size5=48055
elif [[ ${pop2} == "AFR" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=81071; sample_size3=72256; sample_size4=40845; sample_size5=37273
elif [[ ${pop2} == "AMR" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=116404; sample_size3=97169; sample_size4=40172; sample_size5=35457
elif [[ ${pop2} == "AMR" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=79693; sample_size3=94622; sample_size4=40472; sample_size5=25492
elif [[ ${pop2} == "AMR" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=144579; sample_size3=99430; sample_size4=40962; sample_size5=36041
elif [[ ${pop2} == "AMR" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=81071; sample_size3=96341; sample_size4=40845; sample_size5=27955
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

if [[ ${pop2} == "AFR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_UKB_EUR_EAS_subAFR_SAS_AMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/subsample/PRScsx/${trait}_UKB_prune_snplist_${i}_AFR_train_PRScsx_approx${approx}_ratio3.00_repeat${rpt}.txt,data/summary_data/PRScsx/${trait}_SAS_inter_UKB_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_AMR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/summary_result/JointPRS --out_name=${trait}_UKB_EUR_EAS_subAFR_SAS_AMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "AMR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_UKB_EUR_EAS_AFR_SAS_subAMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_AFR_inter_UKB_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_SAS_inter_UKB_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/subsample/PRScsx/${trait}_prune_snplist_${i}_AMR_train_PRScsx_approx${approx}_ratio3.00_repeat${rpt}.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/summary_result/JointPRS --out_name=${trait}_UKB_EUR_EAS_AFR_SAS_subAMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}" >> $job_file
fi

fi


done

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_prune_GLGC_PRS_approx${approx}.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_GLGC_PRS_approx${approx}-$(date +%Y-%m-%d).sh


# Step2: Organize beta by chr pop for each param in each trait
library(data.table)

for (pop2 in c("AFR","AMR")){
for (pop in c("EUR","EAS","AFR","SAS","AMR")){
for (trait in c("HDL","LDL","TC","logTG")){
for (rpt in c(1:4)){

## GWAS_type = subsample_prune
GWAS_type = "subsample_prune"
i = 1
approx = "TRUE"

JointPRS_all <- data.table()
for (chr in 1:22){

    if (pop2 == "AFR"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_UKB_EUR_EAS_subAFR_SAS_AMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    }
    if (pop2 == "AMR"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_UKB_EUR_EAS_AFR_SAS_subAMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    }

    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("SNP","A1","BETA")

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

if (pop2 == "AFR"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_UKB_EUR_EAS_subAFR_SAS_AMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}
if (pop2 == "AMR"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_UKB_EUR_EAS_AFR_SAS_subAMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}

}
}

}
}


# Step3: Clean the chr file
for trait in HDL LDL TC logTG; do

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/
rm -rf ${trait}_*_chr*.txt

done