# JointPRS
# Step1 Obtain PRS beta
# Step1.1 subsample_full [dsq job: 18574581]
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_full_GLGC_PRS.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_full

pop1=EUR

for pop2 in EAS AFR SAS AMR; do
for trait in HDL LDL TC logTG; do
for rpt in {1..4}; do

# sample size
if [[ ${pop2} == "EAS" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=87303; sample_size3=90804; sample_size4=33953; sample_size5=47276
elif [[ ${pop2} == "EAS" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=59770; sample_size3=87759; sample_size4=33658; sample_size5=33989
elif [[ ${pop2} == "EAS" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=108434; sample_size3=92554; sample_size4=34135; sample_size5=48055
elif [[ ${pop2} == "EAS" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=60803; sample_size3=89467; sample_size4=34023; sample_size5=37273
elif [[ ${pop2} == "AFR" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=116404; sample_size3=68103; sample_size4=33953; sample_size5=47276
elif [[ ${pop2} == "AFR" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=79693; sample_size3=65819; sample_size4=33658; sample_size5=33989
elif [[ ${pop2} == "AFR" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=144579; sample_size3=69416; sample_size4=34135; sample_size5=48055
elif [[ ${pop2} == "AFR" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=81071; sample_size3=67100; sample_size4=34023; sample_size5=37273
elif [[ ${pop2} == "SAS" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=116404; sample_size3=90804; sample_size4=25465; sample_size5=47276
elif [[ ${pop2} == "SAS" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=79693; sample_size3=87759; sample_size4=25244; sample_size5=33989
elif [[ ${pop2} == "SAS" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=144579; sample_size3=92554; sample_size4=25601; sample_size5=48055
elif [[ ${pop2} == "SAS" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=81071; sample_size3=89467; sample_size4=25517; sample_size5=37273
elif [[ ${pop2} == "AMR" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=116404; sample_size3=90804; sample_size4=33953; sample_size5=35457
elif [[ ${pop2} == "AMR" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=79693; sample_size3=87759; sample_size4=33658; sample_size5=25492
elif [[ ${pop2} == "AMR" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=144579; sample_size3=92554; sample_size4=34135; sample_size5=36041
elif [[ ${pop2} == "AMR" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=81071; sample_size3=89467; sample_size4=34023; sample_size5=27955
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

if [[ ${pop2} == "EAS" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_SAS_AMR_JointPRS_${GWAS_type}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_full_snplist.txt,data/summary_data/subsample/PRScsx/${trait}_full_snplist_EAS_train_PRScsx_approxFALSE_ratio3.00_repeat${rpt}.txt,data/summary_data/PRScsx/${trait}_AFR_inter_PRScsx_full_snplist.txt,data/summary_data/PRScsx/${trait}_SAS_inter_PRScsx_full_snplist.txt,data/summary_data/PRScsx/${trait}_AMR_inter_PRScsx_full_snplist.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_subEAS_AFR_SAS_AMR_JointPRS_${GWAS_type}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "AFR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_SAS_AMR_JointPRS_${GWAS_type}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_full_snplist.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_full_snplist.txt,data/summary_data/subsample/PRScsx/${trait}_full_snplist_AFR_train_PRScsx_approxFALSE_ratio3.00_repeat${rpt}.txt,data/summary_data/PRScsx/${trait}_SAS_inter_PRScsx_full_snplist.txt,data/summary_data/PRScsx/${trait}_AMR_inter_PRScsx_full_snplist.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_EAS_subAFR_SAS_AMR_JointPRS_${GWAS_type}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "SAS" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_EAS_AFR_subSAS_AMR_JointPRS_${GWAS_type}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_full_snplist.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_full_snplist.txt,data/summary_data/PRScsx/${trait}_AFR_inter_PRScsx_full_snplist.txt,data/summary_data/subsample/PRScsx/${trait}_full_snplist_SAS_train_PRScsx_approxFALSE_ratio3.00_repeat${rpt}.txt,data/summary_data/PRScsx/${trait}_AMR_inter_PRScsx_full_snplist.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_EAS_AFR_subSAS_AMR_JointPRS_${GWAS_type}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "AMR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_EAS_AFR_SAS_subAMR_JointPRS_${GWAS_type}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_full_snplist.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_full_snplist.txt,data/summary_data/PRScsx/${trait}_AFR_inter_PRScsx_full_snplist.txt,data/summary_data/PRScsx/${trait}_SAS_inter_PRScsx_full_snplist.txt,data/summary_data/subsample/PRScsx/${trait}_full_snplist_AMR_train_PRScsx_approxFALSE_ratio3.00_repeat${rpt}.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_EAS_AFR_SAS_subAMR_JointPRS_${GWAS_type}_repeat${rpt}" >> $job_file
fi

fi

done

done
done
done


module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_full_GLGC_PRS.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_full_GLGC_PRS-$(date +%Y-%m-%d).sh


# Step1.2 subsample_prune [dsq job: 18575209; 18575210]
for approx in TRUE FALSE; do

job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_prune_GLGC_PRS_approx${approx}.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune

pop1=EUR

for pop2 in EAS AFR SAS AMR; do
for trait in HDL LDL TC logTG; do
for rpt in {1..4}; do
for i in {1..4}; do

# sample size
if [[ ${pop2} == "EAS" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=87303; sample_size3=90804; sample_size4=33953; sample_size5=47276
elif [[ ${pop2} == "EAS" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=59770; sample_size3=87759; sample_size4=33658; sample_size5=33989
elif [[ ${pop2} == "EAS" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=108434; sample_size3=92554; sample_size4=34135; sample_size5=48055
elif [[ ${pop2} == "EAS" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=60803; sample_size3=89467; sample_size4=34023; sample_size5=37273
elif [[ ${pop2} == "AFR" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=116404; sample_size3=68103; sample_size4=33953; sample_size5=47276
elif [[ ${pop2} == "AFR" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=79693; sample_size3=65819; sample_size4=33658; sample_size5=33989
elif [[ ${pop2} == "AFR" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=144579; sample_size3=69416; sample_size4=34135; sample_size5=48055
elif [[ ${pop2} == "AFR" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=81071; sample_size3=67100; sample_size4=34023; sample_size5=37273
elif [[ ${pop2} == "SAS" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=116404; sample_size3=90804; sample_size4=25465; sample_size5=47276
elif [[ ${pop2} == "SAS" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=79693; sample_size3=87759; sample_size4=25244; sample_size5=33989
elif [[ ${pop2} == "SAS" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=144579; sample_size3=92554; sample_size4=25601; sample_size5=48055
elif [[ ${pop2} == "SAS" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=81071; sample_size3=89467; sample_size4=25517; sample_size5=37273
elif [[ ${pop2} == "AMR" && ${trait} == "HDL" ]]; then
sample_size1=885546; sample_size2=116404; sample_size3=90804; sample_size4=33953; sample_size5=35457
elif [[ ${pop2} == "AMR" && ${trait} == "LDL" ]]; then
sample_size1=840012; sample_size2=79693; sample_size3=87759; sample_size4=33658; sample_size5=25492
elif [[ ${pop2} == "AMR" && ${trait} == "TC" ]]; then
sample_size1=929739; sample_size2=144579; sample_size3=92554; sample_size4=34135; sample_size5=36041
elif [[ ${pop2} == "AMR" && ${trait} == "logTG" ]]; then
sample_size1=860679; sample_size2=81071; sample_size3=89467; sample_size4=34023; sample_size5=27955
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do

if [[ ${pop2} == "EAS" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_SAS_AMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/subsample/PRScsx/${trait}_prune_snplist_${i}_EAS_train_PRScsx_approx${approx}_ratio3.00_repeat${rpt}.txt,data/summary_data/PRScsx/${trait}_AFR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_SAS_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_AMR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_subEAS_AFR_SAS_AMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "AFR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_SAS_AMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/subsample/PRScsx/${trait}_prune_snplist_${i}_AFR_train_PRScsx_approx${approx}_ratio3.00_repeat${rpt}.txt,data/summary_data/PRScsx/${trait}_SAS_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_AMR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_EAS_subAFR_SAS_AMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "SAS" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_EAS_AFR_subSAS_AMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_AFR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/subsample/PRScsx/${trait}_prune_snplist_${i}_SAS_train_PRScsx_approx${approx}_ratio3.00_repeat${rpt}.txt,data/summary_data/PRScsx/${trait}_AMR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_EAS_AFR_subSAS_AMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "AMR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/${trait}_EUR_EAS_AFR_SAS_subAMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_${pop2}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal --sst_file=data/summary_data/PRScsx/${trait}_EUR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_EAS_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_AFR_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/PRScsx/${trait}_SAS_inter_PRScsx_${pop2}_prune_snplist_${i}.txt,data/summary_data/subsample/PRScsx/${trait}_prune_snplist_${i}_AMR_train_PRScsx_approx${approx}_ratio3.00_repeat${rpt}.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/summary_result/JointPRS --out_name=${trait}_EUR_EAS_AFR_SAS_subAMR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}" >> $job_file
fi

fi


done

done
done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/JointPRS/subsample_prune_GLGC_PRS_approx${approx}.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_GLGC_PRS_approx${approx}-$(date +%Y-%m-%d).sh

done

# Step2: Organize beta by chr pop for each param in each trait
library(data.table)

for (pop2 in c("EAS","AFR","SAS","AMR")){
for (pop in c("EUR","EAS","AFR","SAS","AMR")){
for (trait in c("HDL","LDL","TC","logTG")){
for (rpt in c(1:4)){

## GWAS_type = subsample_full
GWAS_type = "subsample_full"

JointPRS_all <- data.table()
for (chr in 1:22){
    
    if (pop2 == "EAS"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_subEAS_AFR_SAS_AMR_JointPRS_",GWAS_type,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    } 
    if (pop2 == "AFR"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_subAFR_SAS_AMR_JointPRS_",GWAS_type,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    }
    if (pop2 == "SAS"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_AFR_subSAS_AMR_JointPRS_",GWAS_type,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    } 
    if (pop2 == "AMR"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_AFR_SAS_subAMR_JointPRS_",GWAS_type,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    }

    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("SNP","A1","BETA")

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

if (pop2 == "EAS"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_subEAS_AFR_SAS_AMR_JointPRS_",GWAS_type,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}
if (pop2 == "AFR"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_subAFR_SAS_AMR_JointPRS_",GWAS_type,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}
if (pop2 == "SAS"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_AFR_subSAS_AMR_JointPRS_",GWAS_type,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}
if (pop2 == "AMR"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_AFR_SAS_subAMR_JointPRS_",GWAS_type,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}

## GWAS_type = subsample_prune
GWAS_type = "subsample_prune"

for (i in c(1:4)){
for (approx in c("TRUE","FALSE")){

JointPRS_all <- data.table()
for (chr in 1:22){

    if (pop2 == "EAS"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_subEAS_AFR_SAS_AMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    } 
    if (pop2 == "AFR"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_subAFR_SAS_AMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    }
    if (pop2 == "SAS"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_AFR_subSAS_AMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    } 
    if (pop2 == "AMR"){
        JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_AFR_SAS_subAMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    }

    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("SNP","A1","BETA")

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

if (pop2 == "EAS"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_subEAS_AFR_SAS_AMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}
if (pop2 == "AFR"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_subAFR_SAS_AMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}
if (pop2 == "SAS"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_AFR_subSAS_AMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}
if (pop2 == "AMR"){
    write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/JointPRS/",trait,"_EUR_EAS_AFR_SAS_subAMR_JointPRS_",GWAS_type,"_",i,"_approx",approx,"_repeat",rpt,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}

}
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