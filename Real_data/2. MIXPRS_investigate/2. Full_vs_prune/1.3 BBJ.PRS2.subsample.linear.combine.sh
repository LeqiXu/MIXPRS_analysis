# Step1: Obtain linear combination weights
# Step1.1: subsample_full [dsq job: ]
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_full_BBJ_linear_weights_full_vs_prune.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_full

pop=EAS
pop1=EUR
pop2=EAS

for trait in WBC NEU LYM MON EOS RBC HCT MCH MCV HB ALT ALP GGT; do
for rpt in {1..4}; do

approx=FALSE
selection_criterion=NO
non_negative_weights=FALSE

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_${pop2}_${GWAS_type}_repeat${rpt}_${pop}_linear_weights_approx${approx}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_full_snplist_${pop}_tune_MIX_approxFALSE_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_repeat${rpt}_beta_${pop2}.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_repeat${rpt}_beta_${pop2}.txt --indep_approx=${approx} --selection_criterion=${selection_criterion} --non_negative_weights=${non_negative_weights} --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}_JointPRS_SDPRX_${pop2}_${GWAS_type}_repeat${rpt}" >> $job_file
fi

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_${pop1}_${pop2}_${GWAS_type}_repeat${rpt}_${pop}_linear_weights_approx${approx}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_full_snplist_${pop}_tune_MIX_approxFALSE_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_repeat${rpt}_beta_${pop1}.txt,result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_repeat${rpt}_beta_${pop2}.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_repeat${rpt}_beta_${pop1}.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_repeat${rpt}_beta_${pop2}.txt --indep_approx=${approx} --selection_criterion=${selection_criterion} --non_negative_weights=${non_negative_weights} --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}_JointPRS_SDPRX_${pop1}_${pop2}_${GWAS_type}_repeat${rpt}" >> $job_file
fi

done
done


module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_full_BBJ_linear_weights_full_vs_prune.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_full_BBJ_linear_weights_full_vs_prune-$(date +%Y-%m-%d).sh


# Step1.2: subsample_prune [dsq job: ]
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_BBJ_linear_weights_full_vs_prune.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune

pop=EAS
pop1=EUR
pop2=EAS

for trait in WBC NEU LYM MON EOS RBC HCT MCH MCV HB ALT ALP GGT; do
for rpt in {1..4}; do

i=1
approx=TRUE
selection_criterion=NO
non_negative_weights=FALSE

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_${pop2}_${GWAS_type}_${i}_repeat${rpt}_${pop}_linear_weights_approx${approx}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_prune_snplist_${i}_${pop}_tune_MIX_approx${approx}_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop2}.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop2}.txt --indep_approx=${approx} --selection_criterion=${selection_criterion} --non_negative_weights=${non_negative_weights} --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}_JointPRS_SDPRX_${pop2}_${GWAS_type}_${i}_repeat${rpt}" >> $job_file
fi

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_${pop1}_${pop2}_${GWAS_type}_${i}_repeat${rpt}_${pop}_linear_weights_approx${approx}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_prune_snplist_${i}_${pop}_tune_MIX_approx${approx}_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop1}.txt,result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop2}.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop1}.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop2}.txt --indep_approx=${approx} --selection_criterion=${selection_criterion} --non_negative_weights=${non_negative_weights} --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}_JointPRS_SDPRX_${pop1}_${pop2}_${GWAS_type}_${i}_repeat${rpt}" >> $job_file
fi

done
done


module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_BBJ_linear_weights_full_vs_prune.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_BBJ_linear_weights_full_vs_prune-$(date +%Y-%m-%d).sh


## Step2: Obtain final MIXPRS weight
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_BBJ_MIXPRS_linear_weights_full_vs_prune.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_full
pop=EAS
pop1=EUR
pop2=EAS

for trait in WBC NEU LYM MON EOS RBC HCT MCH MCV HB ALT ALP GGT; do

pop=${pop2}
approx=FALSE

sst_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt"

JointPRS_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_beta_EUR.txt"
JointPRS_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_beta_EAS.txt"
SDPRX_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_${pop}_beta_EUR.txt"
SDPRX_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_EAS_beta_EAS.txt"

weight_file1="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_repeat1_${pop}_linear_weights_approx${approx}.txt"
weight_file2="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_repeat2_${pop}_linear_weights_approx${approx}.txt"
weight_file3="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_repeat3_${pop}_linear_weights_approx${approx}.txt"
weight_file4="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_repeat4_${pop}_linear_weights_approx${approx}.txt"

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_linear_weights_approx${approx}_${pop}_MIXPRS.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_final_combine.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=${sst_file} --pop=${pop} --prs_beta_file=${JointPRS_EUR},${JointPRS_EAS},${SDPRX_EUR},${SDPRX_EAS} --weight_file=${weight_file1},${weight_file2},${weight_file3},${weight_file4} --indep_approx=${approx} --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}_${GWAS_type}_linear_weights_approx${approx}" >> $job_file
fi

done

GWAS_type=subsample_prune
pop=EAS
pop1=EUR
pop2=EAS

for trait in WBC NEU LYM MON EOS RBC HCT MCH MCV HB ALT ALP GGT; do

pop=${pop2}
approx=TRUE
i=1

sst_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt"

JointPRS_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_beta_EUR.txt"
JointPRS_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_beta_EAS.txt"
SDPRX_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_${pop}_beta_EUR.txt"
SDPRX_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_EAS_beta_EAS.txt"

weight_file1="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_${i}_repeat1_${pop}_linear_weights_approx${approx}.txt"
weight_file2="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_${i}_repeat2_${pop}_linear_weights_approx${approx}.txt"
weight_file3="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_${i}_repeat3_${pop}_linear_weights_approx${approx}.txt"
weight_file4="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_${i}_repeat4_${pop}_linear_weights_approx${approx}.txt"

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${GWAS_type}_${i}_linear_weights_approx${approx}_${pop}_MIXPRS.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_final_combine.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=${sst_file} --pop=${pop} --prs_beta_file=${JointPRS_EUR},${JointPRS_EAS},${SDPRX_EUR},${SDPRX_EAS} --weight_file=${weight_file1},${weight_file2},${weight_file3},${weight_file4} --indep_approx=${approx} --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}_${GWAS_type}_${i}_linear_weights_approx${approx}" >> $job_file
fi

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_BBJ_MIXPRS_linear_weights_full_vs_prune.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_BBJ_MIXPRS_linear_weights_full_vs_prune-$(date +%Y-%m-%d).sh

