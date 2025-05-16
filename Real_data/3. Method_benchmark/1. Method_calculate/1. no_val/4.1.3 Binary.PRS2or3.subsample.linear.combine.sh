## Step1: Linear combination with different snplists
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_Binary_3_final_NNLS_weights.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune
pop1=EUR

for pop2 in EAS AFR; do
for trait in T2D BrC; do
for rpt in {1..4}; do

pop=${pop2}
i=1
approx=TRUE
weight_name="non_negative_linear_weights_approx${approx}"

if [[ ${pop2} == "EAS" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_AFR_${GWAS_type}_${i}_repeat${rpt}_${pop}_${weight_name}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_prune_snplist_${i}_${pop}_tune_MIX_approx${approx}_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EUR.txt,result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EAS.txt,result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_AFR.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EUR.txt,result/summary_result/SDPRX/${trait}_${pop1}_EAS_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EAS.txt,result/summary_result/SDPRX/${trait}_${pop1}_AFR_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_AFR.txt --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}_JointPRS_SDPRX_EUR_EAS_AFR_${GWAS_type}_${i}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "AFR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_AFR_${GWAS_type}_${i}_repeat${rpt}_${pop}_${weight_name}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_prune_snplist_${i}_${pop}_tune_MIX_approx${approx}_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EUR.txt,result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EAS.txt,result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_AFR.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EUR.txt,result/summary_result/SDPRX/${trait}_${pop1}_EAS_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EAS.txt,result/summary_result/SDPRX/${trait}_${pop1}_AFR_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_AFR.txt --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}_JointPRS_SDPRX_EUR_EAS_AFR_${GWAS_type}_${i}_repeat${rpt}" >> $job_file
fi

fi

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_Binary_3_final_NNLS_weights.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_Binary_3_final_NNLS_weights-$(date +%Y-%m-%d).sh

job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_Binary_2_final_NNLS_weights.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune
pop=EAS
pop1=EUR
pop2=EAS

for trait in CAD LuC; do
for rpt in {1..4}; do

i=1
approx=TRUE
weight_name="non_negative_linear_weights_approx${approx}"

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_${pop1}_${pop2}_${GWAS_type}_${i}_repeat${rpt}_${pop}_${weight_name}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_prune_snplist_${i}_${pop}_tune_MIX_approx${approx}_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop1}.txt,result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop2}.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop1}.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop2}.txt --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}_JointPRS_SDPRX_${pop1}_${pop2}_${GWAS_type}_${i}_repeat${rpt}" >> $job_file
fi

done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_Binary_2_final_NNLS_weights.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_Binary_2_final_NNLS_weights-$(date +%Y-%m-%d).sh


## Step2: Obtain final MIXPRS weight
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_Binary_3_MIXPRS_weight.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune
pop1=EUR

for pop2 in EAS AFR; do
for trait in T2D BrC; do

pop=${pop2}
i=1
approx=TRUE
weight_name="non_negative_linear_weights_approx${approx}"

sst_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt"

JointPRS_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_AFR_beta_EUR.txt"
JointPRS_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_AFR_beta_EAS.txt"
JointPRS_AFR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_AFR_beta_AFR.txt"
SDPRX_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_${pop}_beta_EUR.txt"
SDPRX_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_EAS_beta_EAS.txt"
SDPRX_AFR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_AFR_beta_AFR.txt"

weight_file1="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_AFR_${GWAS_type}_${i}_repeat1_${pop}_${weight_name}.txt"
weight_file2="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_AFR_${GWAS_type}_${i}_repeat2_${pop}_${weight_name}.txt"
weight_file3="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_AFR_${GWAS_type}_${i}_repeat3_${pop}_${weight_name}.txt"
weight_file4="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_AFR_${GWAS_type}_${i}_repeat4_${pop}_${weight_name}.txt"

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${pop}_MIXPRS.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_final_combine.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=${sst_file} --pop=${pop} --prs_beta_file=${JointPRS_EUR},${JointPRS_EAS},${JointPRS_AFR},${SDPRX_EUR},${SDPRX_EAS},${SDPRX_AFR} --weight_file=${weight_file1},${weight_file2},${weight_file3},${weight_file4} --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}" >> $job_file
fi

done
done


module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_Binary_3_MIXPRS_weight.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_Binary_3_MIXPRS_weight-$(date +%Y-%m-%d).sh


job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_Binary_2_MIXPRS_weight.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune
pop=EAS
pop1=EUR
pop2=EAS

for trait in CAD LuC; do

pop=${pop2}
i=1
approx=TRUE
weight_name="non_negative_linear_weights_approx${approx}"

sst_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt"

JointPRS_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_beta_EUR.txt"
JointPRS_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_beta_EAS.txt"
SDPRX_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_${pop}_beta_EUR.txt"
SDPRX_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_EAS_beta_EAS.txt"

weight_file1="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_${i}_repeat1_${pop}_${weight_name}.txt"
weight_file2="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_${i}_repeat2_${pop}_${weight_name}.txt"
weight_file3="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_${i}_repeat3_${pop}_${weight_name}.txt"
weight_file4="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_JointPRS_SDPRX_EUR_EAS_${GWAS_type}_${i}_repeat4_${pop}_${weight_name}.txt"

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/Final_weight/no_val/${trait}_${pop}_MIXPRS.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_final_combine.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=${sst_file} --pop=${pop} --prs_beta_file=${JointPRS_EUR},${JointPRS_EAS},${SDPRX_EUR},${SDPRX_EAS} --weight_file=${weight_file1},${weight_file2},${weight_file3},${weight_file4} --out_dir=result/summary_result/Final_weight/no_val --out_name=${trait}" >> $job_file
fi

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_prune_Binary_2_MIXPRS_weight.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_prune_Binary_2_MIXPRS_weight-$(date +%Y-%m-%d).sh
