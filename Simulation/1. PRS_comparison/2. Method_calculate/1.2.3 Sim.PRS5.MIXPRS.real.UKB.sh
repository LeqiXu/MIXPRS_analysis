## Step0: Obtain MIX format GWAS
## subsample PRS
library(data.table)

h2 = 0.4
rhog = 0.8

type = "prune_snplist_1"
approx_list = c("TRUE")

for (sample_size in c("25K","90K")){
for (pop in c("EAS","AFR","SAS","AMR")){
for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

for (rpt in c(1:4)){

for (approx in approx_list){

sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_", type, "_", pop,"_tune_GWAS_approx",approx,"_ratio3.00_repeat",rpt,".txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","BETA","SE","Z","P","N")]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/MIX/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_", type, "_tune_GWAS_approx",approx,"_MIX_repeat",rpt,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

}

}
}

}
}
}

## Step1: Linear combination with different snplists
# Step1: Estimate beta
for subpop in EAS AFR SAS AMR; do

job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/Final_weight_${subpop}.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

sample1=UKB

selection_criterion=NNLS
non_negative_weights=TRUE

type=prune_snplist_1
type_original="prune"
approx_list="TRUE"

for sample2 in 25K 90K; do
for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for rpt in {1..4}; do

for approx in ${approx_list}; do

# JointPRS_name
if [[ "${subpop}" == "EAS" ]]; then
  JointPRS_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_subEAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "AFR" ]]; then
  JointPRS_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_subAFR_SAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "SAS" ]]; then
  JointPRS_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_subSAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "AMR" ]]; then
  JointPRS_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_subAMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
else
  echo "Please provide a valid subpop: EUR, EAS, AFR, SAS, or AMR."
fi

for pop in EUR EAS AFR SAS AMR; do
  declare "${pop}_JointPRS_prs_beta=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/subsample/JointPRS/${JointPRS_name}_beta_${pop}.txt"
done

# SDPRX_name
sample_size="${sample2}"
EUR_SDPRX_prs_beta="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/subsample/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_sub${subpop}_${sample1}_${sample2}_SDPRX_subsample_${type}_approx${approx}_repeat${rpt}_beta_EUR.txt"
for pop in EAS AFR SAS AMR; do
  declare "${pop}_SDPRX_prs_beta=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/subsample/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_sub${pop}_${sample1}_${sample2}_SDPRX_subsample_${type}_approx${approx}_repeat${rpt}_beta_${pop}.txt"
done

prs_beta_file="${EUR_JointPRS_prs_beta},${EAS_JointPRS_prs_beta},${AFR_JointPRS_prs_beta},${SAS_JointPRS_prs_beta},${AMR_JointPRS_prs_beta},${EUR_SDPRX_prs_beta},${EAS_SDPRX_prs_beta},${AFR_SDPRX_prs_beta},${SAS_SDPRX_prs_beta},${AMR_SDPRX_prs_beta}"

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/Final_weight/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_repeat${rpt}_non_negative_linear_weights_approx${approx}.txt" ]]; then

echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/MIX/${subpop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_${type}_tune_GWAS_approx${approx}_MIX_repeat${rpt}.txt --pop=${subpop} --prs_beta_file=${prs_beta_file} --indep_approx=${approx} --selection_criterion=${selection_criterion} --non_negative_weights=${non_negative_weights} --out_dir=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/Final_weight --out_name=sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_repeat${rpt}" >> $job_file

fi

done

done
done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/Final_weight_${subpop}.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-Final_weight_${subpop}-$(date +%Y-%m-%d).sh

done


## Step2: Obtain final MIXPRS weight
for subpop in EAS AFR SAS AMR; do

job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/MIXPRS_${subpop}.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

sample1=UKB

type=prune_snplist_1
type_original="prune"
approx_list="TRUE"

for sample2 in 25K 90K; do
for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do

for approx in ${approx_list}; do

sample_size="${sample2}"

JointPRS_EUR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_EUR.txt"
JointPRS_EAS="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_EAS.txt"
JointPRS_AFR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_AFR.txt"
JointPRS_SAS="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_SAS.txt"
JointPRS_AMR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_AMR.txt"

SDPRX_EUR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_${subpop}_beta_EUR.txt"
SDPRX_EAS="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_EAS_beta_EAS.txt"
SDPRX_AFR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_AFR_beta_AFR.txt"
SDPRX_SAS="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_SAS_beta_SAS.txt"
SDPRX_AMR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${sample2}_SDPRX_real_EUR_AMR_beta_AMR.txt"

prs_file="${JointPRS_EUR},${JointPRS_EAS},${JointPRS_AFR},${JointPRS_SAS},${JointPRS_AMR},${SDPRX_EUR},${SDPRX_EAS},${SDPRX_AFR},${SDPRX_SAS},${SDPRX_AMR}"

weight_file1="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/Final_weight/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_repeat1_${subpop}_non_negative_linear_weights_approx${approx}.txt"
weight_file2="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/Final_weight/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_repeat2_${subpop}_non_negative_linear_weights_approx${approx}.txt"
weight_file3="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/Final_weight/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_repeat3_${subpop}_non_negative_linear_weights_approx${approx}.txt"
weight_file4="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/Final_weight/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_repeat4_${subpop}_non_negative_linear_weights_approx${approx}.txt"
weight_file="${weight_file1},${weight_file2},${weight_file3},${weight_file4}"

sst_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/MIX/${subpop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_MIX_real_all.txt"

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}_${subpop}_MIXPRS.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_final_combine.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=${sst_file} --pop=${subpop} --prs_beta_file=${prs_file} --weight_file=${weight_file} --indep_approx=${approx} --out_dir=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/MIXPRS --out_name=sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_SDPRX_${type}_non_negative_linear_weights_approx${approx}" >> $job_file
fi

done

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/MIXPRS_${subpop}.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-MIXPRS_${subpop}-$(date +%Y-%m-%d).sh

done