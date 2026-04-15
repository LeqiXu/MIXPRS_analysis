## Step0: Prune snplist
library(data.table)

# load weight file
trait = "BMI"

BMI_files <- list(
  Single = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_Single_PGS005198_plink.txt",
  META   = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_META_PGS005199_plink.txt",
  EUR    = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_EUR_PGS005203_plink.txt",
  EAS    = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_EAS_PGS005202_plink.txt",
  AFR    = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_AFR_PGS005200_plink.txt",
  SAS    = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_SAS_PGS005204_plink.txt",
  AMR    = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_AMR_PGS005201_plink.txt"
)

prs_to_inter_prune_plink <- function(file, pop) {
  prs <- fread(file)
  
  ## inter file
  inter_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/",pop,"_inter_snplist.txt"), header = FALSE)
  inter_prs = prs[which(prs$SNP %in% inter_snplist$V1),]
  write.table(inter_prs,paste0(sub("\\_plink.txt$", "", file), "_inter_plink_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
  
  ## inter prune file
  prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/",pop,"/prune_pval1_r20.5_wc250_1.prune.in"), header = FALSE)
  inter_prune_prs = inter_prs[which(inter_prs$SNP %in% prune_snplist$V1),]
  write.table(inter_prune_prs,paste0(sub("\\_plink.txt$", "", file), "_inter_plink_prune_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

  return(list(inter_prs,inter_prune_prs))
}

BMI_PRS_plink_list_EAS <- lapply(BMI_files, prs_to_inter_prune_plink, pop = "EAS")
BMI_PRS_plink_list_AFR <- lapply(BMI_files, prs_to_inter_prune_plink, pop = "AFR")

## Step1: Linear combination with different snplists
trait=BMI
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_MIXPRS_plus_${trait}.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune
pop1=EUR

for pop2 in EAS AFR; do
for rpt in {1..4}; do

pop=${pop2}
i=1
approx=TRUE
weight_name="non_negative_linear_weights_approx${approx}"

if [[ ${pop2} == "EAS" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_EUR_EAS_AFR_${GWAS_type}_${i}_repeat${rpt}_${pop}_${weight_name}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/xz674/MIXPRS/main/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_prune_snplist_${i}_${pop}_tune_MIX_approx${approx}_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EUR.txt,result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EAS.txt,result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_AFR.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EUR.txt,result/summary_result/SDPRX/${trait}_${pop1}_EAS_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EAS.txt,result/summary_result/SDPRX/${trait}_${pop1}_AFR_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_AFR.txt,/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_META_PGS005199_inter_plink_prune_${pop}.txt --out_dir=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus --out_name=${trait}_JointPRS_SDPRX_clinic_EUR_EAS_AFR_${GWAS_type}_${i}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "AFR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_EUR_EAS_AFR_${GWAS_type}_${i}_repeat${rpt}_${pop}_${weight_name}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/xz674/MIXPRS/main/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_prune_snplist_${i}_${pop}_tune_MIX_approx${approx}_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EUR.txt,result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EAS.txt,result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_AFR.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EUR.txt,result/summary_result/SDPRX/${trait}_${pop1}_EAS_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_EAS.txt,result/summary_result/SDPRX/${trait}_${pop1}_AFR_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_AFR.txt,/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_META_PGS005199_inter_plink_prune_${pop}.txt --out_dir=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus --out_name=${trait}_JointPRS_SDPRX_clinic_EUR_EAS_AFR_${GWAS_type}_${i}_repeat${rpt}" >> $job_file
fi

fi

done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_MIXPRS_plus_${trait}.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_MIXPRS_plus_${trait}-$(date +%Y-%m-%d).sh


## Step2: Obtain final MIXPRS weight
trait=BMI
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_MIXPRS_plus_${trait}_2.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune
pop1=EUR

for pop2 in EAS AFR; do

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
PRSCSX_META="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_META_PGS005199_inter_plink_${pop}.txt"

weight_file1="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_EUR_EAS_AFR_${GWAS_type}_${i}_repeat1_${pop}_${weight_name}.txt"
weight_file2="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_EUR_EAS_AFR_${GWAS_type}_${i}_repeat2_${pop}_${weight_name}.txt"
weight_file3="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_EUR_EAS_AFR_${GWAS_type}_${i}_repeat3_${pop}_${weight_name}.txt"
weight_file4="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_EUR_EAS_AFR_${GWAS_type}_${i}_repeat4_${pop}_${weight_name}.txt"

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_plus_${pop}_MIXPRS.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_final_combine.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=${sst_file} --pop=${pop} --prs_beta_file=${JointPRS_EUR},${JointPRS_EAS},${JointPRS_AFR},${SDPRX_EUR},${SDPRX_EAS},${SDPRX_AFR},${PRSCSX_META} --weight_file=${weight_file1},${weight_file2},${weight_file3},${weight_file4} --out_dir=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus --out_name=${trait}_plus" >> $job_file
fi

done


module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_MIXPRS_plus_${trait}_2.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_MIXPRS_plus_${trait}_2-$(date +%Y-%m-%d).sh


## Step 3: Calculate PRS
trait="BMI"

for pop in EAS AFR; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_MIXPRS_plus_prs_${pop}.sscore" ]]; then

module load PLINK/2

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/${pop}_inter_snplist_ukbb.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_plus_${pop}_MIXPRS.txt \
--out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_MIXPRS_plus_prs_${pop}

fi
done