## Step0: Train subsample SBayesRC
library(data.table)

trait = "BMI"

for (pop in c("EAS","AFR")){
for (data_type in c("train","tune")){
for (rpt in c(1:4)){

subsample_sumstat_clean = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/",trait,"_full_snplist_",pop,"_",data_type,"_GWAS_approxFALSE_ratio3.00_repeat",rpt,".txt"))

# SBayesRC
subsample_sumstat_SBayesRC = subsample_sumstat_clean[,c("SNP","A1","A2","A1_Frq","BETA","SE","P","N")]
colnames(subsample_sumstat_SBayesRC) = c("SNP","A1","A2","freq","b","se","p","N")

write.table(subsample_sumstat_SBayesRC, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/SBayesRC/",trait,"_full_snplist_",pop,"_",data_type,"_SBayesRC_approxFALSE_ratio3.00_repeat",rpt,".txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}
}
}

# SBayesRC beta
module load miniconda
conda activate SBayesRC

trait=BMI
data_type=train

for trainpop in EAS AFR; do
for rpt in {1..4}; do

# Variables: need to be fixed
ma_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/SBayesRC/${trait}_full_snplist_${trainpop}_${data_type}_SBayesRC_approxFALSE_ratio3.00_repeat${rpt}.txt" # GWAS summary in COJO format (the only input)
ld_folder="/gpfs/gibbs/pi/zhao/lx94/PRSmap/data/ref_data/SBayesRC/ukb${trainpop}_HM3"                    # LD reference (download from "Resources")
annot="/gpfs/gibbs/pi/zhao/lx94/PRSmap/data/ref_data/SBayesRC/annot_baseline2.2.txt"              # Functional annotation (download from "Resources")
out_prefix="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/SBayesRC/subsample/${trait}_full_snplist_${trainpop}_${data_type}_SBayesRC_approxFALSE_ratio3.00_repeat${rpt}"   # Output prefix, e.g. "./test"
threads=4                       # Number of CPU cores

##############################################
# Code: usually don't need a change in this section
## Note: Flags were documented in the package, use ?function in R to lookup.
## We suggest to run those in multiple jobs (tasks)
export OMP_NUM_THREADS=$threads # Revise the threads

##############################################
# Tidy: optional step, tidy summary data
## "log2file=TRUE" means the messages will be redirected to a log file 
Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
                  output='${out_prefix}_tidy.ma', log2file=TRUE)"
## Best practice: read the log to check issues in your GWAS summary data.  

# Impute: optional step if your summary data doesn't cover the SNP panel
Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='$ld_folder', \
                  output='${out_prefix}_imp.ma', log2file=TRUE)"

# SBayesRC: main function for SBayesRC
Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}_sbrc', annot='$annot', log2file=TRUE)"
done
done

# Clean SBayesRC
## subsample
trait=BMI
data_type=train

for trainpop in EAS AFR; do
for rpt in {1..4}; do
cut -f1-3 /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/SBayesRC/subsample/${trait}_full_snplist_${trainpop}_${data_type}_SBayesRC_approxFALSE_ratio3.00_repeat${rpt}_sbrc.txt \
> /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/SBayesRC/subsample/${trait}_full_snplist_${trainpop}_${data_type}_SBayesRC_repeat${rpt}_sbrc.txt

done
done

## original
trait=BMI

for trainpop in EUR EAS AFR; do

cut -f1-3 /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_${trainpop}_inter_SBayesRC_HM3_sbrc.txt \
> /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_${trainpop}_inter_SBayesRC_sbrc.txt 

done

## Step1: Linear combination with full snplists
trait=BMI
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_MIXPRS_plus_${trait}.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_full
pop1=EUR
data_type=train

for pop2 in EAS AFR; do
for rpt in {1..4}; do

pop=${pop2}
approx=FALSE
weight_name="non_negative_linear_weights_approx${approx}"

if [[ ${pop2} == "EAS" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_SBayesRC_clinic_EUR_EAS_AFR_${GWAS_type}_repeat${rpt}_${pop}_${weight_name}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_full_snplist_${pop}_tune_MIX_approx${approx}_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_repeat${rpt}_beta_EUR.txt,result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_repeat${rpt}_beta_EAS.txt,result/summary_result/JointPRS/${trait}_EUR_subEAS_AFR_JointPRS_${GWAS_type}_repeat${rpt}_beta_AFR.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_repeat${rpt}_beta_${pop1}.txt,result/summary_result/SDPRX/${trait}_EUR_EAS_SDPRX_${GWAS_type}_repeat${rpt}_beta_EAS.txt,result/summary_result/SDPRX/${trait}_EUR_AFR_SDPRX_${GWAS_type}_repeat${rpt}_beta_AFR.txt,result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_EUR_inter_SBayesRC_sbrc.txt,result/summary_result/clinicalPRS/${trait}/SBayesRC/subsample/${trait}_full_snplist_EAS_${data_type}_SBayesRC_repeat${rpt}_sbrc.txt,result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_AFR_inter_SBayesRC_sbrc.txt,/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_META_PGS005199_inter_plink_${pop}.txt --indep_approx=${approx} --out_dir=result/summary_result/clinicalPRS/${trait}/MIXPRS_plus --out_name=${trait}_JointPRS_SDPRX_SBayesRC_clinic_EUR_EAS_AFR_${GWAS_type}_repeat${rpt}" >> $job_file
fi

fi

if [[ ${pop2} == "AFR" ]]; then

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_SBayesRC_clinic_EUR_EAS_AFR_${GWAS_type}_repeat${rpt}_${pop}_${weight_name}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_full_snplist_${pop}_tune_MIX_approx${approx}_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_repeat${rpt}_beta_EUR.txt,result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_repeat${rpt}_beta_EAS.txt,result/summary_result/JointPRS/${trait}_EUR_EAS_subAFR_JointPRS_${GWAS_type}_repeat${rpt}_beta_AFR.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_repeat${rpt}_beta_${pop1}.txt,result/summary_result/SDPRX/${trait}_EUR_EAS_SDPRX_${GWAS_type}_repeat${rpt}_beta_EAS.txt,result/summary_result/SDPRX/${trait}_EUR_AFR_SDPRX_${GWAS_type}_repeat${rpt}_beta_AFR.txt,result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_EUR_inter_SBayesRC_sbrc.txt,result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_EAS_inter_SBayesRC_sbrc.txt,result/summary_result/clinicalPRS/${trait}/SBayesRC/subsample/${trait}_full_snplist_AFR_${data_type}_SBayesRC_repeat${rpt}_sbrc.txt,/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_META_PGS005199_inter_plink_${pop}.txt --indep_approx=${approx} --out_dir=result/summary_result/clinicalPRS/${trait}/MIXPRS_plus --out_name=${trait}_JointPRS_SDPRX_SBayesRC_clinic_EUR_EAS_AFR_${GWAS_type}_repeat${rpt}" >> $job_file
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

GWAS_type=subsample_full
pop1=EUR

for pop2 in EAS AFR; do

pop=${pop2}
i=1
approx=FALSE
weight_name="non_negative_linear_weights_approx${approx}"

sst_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt"

JointPRS_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_AFR_beta_EUR.txt"
JointPRS_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_AFR_beta_EAS.txt"
JointPRS_AFR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_AFR_beta_AFR.txt"
SDPRX_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_${pop}_beta_EUR.txt"
SDPRX_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_EAS_beta_EAS.txt"
SDPRX_AFR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_AFR_beta_AFR.txt"
SBayesRC_EUR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_EUR_inter_SBayesRC_sbrc.txt"
SBayesRC_EAS="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_EAS_inter_SBayesRC_sbrc.txt"
SBayesRC_AFR="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_AFR_inter_SBayesRC_sbrc.txt"
ClinicalPRS_META="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_META_PGS005199_inter_plink_${pop}.txt"

weight_file1="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_SBayesRC_clinic_EUR_EAS_AFR_${GWAS_type}_repeat1_${pop}_${weight_name}.txt"
weight_file2="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_SBayesRC_clinic_EUR_EAS_AFR_${GWAS_type}_repeat2_${pop}_${weight_name}.txt"
weight_file3="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_SBayesRC_clinic_EUR_EAS_AFR_${GWAS_type}_repeat3_${pop}_${weight_name}.txt"
weight_file4="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_SBayesRC_clinic_EUR_EAS_AFR_${GWAS_type}_repeat4_${pop}_${weight_name}.txt"

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_MIXPRS_plus_ClinicalPRS_SBayesRC_${pop}_MIXPRS.txt"

if [[ ! -e ${output_file} ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=MIXPRS_plus_BMI
#SBATCH --output=/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/MIXPRS_plus_BMI_beta.txt

module load miniconda
conda activate py_env
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_final_combine.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=${sst_file} \
--pop=${pop} \
--prs_beta_file=${JointPRS_EUR},${JointPRS_EAS},${JointPRS_AFR},${SDPRX_EUR},${SDPRX_EAS},${SDPRX_AFR},${SBayesRC_EUR},${SBayesRC_EAS},${SBayesRC_AFR},${ClinicalPRS_META} \
--weight_file=${weight_file1},${weight_file2},${weight_file3},${weight_file4} \
--out_dir=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus \
--out_name=${trait}_MIXPRS_plus_ClinicalPRS_SBayesRC

EOT
fi
done

## Step 3: Calculate PRS
trait=BMI

for pop in EAS AFR; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_MIXPRS_plus_ClinicalPRS_SBayesRC_prs_${pop}.sscore" ]]; then

module load PLINK/2

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/${pop}_inter_snplist_ukbb.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_MIXPRS_plus_ClinicalPRS_SBayesRC_${pop}_MIXPRS.txt \
--out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_MIXPRS_plus_ClinicalPRS_SBayesRC_prs_${pop}

fi
done