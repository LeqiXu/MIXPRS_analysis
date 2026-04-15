## Step0: Prune snplist
trait = "CAD"
pop = "EAS"

CAD_PRS_inter_plink = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/CAD/CAD_hg19_META_PGS003725_inter_plink.txt"))
prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/",pop,"/prune_pval1_r20.5_wc250_1.prune.in"), header = FALSE)

CAD_PRS_inter_plink_prune = CAD_PRS_inter_plink[which(CAD_PRS_inter_plink$SNP %in% prune_snplist$V1),]

write.table(CAD_PRS_inter_plink_prune,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/CAD/CAD_hg19_META_PGS003725_inter_plink_prune.txt"),quote=F,sep='\t',row.names=F,col.names=T)


## Step1: Linear combination with different snplists
trait=CAD
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_MIXPRS_plus_${trait}.txt"
> $job_file  # Empty the job file if it already exists

GWAS_type=subsample_prune
pop=EAS
pop1=EUR
pop2=EAS

for rpt in {1..4}; do

i=1
approx=TRUE
weight_name="non_negative_linear_weights_approx${approx}"

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_${pop1}_${pop2}_${GWAS_type}_${i}_repeat${rpt}_${pop}_${weight_name}.txt"

if [[ ! -e ${output_file} ]]; then
    echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_linear_weight.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --sst_file=data/summary_data/subsample/MIX/${trait}_prune_snplist_${i}_${pop}_tune_MIX_approx${approx}_ratio3.00_repeat${rpt}.txt --pop=${pop} --prs_beta_file=result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop1}.txt,result/summary_result/JointPRS/${trait}_${pop1}_${pop2}_JointPRS_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop2}.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop1}.txt,result/summary_result/SDPRX/${trait}_${pop1}_${pop2}_SDPRX_${GWAS_type}_${i}_approx${approx}_repeat${rpt}_beta_${pop2}.txt,/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/CAD/CAD_hg19_META_PGS003725_inter_plink_prune.txt --out_dir=result/summary_result/clinicalPRS/${trait}/MIXPRS_plus --out_name=${trait}_JointPRS_SDPRX_clinic_EUR_EAS_${GWAS_type}_${i}_repeat${rpt}" >> $job_file
fi

done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/subsample_MIXPRS_plus_${trait}.txt --partition=scavenge,day --requeue --mem=30G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-subsample_MIXPRS_plus_${trait}-$(date +%Y-%m-%d).sh


## Step2: Obtain final MIXPRS weight
trait=CAD

GWAS_type=subsample_prune
pop=EAS
pop1=EUR
pop2=EAS

pop=${pop2}
i=1
approx=TRUE
weight_name="non_negative_linear_weights_approx${approx}"

sst_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt"

JointPRS_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_beta_EUR.txt"
JointPRS_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/JointPRS/${trait}_JointPRS_EUR_EAS_beta_EAS.txt"
SDPRX_EUR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_${pop}_beta_EUR.txt"
SDPRX_EAS="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/result/summary_result/Final_weight/no_val/SDPRX/${trait}_SDPRX_EUR_EAS_beta_EAS.txt"
GPSmulti_CAD="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/CAD/CAD_hg19_META_PGS003725_inter_plink.txt"

weight_file1="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_EUR_EAS_${GWAS_type}_${i}_repeat1_${pop}_${weight_name}.txt"
weight_file2="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_EUR_EAS_${GWAS_type}_${i}_repeat2_${pop}_${weight_name}.txt"
weight_file3="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_EUR_EAS_${GWAS_type}_${i}_repeat3_${pop}_${weight_name}.txt"
weight_file4="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_JointPRS_SDPRX_clinic_EUR_EAS_${GWAS_type}_${i}_repeat4_${pop}_${weight_name}.txt"

output_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_plus_${pop}_MIXPRS.txt"

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
#SBATCH --job-name=MIXPRS_plus_CAD
#SBATCH --output=/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/real_data/method/Final_weight/MIXPRS_plus_CAD_beta.txt

module load miniconda
conda activate py_env
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_final_combine.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=${sst_file} \
--pop=${pop} \
--prs_beta_file=${JointPRS_EUR},${JointPRS_EAS},${SDPRX_EUR},${SDPRX_EAS},${GPSmulti_CAD} \
--weight_file=${weight_file1},${weight_file2},${weight_file3},${weight_file4} \
--out_dir=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus \
--out_name=${trait}_plus

EOT
fi

## Step 3: Calculate PRS
trait=CAD
pop=EAS

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_MIXPRS_plus_prs_${pop}.sscore" ]]; then

module load PLINK/2

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/${pop}_inter_snplist_ukbb.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/MIXPRS_plus/${trait}_plus_${pop}_MIXPRS.txt \
--out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_MIXPRS_plus_prs_${pop}

fi