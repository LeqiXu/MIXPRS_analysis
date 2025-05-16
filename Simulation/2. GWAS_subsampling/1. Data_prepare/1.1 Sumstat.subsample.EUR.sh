## Step1: Obtain clean GWAS for MIX
library(data.table)

h2 = 0.4
rhog = 0.8
pop = "EUR"
sample_size = "ukbb"

for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover_validate/clean/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_clean_real_all.txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","BETA","SE","Z","P","N")]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/MIX/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_MIX_real_all.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}

## Step2: Subsample GWAS by MIX
## full snplist
h2=0.4
rhog=0.8
pop=EUR
sample_size=ukbb

for sim_i in {1..5}; do
for p in 0.1 0.01 0.001 5e-04; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_${pop}_tune_GWAS_approxFALSE_ratio3.00_repeat4.txt" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist
#SBATCH --output=out_GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist.txt

module load miniconda
conda activate py_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=data/sim_data/summary_data/discover_validate/MIX/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_MIX_real_all.txt \
--pop=${pop} \
--indep_approx=FALSE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/sim_data/summary_data/subsample_data/clean \
--out_name=${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist

EOT
fi

done
done

## prune snplist
h2=0.4
rhog=0.8
pop=EUR
sample_size=ukbb
i=1

for sim_i in {1..5}; do
for p in 0.1 0.01 0.001 5e-04; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_prune_snplist_${i}_${pop}_tune_GWAS_approxFALSE_ratio3.00_repeat4.txt" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_prune_snplist_${i}
#SBATCH --output=out_GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_prune_snplist_${i}.txt

module load miniconda
conda activate py_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=data/sim_data/summary_data/discover_validate/MIX/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_MIX_real_all.txt \
--pop=${pop} \
--prune_snplist=/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/${pop}_prune_pval1_r20.5_wc250_${i}.snplist \
--indep_approx=TRUE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/sim_data/summary_data/subsample_data/clean \
--out_name=${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_prune_snplist_${i}

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=data/sim_data/summary_data/discover_validate/MIX/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_MIX_real_all.txt \
--pop=${pop} \
--prune_snplist=/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/${pop}_prune_pval1_r20.5_wc250_${i}.snplist \
--indep_approx=FALSE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/sim_data/summary_data/subsample_data/clean \
--out_name=${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_prune_snplist_${i}

EOT
fi

done
done

