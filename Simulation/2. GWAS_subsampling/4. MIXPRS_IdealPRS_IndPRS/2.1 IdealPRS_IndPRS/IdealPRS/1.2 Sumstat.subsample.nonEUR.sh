## Step2: Subsample GWAS by MIX
## full snplist with 100K
h2=0.4
rhog=0.8
sample_size=100K

for pop in EAS AFR SAS AMR; do
for sim_i in {1..5}; do
for p in 0.1 0.01 0.001 5e-04; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_100K_${pop}_tune_GWAS_approxFALSE_ratio3.00_repeat4.txt" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --mem=150G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_100K
#SBATCH --output=out_GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_100K.txt

module load miniconda
conda activate py_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/ref_data/PRScsx/100K \
--sst_file=data/sim_data/summary_data/discover_validate/MIX/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_MIX_real_all.txt \
--pop=${pop} \
--indep_approx=FALSE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/sim_data/summary_data/subsample_data/clean \
--out_name=${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_100K

EOT
fi

done
done
done

## Step2: Subsample GWAS by MIX
## full snplist with 1kg
h2=0.4
rhog=0.8
sample_size=100K

for pop in EAS AFR SAS AMR; do
for sim_i in {1..5}; do
for p in 0.1 0.01 0.001 5e-04; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_1kg_${pop}_tune_GWAS_approxFALSE_ratio3.00_repeat4.txt" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --mem=150G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_1kg
#SBATCH --output=out_GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_1kg.txt

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
--out_name=${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_1kg

EOT
fi

done
done
done