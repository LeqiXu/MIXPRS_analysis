## Step2: Subsample GWAS by MIX
## full snplist with ukbb
h2=0.4
rhog=0.8
pop=EUR
sample_size=ukbb

for sim_i in {1..5}; do
for p in 0.1 0.01 0.001 5e-04; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_ukbb_${pop}_tune_GWAS_approxFALSE_ratio3.00_repeat4.txt" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --mem=150G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_ukbb
#SBATCH --output=out_GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_ukbb.txt

module load miniconda
conda activate py_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/ukbb \
--sst_file=data/sim_data/summary_data/discover_validate/MIX/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_MIX_real_all.txt \
--pop=${pop} \
--indep_approx=FALSE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/sim_data/summary_data/subsample_data/clean \
--out_name=${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist_ukbb

EOT
fi

done
done