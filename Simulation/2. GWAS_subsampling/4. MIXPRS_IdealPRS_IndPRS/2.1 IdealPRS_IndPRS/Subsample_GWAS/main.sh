## full snplist Non-EUR
h2="0.4"
rhog="0.8"
approx="FALSE"

mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/clean

for pop in AFR SAS AMR; do # EUR EAS 
if [[ ${pop} == "EUR" ]]; then
sample_size="ukbb"
mem="50G"
else
sample_size="100K"
mem="20G"
fi
for sim_i in {1..5}; do
for p in 0.1 0.01 0.001 5e-04; do
sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day
#SBATCH --requeue
#SBATCH --mem=${mem}
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist
#SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/0923fullsnp_2/out_GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist.txt

module load miniconda
conda activate py_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/xz674/MIXPRS/main/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/ref_data/discover_validate/100K \
--sst_file=data/sim_data/summary_data/discover_validate/MIX/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_MIX_real_all.txt \
--pop=${pop} \
--indep_approx=${approx} \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/clean \
--out_name=${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist

EOT

done
done
done


# ## full snplist EUR
# h2="0.4"
# rhog="0.8"
# sample_size="ukbb"

# mkdir -p /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/clean

# for pop in EUR; do
# for sim_i in {1..5}; do
# for p in 0.1 0.01 0.001 5e-04; do

# sbatch <<EOT
# #!/bin/bash
# #SBATCH --partition=scavenge,day
# #SBATCH --requeue
# #SBATCH --mem=50G
# #SBATCH --cpus-per-task=1
# #SBATCH --ntasks=1 --nodes=1
# #SBATCH --time=1-00:00:00
# #SBATCH --mail-type=ALL
# #SBATCH --job-name=GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist
# #SBATCH --output=/gpfs/gibbs/pi/zhao/xz674/logs/subsample_GWAS/0923fullsnp/out_GWAS_subsample_generate_${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist.txt

# module load miniconda
# conda activate py_env

# cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

# python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
# --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/ukbb \
# --sst_file=data/sim_data/summary_data/discover_validate/MIX/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_MIX_real_all.txt \
# --pop=${pop} \
# --indep_approx=FALSE \
# --train_tune_ratio=3 \
# --repeat=4 \
# --out_dir=/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/clean \
# --out_name=${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample_size}_full_snplist

# EOT

# done
# done
# done