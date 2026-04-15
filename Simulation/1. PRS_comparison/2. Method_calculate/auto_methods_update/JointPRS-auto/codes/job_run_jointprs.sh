#!/bin/bash
#SBATCH --job-name=JointPRS_array%j
#SBATCH --output=/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/JointPRS/logs/25.10.5/JointPRS_array%j.log
#SBATCH --requeue
#SBATCH --mem=32G
#SBATCH -p scavenge,day,pi_zhao
#SBATCH -t 24:00:00
#SBATCH -c 1

module load miniconda
source /vast/palmer/apps/avx2/software/miniconda/24.9.2/etc/profile.d/conda.sh
conda activate /home/yd357/pi_paths/pi_zhao/.conda/envs/JointPRS 

CODE_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/code/4_methods/JointPRS/codes"
INPUT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/data/sim_data/summary_data/discover_validate/JointPRS"
REF_DIR="/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg"
BIM_PREFIX="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/All/All_test"
RESULT_DIR="/home/yd357/pi_paths/pi_zhao/MIXPRS_Revise/Project1_15K_80K/result/add_benchmark/25.10.4/sim_data/JointPRS"
JOINTPRS_PATH="/gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py"

IFS=',' read -r SIM_I P SAMPLE_SIZE CHR < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${CODE_DIR}/jointprs_params.txt)
RHOG="0.8"

# The EUR population always uses the "80K" GWAS data.
declare -a sst_files
sst_files+=("${INPUT_DIR}/EUR_sim${SIM_I}_p${P}_rho${RHOG}_80K_subcol.txt")
for pop in "EAS" "AFR" "SAS" "AMR"; do
    sst_files+=("${INPUT_DIR}/${pop}_sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE}_subcol.txt")
done

for file in "${sst_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required input file not found: $file"
        exit 1
    fi
done

sst_file_list=$(IFS=, ; echo "${sst_files[*]}")

if [[ "${SAMPLE_SIZE}" == "15K" ]]; then
    sample_size_eur="311600"
    sample_size_other="15000"
elif [[ "${SAMPLE_SIZE}" == "80K" ]]; then
    sample_size_eur="311600"
    sample_size_other="80000"
else
    echo "Error: Unknown sample size specified: ${SAMPLE_SIZE}"
    exit 1
fi
n_gwas_string="${sample_size_eur},${sample_size_other},${sample_size_other},${sample_size_other},${sample_size_other}"

mkdir -p "${RESULT_DIR}"

echo "========================================================"
echo "RUNNING JointPRS with parameters:"
echo "  Simulation ID : ${SIM_I}"
echo "  P-value       : ${P}"
echo "  Sample Scenario : ${SAMPLE_SIZE}"
echo "  Chromosome    : ${CHR}"
echo "  GWAS N values : ${n_gwas_string}"
echo "  Input SST list: ${sst_file_list}"
echo "========================================================"

python $JOINTPRS_PATH \
    --ref_dir=$REF_DIR \
    --bim_prefix=$BIM_PREFIX \
    --sst_file=${sst_file_list} \
    --rho_cons=1,1,1,1,1 \
    --n_gwas=${n_gwas_string} \
    --chrom=${CHR} \
    --pop=EUR,EAS,AFR,SAS,AMR \
    --out_dir=$RESULT_DIR \
    --out_name=sim${SIM_I}_p${P}_rho${RHOG}_${SAMPLE_SIZE}_JointPRS

echo "JointPRS analysis complete for sim${SIM_I}, p=${P}, sample_size=${SAMPLE_SIZE}, chr=${CHR}"