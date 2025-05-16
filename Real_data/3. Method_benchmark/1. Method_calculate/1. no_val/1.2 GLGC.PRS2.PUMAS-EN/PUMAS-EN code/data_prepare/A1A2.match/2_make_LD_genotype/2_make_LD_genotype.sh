#!/bin/bash

# --- Population & Study Setup ---
declare -A study_pops
study_pops[GLGC]="0 1 2 3 4"
study_pops[PAGE]="0 1 2"
study_pops[BBJ]="0 1"

pops=("EUR" "EAS" "AFR" "AMR" "SAS")
popnames=("eur" "eas" "afr" "amr" "sas")

input_base="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/A1A2.match/1_extract_samples&snps"
output_base="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/LD_ref/A1A2.match/2_make_LD_genotype"
geno_base="/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/geno_data"
plink_path="/vast/palmer/apps/avx2/software/PLINK/1.9b_6.21-x86_64/plink"

# --- Create Output Directories ---
for study in "${!study_pops[@]}"; do
    for idx in ${study_pops[$study]}; do
        pop=${pops[$idx]}
        mkdir -p "${output_base}/${study}/${pop}/qced"
    done
done

# --- Run PLINK for Each Study/Population ---
for study in "${!study_pops[@]}"; do
    for idx in ${study_pops[$study]}; do
        pop=${pops[$idx]}
        popname=${popnames[$idx]}

        ${plink_path} --bfile "${geno_base}/${pop}" \
                      --keep "${input_base}/${study}/tmp/${pop}_LD500.samples.txt" \
                      --extract "${input_base}/${study}/PRScs_${popname}_1kg.snplist" \
                      --a1-allele "${input_base}/${study}/PRScs_${popname}_1kg_A1.txt" 2 1 \
                      --fill-missing-a2 \
                      --make-bed \
                      --out "${output_base}/${study}/${pop}/qced/1KG.LD500.PRSCS.1kg.A1A2match"

        ${plink_path} --bfile "${output_base}/${study}/${pop}/qced/1KG.LD500.PRSCS.1kg.A1A2match" \
                      --freq \
                      --keep-allele-order \
                      --out "${output_base}/${study}/${pop}/qced/1KG.LD500.PRSCS.1kg.A1A2match.freq"
    done
done
