#!/bin/bash
#SBATCH --mem=32G
#SBATCH -p scavenge,day,week,pi_zhao
#SBATCH -t 0:30:00
#SBATCH -c 1

BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/PRS-CS"
CODE_DIR="${BASE_DIR}/codes/base"
PARAM_DIR="${CODE_DIR}"
job_param_file="${PARAM_DIR}/format_params.txt"

# Read parameters for this specific task
IFS=',' read -r TRAIT ANCESTRY PUMAS ITER < <(sed -n "${SLURM_ARRAY_TASK_ID}p" ${job_param_file})

TRAIT=$(echo $TRAIT | tr -d '"')
ANCESTRY=$(echo $ANCESTRY | tr -d '"')
PUMAS=$(echo $PUMAS | tr -d '"')
ITER=$(echo $ITER | tr -d '"') # Will be empty for PUMAS=F

echo "Parameters:"
echo "  TRAIT:    ${TRAIT}"
echo "  ANCESTRY: ${ANCESTRY}"
echo "  PUMAS:    ${PUMAS}"
echo "  ITER:     ${ITER:-'N/A'}"

# Define paths based on PUMAS mode
if [ "$PUMAS" == "T" ]; then
    # PUMAS paths
    GWAS_IN_PATH="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/subsampling/PUMA-ensemble.subsampling/EN/${TRAIT}_${ANCESTRY}_inter"
    GWAS_RAW_FILE="${GWAS_IN_PATH}/${TRAIT}_${ANCESTRY}_inter.gwas.omnibus.ite${ITER}.txt"
    SUBSET_GWAS_PATH="${BASE_DIR}/results/sub_gwas/PUMAS-ite${ITER}"
    GWAS_FORMATTED_OUT="${SUBSET_GWAS_PATH}/${TRAIT}_${ANCESTRY}.prscs.txt"
else
    # Regular GWAS paths (PUMAS=F)
    GWAS_IN_PATH="/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean"
    GWAS_RAW_FILE="${GWAS_IN_PATH}/${TRAIT}_${ANCESTRY}_inter_clean.txt"
    SUBSET_GWAS_PATH="${BASE_DIR}/results/sub_gwas/Regular"
    GWAS_FORMATTED_OUT="${SUBSET_GWAS_PATH}/${TRAIT}_${ANCESTRY}.prscs.txt"
fi

echo "Paths:"
echo "  Raw Input File: ${GWAS_RAW_FILE}"
echo "  Output Dir:     ${SUBSET_GWAS_PATH}"
echo "  Formatted File: ${GWAS_FORMATTED_OUT}"

mkdir -p "${SUBSET_GWAS_PATH}"
echo "Ensured output directory exists: ${SUBSET_GWAS_PATH}"

# Check if the raw input GWAS file exists
if [ ! -f "${GWAS_RAW_FILE}" ]; then
    echo "ERROR: Raw input GWAS file does not exist: ${GWAS_RAW_FILE}"
    echo "Skipping formatting for this task."
    exit 1
fi

# Check if a valid formatted file already exists
format_needed=1
if [ -f "${GWAS_FORMATTED_OUT}" ]; then
    echo "Formatted file already exists: ${GWAS_FORMATTED_OUT}"
    if [ $(awk 'END {print NR}' "${GWAS_FORMATTED_OUT}") -ge 6 ]; then
        INVALID_LINES=$(awk 'NR>1 && NF!=5 {count++} END {print count+0}' "${GWAS_FORMATTED_OUT}")
        if [ "${INVALID_LINES}" -eq 0 ]; then
            echo "Existing file appears valid (>=6 lines, all data lines have 5 fields). Skipping formatting."
            format_needed=0
        else
            echo "WARNING: Existing file has ${INVALID_LINES} invalid lines (NF!=5). Will re-format."
        fi
    else
        echo "WARNING: Existing file has < 6 lines. Will re-format."
    fi
fi

# Perform formatting if needed
if [ $format_needed -eq 1 ]; then
    echo "Formatting GWAS file: ${GWAS_RAW_FILE}"

    echo -e "SNP\tA1\tA2\tBETA\tP" > "${GWAS_FORMATTED_OUT}.tmp"

    awk 'BEGIN {FS="\t|,| "; OFS="\t"; snp=0; a1=0; a2=0; beta=0; p=0}
         NR==1 {
             for(i=1;i<=NF;i++) {
                 if($i=="SNP" || $i=="snp") snp=i;
                 if($i=="A1" || $i=="a1") a1=i;
                 if($i=="A2" || $i=="a2") a2=i;
                 if($i=="BETA" || $i=="beta") beta=i;
                 if($i=="P" || $i=="p") p=i;
             }
             # Check if all required columns were found
             if (snp==0 || a1==0 || a2==0 || beta==0 || p==0) {
                 print "ERROR: Missing required columns in header of " FILENAME > "/dev/stderr"; exit 1
             }
             next # Skip header row in output appending step
         }
         # Process data rows: print only if all required fields are present and non-empty
         (snp && a1 && a2 && beta && p && $snp!="" && $a1!="" && $a2!="" && $beta!="" && $p!="") {
             print $snp, $a1, $a2, $beta, $p
         }' "${GWAS_RAW_FILE}" | awk 'NF==5' >> "${GWAS_FORMATTED_OUT}.tmp"

    if [ $? -eq 0 ] && [ -f "${GWAS_FORMATTED_OUT}.tmp" ]; then
        mv "${GWAS_FORMATTED_OUT}.tmp" "${GWAS_FORMATTED_OUT}"
        echo "Successfully formatted GWAS file. Output saved to: ${GWAS_FORMATTED_OUT}"
        final_lines=$(wc -l < "${GWAS_FORMATTED_OUT}")
        echo "Formatted file contains ${final_lines} lines (including header)."
    else
        echo "ERROR: Formatting failed for ${GWAS_RAW_FILE}. Check awk errors or temporary file."
        rm -f "${GWAS_FORMATTED_OUT}.tmp"
        exit 1
    fi

else
    echo "Using previously formatted file: ${GWAS_FORMATTED_OUT}"
fi

echo "Formatting job task ${SLURM_ARRAY_TASK_ID} finished."
date