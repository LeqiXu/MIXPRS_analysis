#!/bin/bash
#SBATCH --mem=4G
#SBATCH --partition=scavenge,day
#SBATCH -t 1:00:00
#SBATCH -c 1

# --- Configuration ---
BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v3/concatenate"
FILTER_STEP_DIR="${BASE_DIR}/codes/array_job/step4_filtered_pumas_eval"
SEPARATED_BASE_DIR="${BASE_DIR}/results/separate"
MODEL_LIST_FULL="${SEPARATED_BASE_DIR}/model_list.txt"
SUB_BETA_SRC_DIR_BASE="${SEPARATED_BASE_DIR}/sub_beta"
FULL_BETA_SRC_DIR_BASE="${SEPARATED_BASE_DIR}/full_beta"
FILTERED_LIST_DIR="${FILTER_STEP_DIR}/filtered_model_lists"
JOB_PARAMS_FILE="${BASE_DIR}/codes/array_job/step4_pumas_eval/job_params_step4_pumas_eval.txt"
K_FOLDS=4

mkdir -p "${FILTERED_LIST_DIR}" "${FILTER_STEP_DIR}/logs"

# --- Process Each Trait/Population ---
while IFS=',' read -r TARGET_TRAIT TARGET_POP; do
    [ -z "$TARGET_TRAIT" ] || [ -z "$TARGET_POP" ] && continue

    TRAIT_TAG="${TARGET_TRAIT}_${TARGET_POP}_inter"
    SUB_BETA_SRC_DIR="${SUB_BETA_SRC_DIR_BASE}/${TRAIT_TAG}"
    FULL_BETA_SRC_DIR="${FULL_BETA_SRC_DIR_BASE}/${TRAIT_TAG}"
    ZERO_MODELS_FILE="${FILTERED_LIST_DIR}/${TARGET_TRAIT}_${TARGET_POP}.zero_models.txt"
    FILTERED_MODELS_FILE="${FILTERED_LIST_DIR}/${TARGET_TRAIT}_${TARGET_POP}.filtered_models.txt"

    > "${ZERO_MODELS_FILE}"
    zero_model_count=0

    # --- Check Each Model ---
    while IFS= read -r model_name; do
        [ -z "$model_name" ] && continue
        is_all_zero=0

        files_to_check=("${FULL_BETA_SRC_DIR}/${TRAIT_TAG}.${model_name}.ite1.txt")
        for i in $(seq 1 $K_FOLDS); do
            files_to_check+=("${SUB_BETA_SRC_DIR}/${TRAIT_TAG}.${model_name}.ite${i}.txt")
        done

        for file in "${files_to_check[@]}"; do
            [ ! -f "$file" ] && continue
            awk 'BEGIN{all_zero=1} $3!=0 {all_zero=0; exit 1} END{exit !all_zero}' "$file"
            if [ $? -eq 0 ]; then
                is_all_zero=1
                break
            fi
        done

        [ "$is_all_zero" -eq 1 ] && echo "$model_name" >> "${ZERO_MODELS_FILE}" && ((zero_model_count++))
    done < "${MODEL_LIST_FULL}"

    # --- Create Filtered List ---
    if [ -s "${ZERO_MODELS_FILE}" ]; then
        grep -v -x -f "${ZERO_MODELS_FILE}" "${MODEL_LIST_FULL}" > "${FILTERED_MODELS_FILE}"
    else
        cp "${MODEL_LIST_FULL}" "${FILTERED_MODELS_FILE}"
    fi

done < "${JOB_PARAMS_FILE}"
