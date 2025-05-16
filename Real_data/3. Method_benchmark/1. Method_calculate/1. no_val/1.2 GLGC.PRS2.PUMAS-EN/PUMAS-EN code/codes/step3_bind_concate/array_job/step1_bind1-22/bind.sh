#!/bin/bash
#SBATCH --mem=4G
#SBATCH --partition=scavenge,day,week,pi_zhao
#SBATCH -t 4:00:00
#SBATCH -c 1

# --- Configuration ---
POP_INFO="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/PUMAS-EN/data/basic_pop_info/basic_pop_info.csv"
CHROMSOMES=( $(seq 1 22) )
MAX_ITER=4
CONCAT_BASE_OUTPUT_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v3/concatenate/results"

# Lassosum Parameters
LASSOSUM_S_PARAMS=("0.2" "0.5" "0.9")
LASSOSUM_LAMBDA_PARAMS=("0.005" "0.01")

get_trait_ancestry() {
    local line="$1"
    trait=$(echo "$line" | awk -F, '{gsub(/^["[:space:]]+|["[:space:]]+$/, "", $1); print $1}')
    ancestry=$(echo "$line" | awk -F, '{gsub(/^["[:space:]]+|["[:space:]]+$/, "", $2); print $2}')
}

# --- Process Lassosum ---
METHOD_NAME="Lassosum"
LASSOSUM_SOURCE_BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/Lassosum"

tail -n +2 "${POP_INFO}" | while IFS= read -r line; do
    get_trait_ancestry "$line"
    [[ -z "$trait" || -z "$ancestry" ]] && continue
    trait_anc_file_prefix="${trait}_${ancestry}"
    final_output_path="${CONCAT_BASE_OUTPUT_DIR}/${METHOD_NAME}/${trait_anc_file_prefix}"

    # --- PUMAS=T ---
    for iter in $(seq 1 $MAX_ITER); do
        for s in "${LASSOSUM_S_PARAMS[@]}"; do
            for lambda in "${LASSOSUM_LAMBDA_PARAMS[@]}"; do
                case "$lambda" in
                    "0.005") tag="lambda_0p005" ;;
                    "0.01")  tag="lambda_0p01" ;;
                    *) continue ;;
                esac
                src_dir="${LASSOSUM_SOURCE_BASE_DIR}/results/PRSweights/PUMAS-ite${iter}/${trait_anc_file_prefix}"
                out_file="${final_output_path}/${trait_anc_file_prefix}_PUMASite${iter}_ALLCHR_sX${s}_${tag}.lassosum.betas.txt"
                files_to_cat=(); all_exist=true
                for chr in "${CHROMSOMES[@]}"; do
                    f="${src_dir}/${trait_anc_file_prefix}_chr${chr}_sX${s}_${tag}.lassosum.betas.txt"
                    [[ -f "$f" ]] && files_to_cat+=("$f") || { all_exist=false; break; }
                done
                if $all_exist && [ ${#files_to_cat[@]} -eq 22 ]; then
                    mkdir -p "${final_output_path}"
                    head -n 1 "${files_to_cat[0]}" > "$out_file"
                    for f in "${files_to_cat[@]}"; do tail -n +2 "$f" >> "$out_file"; done
                fi
            done
        done
    done

    # --- PUMAS=F ---
    for s in "${LASSOSUM_S_PARAMS[@]}"; do
        for lambda in "${LASSOSUM_LAMBDA_PARAMS[@]}"; do
            case "$lambda" in
                "0.005") tag="lambda_0p005" ;;
                "0.01")  tag="lambda_0p01" ;;
                *) continue ;;
            esac
            src_dir="${LASSOSUM_SOURCE_BASE_DIR}/results/PRSweights/Regular/${trait_anc_file_prefix}"
            out_file="${final_output_path}/${trait_anc_file_prefix}_Regular_ALLCHR_sX${s}_${tag}.lassosum.betas.txt"
            files_to_cat=(); all_exist=true
            for chr in "${CHROMSOMES[@]}"; do
                f="${src_dir}/${trait_anc_file_prefix}_chr${chr}_sX${s}_${tag}.lassosum.betas.txt"
                [[ -f "$f" ]] && files_to_cat+=("$f") || { all_exist=false; break; }
            done
            if $all_exist && [ ${#files_to_cat[@]} -eq 22 ]; then
                mkdir -p "${final_output_path}"
                head -n 1 "${files_to_cat[0]}" > "$out_file"
                for f in "${files_to_cat[@]}"; do tail -n +2 "$f" >> "$out_file"; done
            fi
        done
    done
done


# --- Process LDPred2 ---
METHOD_NAME="LDPred2"
echo "[INFO] Processing ${METHOD_NAME} files..."
LDPRED2_SOURCE_BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/LDPred2"
LDPRED2_SOURCE_RESULTS_DIR="${LDPRED2_SOURCE_BASE_DIR}/results/PRSweights"
LDPRED2_P_OPTIONS=("0.001" "0.01" "0.1")
LDPRED2_H2_RATIOS=("0.1" "0.3")
LDPRED2_SPARSE_OPTIONS=("T" "F")
LDPRED2_SHRINK_CORR_AUTO="0.5"

tail -n +2 "${POP_INFO}" | while IFS= read -r line; do
    get_trait_ancestry "$line"
    if [ -z "$trait" ] || [ -z "$ancestry" ]; then continue; fi
    trait_anc_dir_name="${trait}_${ancestry}"
    trait_anc_file_prefix="${trait}_${ancestry}"

    final_output_path="${CONCAT_BASE_OUTPUT_DIR}/${METHOD_NAME}/${trait_anc_dir_name}"
    echo "  [${METHOD_NAME}] Trait: ${trait}, Ancestry: ${ancestry} (Output SubDir: ${final_output_path})"

    # PUMAS=T
    source_pumas_t_segment="PUMAS/ite" 
    for iter in $(seq 1 $MAX_ITER); do
        pumas_filename_tag_for_output="_PUMASite${iter}"
        current_source_pumas_dir="${LDPRED2_SOURCE_RESULTS_DIR}/${source_pumas_t_segment}${iter}/${trait_anc_dir_name}"
        # Grid mode
        for p_param in "${LDPRED2_P_OPTIONS[@]}"; do
            for h2_ratio in "${LDPRED2_H2_RATIOS[@]}"; do
                for sparse_opt in "${LDPRED2_SPARSE_OPTIONS[@]}"; do
                    target_filename_base="${trait_anc_file_prefix}${pumas_filename_tag_for_output}_ALLCHR_grid_p${p_param}_h${h2_ratio}_sp${sparse_opt}.ldpred2.txt"
                    target_concat_file="${final_output_path}/${target_filename_base}"
                    files_to_cat=()
                    all_files_exist=true
                    for chr_num in "${CHROMSOMES[@]}"; do
                        src_file="${current_source_pumas_dir}/${trait_anc_file_prefix}_chr${chr_num}_grid_p${p_param}_h${h2_ratio}_sp${sparse_opt}.ldpred2.txt"
                        if [ -f "${src_file}" ]; then files_to_cat+=("${src_file}"); else all_files_exist=false; echo "    [${METHOD_NAME} WARNING] PUMAS=T grid, iter=${iter}, p=${p_param}, h=${h2_ratio}, sp=${sparse_opt}: Source file missing: ${src_file}"; break; fi
                    done
                    if $all_files_exist && [ ${#files_to_cat[@]} -eq 22 ]; then
                        echo "    [${METHOD_NAME}] Concatenating PUMAS=T grid, iter=${iter}, params... to ${target_concat_file}"
                        mkdir -p "${final_output_path}"; head -n 1 "${files_to_cat[0]}" > "${target_concat_file}"; 
                        for f_idx in $(seq 0 $((${#files_to_cat[@]} - 1)) ); do tail -n +2 "${files_to_cat[$f_idx]}" >> "${target_concat_file}"; done
                    fi
                done 
            done 
        done 
        # Auto mode
        target_filename_base="${trait_anc_file_prefix}${pumas_filename_tag_for_output}_ALLCHR_auto_sc${LDPRED2_SHRINK_CORR_AUTO}.ldpred2.txt"
        target_concat_file="${final_output_path}/${target_filename_base}"
        files_to_cat=()
        all_files_exist=true
        for chr_num in "${CHROMSOMES[@]}"; do
            src_file="${current_source_pumas_dir}/${trait_anc_file_prefix}_chr${chr_num}_auto_sc${LDPRED2_SHRINK_CORR_AUTO}.ldpred2.txt"
            if [ -f "${src_file}" ]; then files_to_cat+=("${src_file}"); else all_files_exist=false; echo "    [${METHOD_NAME} WARNING] PUMAS=T auto, iter=${iter}: Source file missing: ${src_file}"; break; fi
        done
        if $all_files_exist && [ ${#files_to_cat[@]} -eq 22 ]; then
            echo "    [${METHOD_NAME}] Concatenating PUMAS=T auto, iter=${iter} to ${target_concat_file}"
            mkdir -p "${final_output_path}"; head -n 1 "${files_to_cat[0]}" > "${target_concat_file}"; 
            for f_idx in $(seq 0 $((${#files_to_cat[@]} - 1)) ); do tail -n +2 "${files_to_cat[$f_idx]}" >> "${target_concat_file}"; done
        fi
    done # iter

    # PUMAS=F
    pumas_filename_tag_for_output="_Regular"
    current_source_regular_dir="${LDPRED2_SOURCE_RESULTS_DIR}/Regular/${trait_anc_dir_name}"
    # Grid mode
    for p_param in "${LDPRED2_P_OPTIONS[@]}"; do
        for h2_ratio in "${LDPRED2_H2_RATIOS[@]}"; do
            for sparse_opt in "${LDPRED2_SPARSE_OPTIONS[@]}"; do
                target_filename_base="${trait_anc_file_prefix}${pumas_filename_tag_for_output}_ALLCHR_grid_p${p_param}_h${h2_ratio}_sp${sparse_opt}.ldpred2.txt"
                target_concat_file="${final_output_path}/${target_filename_base}"
                files_to_cat=()
                all_files_exist=true
                for chr_num in "${CHROMSOMES[@]}"; do
                    src_file="${current_source_regular_dir}/${trait_anc_file_prefix}_chr${chr_num}_grid_p${p_param}_h${h2_ratio}_sp${sparse_opt}.ldpred2.txt"
                    if [ -f "${src_file}" ]; then files_to_cat+=("${src_file}"); else all_files_exist=false; echo "    [${METHOD_NAME} WARNING] PUMAS=F grid, params...: Source file missing: ${src_file}"; break; fi
                done
                if $all_files_exist && [ ${#files_to_cat[@]} -eq 22 ]; then
                    echo "    [${METHOD_NAME}] Concatenating PUMAS=F grid, params... to ${target_concat_file}"
                    mkdir -p "${final_output_path}"; head -n 1 "${files_to_cat[0]}" > "${target_concat_file}"; 
                    for f_idx in $(seq 0 $((${#files_to_cat[@]} - 1)) ); do tail -n +2 "${files_to_cat[$f_idx]}" >> "${target_concat_file}"; done
                fi
            done 
        done 
    done 
    # Auto mode
    target_filename_base="${trait_anc_file_prefix}${pumas_filename_tag_for_output}_ALLCHR_auto_sc${LDPRED2_SHRINK_CORR_AUTO}.ldpred2.txt"
    target_concat_file="${final_output_path}/${target_filename_base}"
    files_to_cat=()
    all_files_exist=true
    for chr_num in "${CHROMSOMES[@]}"; do
        src_file="${current_source_regular_dir}/${trait_anc_file_prefix}_chr${chr_num}_auto_sc${LDPRED2_SHRINK_CORR_AUTO}.ldpred2.txt"
        if [ -f "${src_file}" ]; then files_to_cat+=("${src_file}"); else all_files_exist=false; echo "    [${METHOD_NAME} WARNING] PUMAS=F auto: Source file missing: ${src_file}"; break; fi
    done
    if $all_files_exist && [ ${#files_to_cat[@]} -eq 22 ]; then
        echo "    [${METHOD_NAME}] Concatenating PUMAS=F auto to ${target_concat_file}"
        mkdir -p "${final_output_path}"; head -n 1 "${files_to_cat[0]}" > "${target_concat_file}"; 
        for f_idx in $(seq 0 $((${#files_to_cat[@]} - 1)) ); do tail -n +2 "${files_to_cat[$f_idx]}" >> "${target_concat_file}"; done
    fi
done < <(tail -n +2 "${POP_INFO}")
echo "[INFO] ${METHOD_NAME} processing finished."
echo "-----------------------------------------------------"


# --- Process PRS-CS ---
METHOD_NAME="PRS-CS"
echo "[INFO] Processing ${METHOD_NAME} files..."
PRSCS_SOURCE_BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/PRS-CS"
PRSCS_PHI_OPTIONS=("1e-4" "1e-6" "auto") 
PRSCS_HEADER="CHR\tSNP\tPOS\tA1\tA2\tBETA"

tail -n +2 "${POP_INFO}" | while IFS= read -r line; do
    get_trait_ancestry "$line"
    if [ -z "$trait" ] || [ -z "$ancestry" ]; then continue; fi
    trait_anc_dir_name="${trait}_${ancestry}"
    trait_anc_file_prefix="${trait}_${ancestry}"

    final_output_path="${CONCAT_BASE_OUTPUT_DIR}/${METHOD_NAME}/${trait_anc_dir_name}"
    echo "  [${METHOD_NAME}] Trait: ${trait}, Ancestry: ${ancestry} (Output SubDir: ${final_output_path})"
    
    source_pumas_mode_prefix="PUMAS-ite"
    # PUMAS=T
    for iter in $(seq 1 $MAX_ITER); do
        pumas_filename_tag_for_output="_PUMASite${iter}"
        current_source_pumas_dir="${PRSCS_SOURCE_BASE_DIR}/results/PRSweights/${source_pumas_mode_prefix}${iter}/${trait_anc_dir_name}"
        for phi_val in "${PRSCS_PHI_OPTIONS[@]}"; do
            phi_fmt_fname=""
            case "$phi_val" in
                "1e-4") phi_fmt_fname="1e-04" ;;
                "1e-6") phi_fmt_fname="1e-06" ;;
                "auto") phi_fmt_fname="auto" ;;
            esac
            
            source_file_base_name_part="${trait_anc_file_prefix}_pst_eff_a1_b0.5_phi${phi_fmt_fname}"
            target_filename_base="${trait_anc_file_prefix}${pumas_filename_tag_for_output}_pst_eff_a1_b0.5_phi${phi_fmt_fname}_ALLCHR.txt"
            target_concat_file="${final_output_path}/${target_filename_base}"
            
            source_input_file_prefix_for_chr="${source_file_base_name_part}_chr"

            files_to_cat=()
            all_files_exist=true
            for chr_num in "${CHROMSOMES[@]}"; do
                src_file="${current_source_pumas_dir}/${source_input_file_prefix_for_chr}${chr_num}.txt"
                if [ -f "${src_file}" ]; then files_to_cat+=("${src_file}"); else all_files_exist=false; echo "    [${METHOD_NAME} WARNING] PUMAS=T, iter=${iter}, phi=${phi_fmt_fname}: Source file missing: ${src_file}"; break; fi
            done

            if $all_files_exist && [ ${#files_to_cat[@]} -eq 22 ]; then
                echo "    [${METHOD_NAME}] Concatenating PUMAS=T, iter=${iter}, phi=${phi_fmt_fname} to ${target_concat_file}"
                mkdir -p "${final_output_path}"
                echo -e "${PRSCS_HEADER}" > "${target_concat_file}" 
                for f in "${files_to_cat[@]}"; do
                    cat "${f}" >> "${target_concat_file}" 
                done
            fi
        done # phi_val
    done # iter

    # PUMAS=F
    pumas_filename_tag_for_output="_Regular"
    current_source_regular_dir="${PRSCS_SOURCE_BASE_DIR}/results/PRSweights/Regular/${trait_anc_dir_name}"
    for phi_val in "${PRSCS_PHI_OPTIONS[@]}"; do
        phi_fmt_fname=""
        case "$phi_val" in
            "1e-4") phi_fmt_fname="1e-04" ;;
            "1e-6") phi_fmt_fname="1e-06" ;;
            "auto") phi_fmt_fname="auto" ;;
        esac

        source_file_base_name_part="${trait_anc_file_prefix}_pst_eff_a1_b0.5_phi${phi_fmt_fname}"
        target_filename_base="${trait_anc_file_prefix}${pumas_filename_tag_for_output}_pst_eff_a1_b0.5_phi${phi_fmt_fname}_ALLCHR.txt"
        target_concat_file="${final_output_path}/${target_filename_base}"

        source_input_file_prefix_for_chr="${source_file_base_name_part}_chr"
        
        files_to_cat=()
        all_files_exist=true
        for chr_num in "${CHROMSOMES[@]}"; do
            src_file="${current_source_regular_dir}/${source_input_file_prefix_for_chr}${chr_num}.txt"
            if [ -f "${src_file}" ]; then files_to_cat+=("${src_file}"); else all_files_exist=false; echo "    [${METHOD_NAME} WARNING] PUMAS=F, phi=${phi_fmt_fname}: Source file missing: ${src_file}"; break; fi
        done

        if $all_files_exist && [ ${#files_to_cat[@]} -eq 22 ]; then
            echo "    [${METHOD_NAME}] Concatenating PUMAS=F, phi=${phi_fmt_fname} to ${target_concat_file}"
            mkdir -p "${final_output_path}"
            echo -e "${PRSCS_HEADER}" > "${target_concat_file}" 
            for f in "${files_to_cat[@]}"; do
                cat "${f}" >> "${target_concat_file}" 
            done
        fi
    done # phi_val
done < <(tail -n +2 "${POP_INFO}")
echo "[INFO] ${METHOD_NAME} processing finished."
echo "-----------------------------------------------------"


# --- Process SBLUP ---
METHOD_NAME="SBLUP"
echo "[INFO] Processing ${METHOD_NAME} files..."
SBLUP_SOURCE_BASE_DIR="/home/yd357/pi_paths/pi_zhao/GWAS_Subsample_JointPRS/v2_PUMAS-EN/single-pop/SBLUP"

tail -n +2 "${POP_INFO}" | while IFS= read -r line; do
    get_trait_ancestry "$line"
    if [ -z "$trait" ] || [ -z "$ancestry" ]; then continue; fi
    trait_anc_dir_name="${trait}_${ancestry}"
    trait_anc_file_prefix="${trait}_${ancestry}" 

    final_output_path="${CONCAT_BASE_OUTPUT_DIR}/${METHOD_NAME}/${trait_anc_dir_name}"
    echo "  [${METHOD_NAME}] Trait: ${trait}, Ancestry: ${ancestry} (Output SubDir: ${final_output_path})"
    
    source_pumas_mode_prefix="PUMAS-ite"
    # PUMAS=T
    for iter in $(seq 1 $MAX_ITER); do
        pumas_filename_tag_for_output="_PUMASite${iter}"
        current_source_pumas_dir="${SBLUP_SOURCE_BASE_DIR}/results/PRSweights/${source_pumas_mode_prefix}${iter}/${trait_anc_dir_name}"
        
        target_filename_base="${trait_anc_file_prefix}${pumas_filename_tag_for_output}_ALLCHR.sblup.cojo" 
        target_concat_file="${final_output_path}/${target_filename_base}"
        
        files_to_cat=()
        all_files_exist=true
        for chr_num in "${CHROMSOMES[@]}"; do
            src_file="${current_source_pumas_dir}/chr${chr_num}.sblup.cojo" 
            if [ -f "${src_file}" ]; then files_to_cat+=("${src_file}"); else all_files_exist=false; echo "    [${METHOD_NAME} WARNING] PUMAS=T, iter=${iter}: Source file missing: ${src_file}"; break; fi
        done

        if $all_files_exist && [ ${#files_to_cat[@]} -eq 22 ]; then
            echo "    [${METHOD_NAME}] Concatenating PUMAS=T, iter=${iter} to ${target_concat_file}"
            mkdir -p "${final_output_path}"
            head -n 1 "${files_to_cat[0]}" > "${target_concat_file}"
            for f_idx in $(seq 0 $((${#files_to_cat[@]} - 1)) ); do tail -n +2 "${files_to_cat[$f_idx]}" >> "${target_concat_file}"; done
        fi
    done # iter

    # PUMAS=F
    pumas_filename_tag_for_output="_Regular"
    current_source_regular_dir="${SBLUP_SOURCE_BASE_DIR}/results/PRSweights/Regular/${trait_anc_dir_name}"
    
    target_filename_base="${trait_anc_file_prefix}${pumas_filename_tag_for_output}_ALLCHR.sblup.cojo" 
    target_concat_file="${final_output_path}/${target_filename_base}"

    files_to_cat=()
    all_files_exist=true
    for chr_num in "${CHROMSOMES[@]}"; do
        src_file="${current_source_regular_dir}/chr${chr_num}.sblup.cojo" 
        if [ -f "${src_file}" ]; then files_to_cat+=("${src_file}"); else all_files_exist=false; echo "    [${METHOD_NAME} WARNING] PUMAS=F: Source file missing: ${src_file}"; break; fi
    done

    if $all_files_exist && [ ${#files_to_cat[@]} -eq 22 ]; then
        echo "    [${METHOD_NAME}] Concatenating PUMAS=F to ${target_concat_file}"
        mkdir -p "${final_output_path}"
        head -n 1 "${files_to_cat[0]}" > "${target_concat_file}"
        for f_idx in $(seq 0 $((${#files_to_cat[@]} - 1)) ); do tail -n +2 "${files_to_cat[$f_idx]}" >> "${target_concat_file}"; done
    fi
done < <(tail -n +2 "${POP_INFO}")
echo "[INFO] ${METHOD_NAME} processing finished."
echo "-----------------------------------------------------"
echo "All PRS results concatenation processes finished."