# Clean data
library(data.table)

# load bim map
map_1kg <- fread("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg/snpinfo_mult_1kg_hm3")
map_1kg <- map_1kg[, .(SNP, CHR, BP)]
setnames(map_1kg, c("SNP", "CHR", "BP"), c("SNP", "CHR", "POS"))

# load weight file
trait = "BMI"

BMI_files <- list(
  Single = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_Single_PGS005198.txt",
  META   = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_META_PGS005199.txt",
  EUR    = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_EUR_PGS005203.txt",
  EAS    = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_EAS_PGS005202.txt",
  AFR    = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_AFR_PGS005200.txt",
  SAS    = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_SAS_PGS005204.txt",
  AMR    = "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/BMI/BMI_hg19_AMR_PGS005201.txt"
)

prs_to_plink <- function(file, map_1kg) {
  prs <- fread(file)
  # keep + rename columns
  prs <- prs[, .(chr_name, chr_position, effect_allele, other_allele, effect_weight)]
  setnames(prs, c("CHR", "POS", "A1", "A2", "BETA"))
  prs[, POS := as.integer(POS)]
  
  # merge with 1KG map by CHR, POS
  prs_map <- merge(prs, map_1kg, by = c("CHR", "POS"))
  
  # keep PLINK-score columns
  prs_plink <- prs_map[, .(SNP, A1, BETA)]
  write.table(prs_plink,paste0(sub("\\.txt$", "", file), "_plink.txt"),quote=F,sep='\t',row.names=F,col.names=T)

  return(prs_plink)
}

BMI_PRS_plink_list <- lapply(BMI_files, prs_to_plink, map_1kg = map_1kg)

# ClinicalPRS_BMI
trait=BMI

for type in Single META EUR EAS AFR SAS AMR; do

score_file=$(ls /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/${trait}/${trait}_hg19_${type}_PGS*_plink.txt)

for pop in EAS AFR; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_ClinicalPRS_${type}_prs_${pop}.sscore" ]]; then

module load PLINK/2

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/${pop}_inter_snplist_ukbb.txt \
--score ${score_file} \
--out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_ClinicalPRS_${type}_prs_${pop}

fi
done
done