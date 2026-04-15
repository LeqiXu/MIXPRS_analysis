# Clean data
library(data.table)

trait = "CAD"

for (pop in c("EUR","EAS")){
sumstat_clean = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/", trait, "_", pop, "_inter_clean.txt"))

# SBayesRC
sumstat_SBayesRC = sumstat_clean[,c("SNP","A1","A2","MAF","BETA","SE","P","N")]
colnames(sumstat_SBayesRC) = c("SNP","A1","A2","freq","b","se","p","N")
write.table(sumstat_SBayesRC, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/SBayesRC/",trait,"_",pop,"_inter_SBayesRC.txt"), row.names=F, col.names=T, quote=F, append=F, sep = "\t")

}

# SBayesRC beta
module load miniconda
conda activate SBayesRC

trait=CAD

for trainpop in EUR EAS; do

# Variables: need to be fixed
ma_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/SBayesRC/${trait}_${trainpop}_inter_SBayesRC.txt" # GWAS summary in COJO format (the only input)
ld_folder="/gpfs/gibbs/pi/zhao/lx94/PRSmap/data/ref_data/SBayesRC/ukb${trainpop}_HM3"                    # LD reference (download from "Resources")
annot="/gpfs/gibbs/pi/zhao/lx94/PRSmap/data/ref_data/SBayesRC/annot_baseline2.2.txt"              # Functional annotation (download from "Resources")
out_prefix="/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_${trainpop}_inter_SBayesRC_HM3"   # Output prefix, e.g. "./test"
threads=4                       # Number of CPU cores

##############################################
# Code: usually don't need a change in this section
## Note: Flags were documented in the package, use ?function in R to lookup.
## We suggest to run those in multiple jobs (tasks)
export OMP_NUM_THREADS=$threads # Revise the threads

##############################################
# Tidy: optional step, tidy summary data
## "log2file=TRUE" means the messages will be redirected to a log file 
Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
                  output='${out_prefix}_tidy.ma', log2file=TRUE)"
## Best practice: read the log to check issues in your GWAS summary data.  

# Impute: optional step if your summary data doesn't cover the SNP panel
Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='$ld_folder', \
                  output='${out_prefix}_imp.ma', log2file=TRUE)"

# SBayesRC: main function for SBayesRC
Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}_sbrc', annot='$annot', log2file=TRUE)"
done


# SBayesRC score
trait=CAD
targetpop=EAS

for trainpop in EUR EAS; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_${trainpop}_inter_SBayesRC_HM3_prs_${targetpop}.sscore" ]]; then

module load PLINK/2

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${targetpop} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/${targetpop}_inter_snplist_ukbb.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/SBayesRC/${trait}_${trainpop}_inter_SBayesRC_HM3_sbrc.txt 1 2 3 \
--out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_${trainpop}_inter_SBayesRC_HM3_prs_${targetpop}

fi
done