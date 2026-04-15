# Clean data
library(data.table)
trait = "CAD"
pop = "EAS"

# load weight file
CAD_PRS = fread("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/CAD/CAD_hg19_META_PGS003725.txt")
CAD_PRS = CAD_PRS[,c("chr_name","chr_position","effect_allele","other_allele","effect_weight")]
colnames(CAD_PRS) = c("CHR","POS","A1","A2","BETA")
CAD_PRS$POS = as.integer(CAD_PRS$POS)

# load bim map
map_1kg = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg/ldblk_1kg_",tolower(pop),"/snpinfo_1kg_hm3"))
map_1kg = map_1kg[,c("SNP","CHR","BP")]
colnames(map_1kg) = c("SNP","CHR","POS")

# prs_score file
CAD_PRS_map = merge(CAD_PRS,map_1kg,by=c("CHR","POS"))
CAD_PRS_plink = CAD_PRS_map[,c("SNP","A1","BETA")]

write.table(CAD_PRS_plink,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/CAD/CAD_hg19_META_PGS003725_plink.txt"),quote=F,sep='\t',row.names=F,col.names=T)

# inter prs_score file
inter_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/",pop,"_inter_snplist.txt"), header = FALSE)
CAD_PRS_inter_plink = CAD_PRS_plink[which(CAD_PRS_plink$SNP %in% inter_snplist$V1),]
write.table(CAD_PRS_inter_plink,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/CAD/CAD_hg19_META_PGS003725_inter_plink.txt"),quote=F,sep='\t',row.names=F,col.names=T)

# ClinicalPRS_CAD
trait=CAD
pop=EAS

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_ClinicalPRS_prs_${pop}.sscore" ]]; then

module load PLINK/2

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
--double-id \
--threads 1 \
--extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/snplist_data/${pop}_inter_snplist_ukbb.txt \
--score /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/PGS_data/CAD/CAD_hg19_META_PGS003725_inter_plink.txt \
--out /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/summary_result/clinicalPRS/${trait}/UKB_${trait}_ClinicalPRS_prs_${pop}

fi
