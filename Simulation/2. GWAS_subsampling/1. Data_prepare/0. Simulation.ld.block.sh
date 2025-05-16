## 1. Divide LD blocks for each population
library(data.table)

print("Start!")

for (pop in c("EAS","AFR","SAS","AMR")){

  if (pop == "EUR"){
    blk = read.table("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_info/EUR_block.txt", header=T)
  }else if (pop %in% c("SAS", "EAS")){
    blk = read.table("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_info/EAS_block.txt", header=T)
  }else if(pop %in% c("AFR", "AMR")){
    blk = read.table("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_info/AFR_block.txt", header=T)
  }
  
  for (geno in c("100K")){
    print(paste0("Start: pop=", pop, ", geno=", geno))
    
    data = read.table(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/",pop,"/discover_validate/",pop,"_",geno,"_all.bim"))
    colnames(data) = c("chr","snpid","nonthing","position","minor","major")
    
    size = rep(0, nrow(blk))
    chrom = rep(0, nrow(blk))
    for (i in 1:nrow(blk)) {
      chr = as.numeric(sub("chr","",blk[i,1]))
      idx = which(data$position >= as.numeric(blk[i,2]) & data$position < as.numeric(blk[i,3]) & data$chr == chr)
      if (length(idx) > 0) {
        write.table(data[idx,2], file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_snplist/",geno,"/",pop,"/",i), 
                    row.names=F, col.names=F, quote=F, append=F)
      }
      size[i] = length(idx)
      chrom[i] = chr
    }
    
    write.table(size, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_info/",pop,"_",geno,"_blk_size.txt"), 
                row.names=F, col.names=F, quote=F, append=F)
    
    write.table(chrom, file=paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_info/",pop,"_",geno,"_blk_chr.txt"), 
                row.names=F, col.names=F, quote=F, append=F)
  }

print("Finish!")
}


## 2. Calculate LD for each population
geno=100K

for pop in EAS AFR SAS AMR; do
mkdir -p /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk/${geno}/${pop}

# Iterate through each file in the block-wise snplist directory
module load PLINK/1.9b_6.21-x86_64

for file in "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_snplist/${geno}/${pop}"/*; do
# Extract the block number from the filename
blk=$(basename "${file}")

plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/discover_validate/${pop}_${geno}_all \
--keep-allele-order \
--extract ${file} \
--r square \
--out /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk/${geno}/${pop}/ldblk${blk}

done
done

## 3. Write LD for each population
module load miniconda
conda activate py_env

split=discover_validate
geno=100K
fold=0

for pop in EAS AFR SAS AMR; do

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk_func/ld_write.py \
--pop=${pop} \
--split=${split} \
--geno=${geno} \
--fold=${ff}

done

## 4. Generate LD snp info for each population
module load PLINK/1

split=discover_validate
geno_type=100K
cv=all

for pop in EAS AFR SAS AMR; do

# Calculate MAF using PLINK
plink --bfile /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/${split}/${pop}_${geno}_${cv} \
--allow-extra-chr \
--freq \
--out /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk/${geno_type}/${pop}/${pop}_${geno}_maf_${cv}

if [[ ${cv} == "all" ]]; then
# Create output file and write header
echo -e "CHR\tSNP\tBP\tA1\tA2\tMAF" > /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/${split}/${geno}/ldblk_${geno_type}_${pop,,}/snpinfo_${geno_type}_hm3
# Merge MAF data with .bim data
awk 'NR==FNR{a1[$2]=$3; a2[$2]=$4; maf[$2]=$5; next; next} {print $1 "\t" $2 "\t" $4 "\t" a1[$2] "\t" a2[$2] "\t" maf[$2]}' /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/ldblk/${geno_type}/${pop}/${pop}_${geno}_maf_${cv}.frq /gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/geno_data/${pop}/${split}/${pop}_${geno}_${cv}.bim >> /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/ref_data/${split}/${geno}/ldblk_${geno_type}_${pop,,}/snpinfo_${geno_type}_hm3
fi

done