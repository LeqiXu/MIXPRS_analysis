## Step1: Obtain clean GWAS for MIX
library(data.table)

for (trait in c("T2D","BrC")){
for (pop in c("EAS","AFR")){

sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/", trait, "_", pop, "_inter_clean.txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","BETA","SE","Z","P","N")]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/MIX/",trait, "_", pop, "_inter_MIX.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}

library(data.table)

for (trait in c("CAD","LuC")){
for (pop in c("EAS")){

sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/", trait, "_", pop, "_inter_clean.txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","BETA","SE","Z","P","N")]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/MIX/",trait, "_", pop, "_inter_MIX.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}


## Step2: Subsample GWAS by MIX
## full snplist
for trait in T2D BrC; do
for pop in EAS AFR; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/${trait}_full_snplist_${pop}_tune_GWAS_approxFALSE_ratio3.00_repeat4.txt" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${trait}_full_snplist
#SBATCH --output=out_GWAS_subsample_generate_${trait}_full_snplist.txt

module load miniconda
conda activate py_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt \
--pop=${pop} \
--indep_approx=FALSE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/summary_data/subsample/clean \
--out_name=${trait}_full_snplist

EOT
fi

done
done

for trait in CAD LuC; do
for pop in EAS; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/${trait}_full_snplist_${pop}_tune_GWAS_approxFALSE_ratio3.00_repeat4.txt" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${trait}_full_snplist
#SBATCH --output=out_GWAS_subsample_generate_${trait}_full_snplist.txt

module load miniconda
conda activate py_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt \
--pop=${pop} \
--indep_approx=FALSE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/summary_data/subsample/clean \
--out_name=${trait}_full_snplist

EOT
fi

done
done

## prune snplist
for trait in T2D BrC; do
for pop in EAS AFR; do
for i in {1..4}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/${trait}_prune_snplist_${i}_${pop}_tune_GWAS_approxTRUE_ratio3.00_repeat4.txt" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${trait}_prune_snplist_${i}
#SBATCH --output=out_GWAS_subsample_generate_${trait}_prune_snplist_${i}.txt

module load miniconda
conda activate py_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt \
--pop=${pop} \
--prune_snplist=/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/${pop}_prune_pval1_r20.5_wc250_${i}.snplist \
--indep_approx=FALSE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/summary_data/subsample/clean \
--out_name=${trait}_prune_snplist_${i}

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt \
--pop=${pop} \
--prune_snplist=/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/${pop}_prune_pval1_r20.5_wc250_${i}.snplist \
--indep_approx=TRUE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/summary_data/subsample/clean \
--out_name=${trait}_prune_snplist_${i}

EOT
fi

done
done
done

for trait in CAD LuC; do
for pop in EAS; do
for i in {1..4}; do

if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/subsample/clean/${trait}_prune_snplist_${i}_${pop}_tune_GWAS_approxTRUE_ratio3.00_repeat4.txt" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=GWAS_subsample_generate_${trait}_prune_snplist_${i}
#SBATCH --output=out_GWAS_subsample_generate_${trait}_prune_snplist_${i}.txt

module load miniconda
conda activate py_env

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt \
--pop=${pop} \
--prune_snplist=/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/${pop}_prune_pval1_r20.5_wc250_${i}.snplist \
--indep_approx=FALSE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/summary_data/subsample/clean \
--out_name=${trait}_prune_snplist_${i}

python /gpfs/gibbs/pi/zhao/lx94/SWIFT/method/MIXPRS/MIX_subsample2.py \
--ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg \
--sst_file=data/summary_data/MIX/${trait}_${pop}_inter_MIX.txt \
--pop=${pop} \
--prune_snplist=/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/${pop}_prune_pval1_r20.5_wc250_${i}.snplist \
--indep_approx=TRUE \
--train_tune_ratio=3 \
--repeat=4 \
--out_dir=data/summary_data/subsample/clean \
--out_name=${trait}_prune_snplist_${i}

EOT
fi

done
done
done

## Step3: Obtain GWAS with the corresponding snplist
library(data.table)

for (trait in c("T2D","BrC")){
for (pop in c("EUR","EAS","AFR")){

sumstat_data_full_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/", trait, "_", pop, "_inter_clean.txt"))

write.table(sumstat_data_full_snplist,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/clean/",trait, "_", pop, "_inter_clean_full_snplist.txt"),quote=F,sep='\t',row.names=F,col.names=T)

for (i in c(1:4)){
for (pop2 in c("EAS","AFR")){

prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/",pop2,"_prune_pval1_r20.5_wc250_",i,".snplist"), header = FALSE)
sumstat_data_prune_snplist = sumstat_data_full_snplist[which(sumstat_data_full_snplist$SNP %in% prune_snplist$V1),]

write.table(sumstat_data_full_snplist,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/clean/",trait, "_", pop, "_inter_clean_",pop2,"_prune_snplist_",i,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}

}
}

library(data.table)

for (trait in c("CAD","LuC")){
for (pop in c("EUR")){

sumstat_data_full_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/summary_data/clean/", trait, "_", pop, "_inter_clean.txt"))

write.table(sumstat_data_full_snplist,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/clean/",trait, "_", pop, "_inter_clean_full_snplist.txt"),quote=F,sep='\t',row.names=F,col.names=T)

for (i in c(1:4)){
pop2 = "EAS"

prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/",pop2,"_prune_pval1_r20.5_wc250_",i,".snplist"), header = FALSE)
sumstat_data_prune_snplist = sumstat_data_full_snplist[which(sumstat_data_full_snplist$SNP %in% prune_snplist$V1),]

write.table(sumstat_data_full_snplist,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/summary_data/clean/",trait, "_", pop, "_inter_clean_",pop2,"_prune_snplist_",i,".txt"),quote=F,sep='\t',row.names=F,col.names=T)
}

}
}
