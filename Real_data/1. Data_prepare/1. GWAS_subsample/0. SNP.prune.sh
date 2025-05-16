## Step1: Pruning with different snplist
## prune_1
for i in 1; do
for pval in 1; do
for r2 in 0.5; do
for wc in 250; do
for pop in EUR EAS AFR SAS AMR; do

module load PLINK/2

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/geno_data/${pop} \
--indep-pairwise ${wc} 5 ${r2} \
--out /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/${pop}/prune_pval${pval}_r2${r2}_wc${wc}_${i}

done
done
done
done
done

## prune_i (i >=2)
library(data.table)

prune_i=10

for (pop in c("EUR","EAS","AFR","SAS","AMR")){

all_snplist_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/geno_data/",pop,".bim"), header = FALSE)
all_snplist = all_snplist_table$V2

for (i in c(1:(prune_i - 1))){
prune_snplist_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/",pop,"/prune_pval1_r20.5_wc250_",i,".prune.in"), header = FALSE)
prune_snplist = prune_snplist_table$V1

if (i == 1){

filtered_snplist = setdiff(all_snplist, prune_snplist)

} else{

filtered_snplist = setdiff(filtered_snplist, prune_snplist)

}
}

filtered_snplist_table = data.table(SNP = filtered_snplist)

write.table(filtered_snplist_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/",pop,"/prune_pval1_r20.5_wc250_",prune_i,".preferred"),quote=F,sep='\t',row.names=F,col.names=T)
}

prune_i=10

for i in ${prune_i}; do
for pval in 1; do
for r2 in 0.5; do
for wc in 250; do
for pop in EUR EAS AFR SAS AMR; do

module load PLINK/2

plink2 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/1000g_phase3_data/geno_data/${pop} \
--indep-pairwise ${wc} 5 ${r2} \
--indep-preferred /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/${pop}/prune_pval1_r20.5_wc250_${i}.preferred \
--out /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/${pop}/prune_pval${pval}_r2${r2}_wc${wc}_${i}

done
done
done
done
done

## Step2: Clean prune snplist
library(data.table)

for (pop in c("EUR","EAS","AFR","SAS","AMR")){

for (i in c(1:10)){
prune_snplist_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/",pop,"/prune_pval1_r20.5_wc250_",i,".prune.in"), header = FALSE)
prune_snplist = prune_snplist_table$V1


snplist_table = data.table(SNP = prune_snplist)

write.table(snplist_table,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/",pop,"_prune_pval1_r20.5_wc250_",i,".snplist"),quote=F,sep='\t',row.names=F,col.names=F)

}
}

## Step4: Check the union of the first four snp list
library(data.table)

for (pop in c("EUR","EAS","AFR","SAS","AMR")){

prune_snplist_union = c()

for (i in c(1:4)){
prune_snplist_table = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/",pop,"_prune_pval1_r20.5_wc250_",i,".snplist"), header = FALSE)
prune_snplist = prune_snplist_table$V1

prune_snplist_union = c(prune_snplist_union, prune_snplist)

}

print(paste0(pop,": first four prune snplist total number: ", length(unique(prune_snplist_union))))
}