## 1. Obtain popcorn format GWAS
library(data.table)

h2 = 0.4
rhog = 0.8
pop = "EUR"
sample_size = "UKB"

for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/summary_data/",pop,"/discover/clean/",pop,"_sim",sim_i,"_p",p,"_rho",rhog,"_",sample_size,"_clean_real.txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","MAF","N","BETA","SE","Z")]
colnames(sumstat_data) = c("SNP","a1","a2","af","N","beta","SE","Z")

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/popcorn/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_popcorn_real_all.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}

library(data.table)

h2 = 0.4
rhog = 0.8


for (sample_size in c("25K","90K")){
for (pop in c("EAS","AFR","SAS","AMR")){
for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/summary_data/",pop,"/discover_validate/clean/",pop,"_sim",sim_i,"_p",p,"_rho",rhog,"_",sample_size,"_clean_real.txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","MAF","N","BETA","SE","Z")]
colnames(sumstat_data) = c("SNP","a1","a2","af","N","beta","SE","Z")

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/popcorn/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_popcorn_real_all.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}
}
}

## 2. Obtain genetic correlation
module load miniconda
conda activate popcorn

cd /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/Popcorn

h2=0.4
rhog=0.8

pop1=EUR
sample1=UKB

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for pop2 in EAS AFR SAS AMR; do
for sample2 in 25K 90K; do

popcorn fit -v 0 \
--cfile ./ref/${pop1}_${pop2}_all_gen_eff.cscore \
--gen_effect \
--sfile1 /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/popcorn/${pop1}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_popcorn_real_all.txt \
--sfile2 /gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/popcorn/${pop2}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_popcorn_real_all.txt \
/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/popcorn/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_${pop2}_${sample1}_${sample2}_popcorn_real_corr.txt

done
done
done
done

## 2. Check genetic correlation
h2=0.4
rhog=0.8

pop1=EUR
sample1=UKB

for sample2 in 25K 90K; do
for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do

line="## sample=${sample1}_${sample2} p=${p} rhog=${rhog} i=${sim_i}"

# For each target pop2
for pop2 in EAS AFR SAS AMR; do
      
# Construct the file name from which to extract pge
file=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/popcorn/\
sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_${pop2}_${sample1}_${sample2}_popcorn_real_corr.txt

pge_val=$(grep '^pge' "${file}" | awk '{printf "%.2f", $2}')
      
# Add it to our line under the heading "EUR_POP:"
line="${line}  ${pop1}_${pop2}:${pge_val}"

done

echo "${line}"

done
done
done


## sample=UKB_25K p=0.001 rhog=0.8 i=1  EUR_EAS:0.80  EUR_AFR:0.81  EUR_SAS:0.75  EUR_AMR:0.79
## sample=UKB_25K p=0.01  rhog=0.8 i=1  EUR_EAS:0.77  EUR_AFR:0.80  EUR_SAS:0.73  EUR_AMR:0.74
## sample=UKB_25K p=5e-04 rhog=0.8 i=1  EUR_EAS:0.79  EUR_AFR:0.96  EUR_SAS:0.77  EUR_AMR:0.81
## sample=UKB_25K p=0.1   rhog=0.8 i=1  EUR_EAS:0.77  EUR_AFR:0.80  EUR_SAS:0.81  EUR_AMR:0.81
## sample=UKB_25K p=0.001 rhog=0.8 i=2  EUR_EAS:0.71  EUR_AFR:0.65  EUR_SAS:0.75  EUR_AMR:0.68
## sample=UKB_25K p=0.01  rhog=0.8 i=2  EUR_EAS:0.91  EUR_AFR:0.77  EUR_SAS:0.80  EUR_AMR:0.88
## sample=UKB_25K p=5e-04 rhog=0.8 i=2  EUR_EAS:0.78  EUR_AFR:0.79  EUR_SAS:0.91  EUR_AMR:0.83
## sample=UKB_25K p=0.1   rhog=0.8 i=2  EUR_EAS:0.85  EUR_AFR:0.81  EUR_SAS:0.79  EUR_AMR:0.77
## sample=UKB_25K p=0.001 rhog=0.8 i=3  EUR_EAS:0.77  EUR_AFR:0.69  EUR_SAS:0.81  EUR_AMR:0.85
## sample=UKB_25K p=0.01  rhog=0.8 i=3  EUR_EAS:0.83  EUR_AFR:0.80  EUR_SAS:0.70  EUR_AMR:0.70
## sample=UKB_25K p=5e-04 rhog=0.8 i=3  EUR_EAS:0.84  EUR_AFR:0.67  EUR_SAS:0.74  EUR_AMR:0.76
## sample=UKB_25K p=0.1   rhog=0.8 i=3  EUR_EAS:0.86  EUR_AFR:0.74  EUR_SAS:0.73  EUR_AMR:0.79
## sample=UKB_25K p=0.001 rhog=0.8 i=4  EUR_EAS:0.85  EUR_AFR:0.77  EUR_SAS:0.76  EUR_AMR:0.80
## sample=UKB_25K p=0.01  rhog=0.8 i=4  EUR_EAS:0.79  EUR_AFR:0.77  EUR_SAS:0.78  EUR_AMR:0.77
## sample=UKB_25K p=5e-04 rhog=0.8 i=4  EUR_EAS:0.65  EUR_AFR:0.71  EUR_SAS:0.73  EUR_AMR:0.86
## sample=UKB_25K p=0.1   rhog=0.8 i=4  EUR_EAS:0.78  EUR_AFR:0.82  EUR_SAS:0.74  EUR_AMR:0.82
## sample=UKB_25K p=0.001 rhog=0.8 i=5  EUR_EAS:0.78  EUR_AFR:0.75  EUR_SAS:0.72  EUR_AMR:0.77
## sample=UKB_25K p=0.01  rhog=0.8 i=5  EUR_EAS:0.82  EUR_AFR:0.71  EUR_SAS:0.80  EUR_AMR:0.82
## sample=UKB_25K p=5e-04 rhog=0.8 i=5  EUR_EAS:0.87  EUR_AFR:0.92  EUR_SAS:0.83  EUR_AMR:0.83
## sample=UKB_25K p=0.1   rhog=0.8 i=5  EUR_EAS:0.81  EUR_AFR:0.76  EUR_SAS:0.72  EUR_AMR:0.80
## sample=UKB_90K p=0.001 rhog=0.8 i=1  EUR_EAS:0.79  EUR_AFR:0.83  EUR_SAS:0.69  EUR_AMR:0.79
## sample=UKB_90K p=0.01  rhog=0.8 i=1  EUR_EAS:0.79  EUR_AFR:0.78  EUR_SAS:0.74  EUR_AMR:0.77
## sample=UKB_90K p=5e-04 rhog=0.8 i=1  EUR_EAS:0.82  EUR_AFR:0.87  EUR_SAS:0.83  EUR_AMR:0.76
## sample=UKB_90K p=0.1   rhog=0.8 i=1  EUR_EAS:0.78  EUR_AFR:0.78  EUR_SAS:0.78  EUR_AMR:0.82
## sample=UKB_90K p=0.001 rhog=0.8 i=2  EUR_EAS:0.72  EUR_AFR:0.66  EUR_SAS:0.73  EUR_AMR:0.69
## sample=UKB_90K p=0.01  rhog=0.8 i=2  EUR_EAS:0.85  EUR_AFR:0.82  EUR_SAS:0.78  EUR_AMR:0.84
## sample=UKB_90K p=5e-04 rhog=0.8 i=2  EUR_EAS:0.81  EUR_AFR:0.79  EUR_SAS:0.80  EUR_AMR:0.84
## sample=UKB_90K p=0.1   rhog=0.8 i=2  EUR_EAS:0.84  EUR_AFR:0.76  EUR_SAS:0.75  EUR_AMR:0.80
## sample=UKB_90K p=0.001 rhog=0.8 i=3  EUR_EAS:0.83  EUR_AFR:0.78  EUR_SAS:0.78  EUR_AMR:0.82
## sample=UKB_90K p=0.01  rhog=0.8 i=3  EUR_EAS:0.84  EUR_AFR:0.83  EUR_SAS:0.75  EUR_AMR:0.75
## sample=UKB_90K p=5e-04 rhog=0.8 i=3  EUR_EAS:0.84  EUR_AFR:0.73  EUR_SAS:0.73  EUR_AMR:0.70
## sample=UKB_90K p=0.1   rhog=0.8 i=3  EUR_EAS:0.81  EUR_AFR:0.75  EUR_SAS:0.75  EUR_AMR:0.79
## sample=UKB_90K p=0.001 rhog=0.8 i=4  EUR_EAS:0.83  EUR_AFR:0.71  EUR_SAS:0.75  EUR_AMR:0.82
## sample=UKB_90K p=0.01  rhog=0.8 i=4  EUR_EAS:0.84  EUR_AFR:0.79  EUR_SAS:0.80  EUR_AMR:0.74
## sample=UKB_90K p=5e-04 rhog=0.8 i=4  EUR_EAS:0.70  EUR_AFR:0.73  EUR_SAS:0.70  EUR_AMR:0.82
## sample=UKB_90K p=0.1   rhog=0.8 i=4  EUR_EAS:0.79  EUR_AFR:0.78  EUR_SAS:0.78  EUR_AMR:0.85
## sample=UKB_90K p=0.001 rhog=0.8 i=5  EUR_EAS:0.80  EUR_AFR:0.72  EUR_SAS:0.78  EUR_AMR:0.77
## sample=UKB_90K p=0.01  rhog=0.8 i=5  EUR_EAS:0.75  EUR_AFR:0.72  EUR_SAS:0.78  EUR_AMR:0.83
## sample=UKB_90K p=5e-04 rhog=0.8 i=5  EUR_EAS:0.83  EUR_AFR:0.88  EUR_SAS:0.83  EUR_AMR:0.86
## sample=UKB_90K p=0.1   rhog=0.8 i=5  EUR_EAS:0.79  EUR_AFR:0.78  EUR_SAS:0.74  EUR_AMR:0.79