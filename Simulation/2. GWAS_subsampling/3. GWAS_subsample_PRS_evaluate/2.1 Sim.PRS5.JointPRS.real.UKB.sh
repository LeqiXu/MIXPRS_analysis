## Step0: Obtain PRScsx format GWAS
library(data.table)

h2 = 0.4
rhog = 0.8
pop = "EUR"
sample_size = "ukbb"

for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover_validate/clean/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_clean_real_all.txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","BETA","P")]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/PRScsx/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_PRScsx_real_all.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}

library(data.table)

h2 = 0.4
rhog = 0.8
sample_size = "100K"

for (pop in c("EAS","AFR","SAS","AMR")){
for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/xz674/data/sim_data/subsample_GWAS/summary_data/",pop,"/discover_validate/clean/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_clean_real_all.txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","BETA","P")]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/PRScsx/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_PRScsx_real_all.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}
}


# Step1: Estimate beta
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/JointPRS_scenario1.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

sample1=ukbb
sample2=100K

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do

# sample size
if [[ ${sample2} == "100K" ]]; then
sample_size1=311600; sample_size2=100000
else
echo "Please provide the available phenotype"
fi

for chr in {1..22}; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_real_AMR_pst_eff_a1_b0.5_phiauto_chr${chr}.txt" ]]; then

echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/All/All_test --sst_file=data/sim_data/summary_data/discover_validate/PRScsx/EUR_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_PRScsx_real_all.txt,data/sim_data/summary_data/discover_validate/PRScsx/EAS_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_PRScsx_real_all.txt,data/sim_data/summary_data/discover_validate/PRScsx/AFR_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_PRScsx_real_all.txt,data/sim_data/summary_data/discover_validate/PRScsx/SAS_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_PRScsx_real_all.txt,data/sim_data/summary_data/discover_validate/PRScsx/AMR_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_PRScsx_real_all.txt --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size2},${sample_size2},${sample_size2} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/sim_result/JointPRS --out_name=sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_real" >> $job_file


fi

done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/JointPRS_scenario1.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-JointPRS_scenario1-$(date +%Y-%m-%d).sh


# Step2: Organize beta by chr pop for each param in each scenario
library(data.table)

h2=0.4
rhog=0.8 

sample1="ukbb"
sample2="100K"

for (sim_i in c(1:5)){
for (p in c(0.001,0.01,5e-04,0.1)){
for (pop in c("EUR","EAS","AFR","SAS","AMR")){

JointPRS_all <- data.table()
for(chr in 1:22){
    
    JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_EUR_EAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_real_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("rsID","A1",pop)

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_",sample1,"_",sample2,"_JointPRS_real_EUR_EAS_AFR_SAS_AMR_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}
}


# Step3: Clean the previous result
h2=0.4
rhog=0.8

sample1=ukbb
sample2=100K

for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for chr in {1..22}; do
for pop in EUR EAS AFR SAS AMR; do

rm -rf /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/JointPRS/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_real_${pop}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt

done
done
done
done