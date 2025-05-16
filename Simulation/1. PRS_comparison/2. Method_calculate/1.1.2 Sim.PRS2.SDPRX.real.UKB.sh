## Step0: Obtain SDPRX format GWAS
library(data.table)

h2 = 0.4
rhog = 0.8
pop = "EUR"
sample_size = "UKB"

for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/summary_data/",pop,"/discover/clean/",pop,"_sim",sim_i,"_p",p,"_rho",rhog,"_",sample_size,"_clean_real.txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","Z","P","N")]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/SDPRX/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_SDPRX_real_all.txt"),quote=F,sep='\t',row.names=F,col.names=T)

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
sumstat_data = sumstat_data[,c("SNP","A1","A2","Z","P","N")]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/SDPRX/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_SDPRX_real_all.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}
}
}

# Step1: Estimate beta
job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/SDPRX_scenario1.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

pop1=EUR
sample1=UKB

for sample2 in 25K 90K; do
for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for pop2 in EAS AFR SAS AMR; do

# sample size
if [[ ${sample2} == "25K" ]]; then
sample_size1=311600; sample_size2=25000
elif [[ ${sample2} == "90K" ]]; then
sample_size1=311600; sample_size2=90000
else
echo "Please provide the available phenotype"
fi

file=/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/popcorn/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_${pop2}_${sample1}_${sample2}_popcorn_real_corr.txt
rho_est=$(grep '^pge' "${file}" | awk '{printf "%.2f", $2}')

for chr in {1..22}; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_${pop2}_${sample1}_${sample2}_SDPRX_real_chr${chr}_2.txt" ]]; then

echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/SDPRX/SDPRX.py --load_ld /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/SDPRX/EUR_${pop2} --valid /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/All/All_test.bim --ss1 data/sim_data/summary_data/discover_validate/SDPRX/${pop1}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_SDPRX_real_all.txt --ss2 data/sim_data/summary_data/discover_validate/SDPRX/${pop2}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_SDPRX_real_all.txt --N1 ${sample_size1} --N2 ${sample_size2} --mcmc_samples 2000 --burn 1000 --force_shared True --chr ${chr} --rho ${rho_est} --out result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_${pop2}_${sample1}_${sample2}_SDPRX_real_chr${chr}" >> $job_file

fi
done
done
done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/SDPRX_scenario1.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-SDPRX_scenario1-$(date +%Y-%m-%d).sh


# Step2: Organize beta by chr pop for each param in each scenario
library(data.table)

h2=0.4
rhog=0.8

sample1="UKB"

for (sample2 in c("25K","90K")){
for (sim_i in c(1:5)){
for (p in c(0.001,0.01,5e-04,0.1)){
for (pop in c("EAS","AFR","SAS","AMR")){

SDPRX_all <- data.table()
for(chr in 1:22){
    
    SDPRX_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_EUR_",pop,"_",sample1,"_",sample2,"_SDPRX_real_chr",chr,"_2.txt"))
    names(SDPRX_pop_chr) = c("rsID","A1",pop)

    SDPRX_all = rbind(SDPRX_all,SDPRX_pop_chr)
    
}

write.table(SDPRX_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog,"_",sample1,"_",sample2,"_SDPRX_real_EUR_",pop,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)


SDPRX_all <- data.table()
for(chr in 1:22){
    
    SDPRX_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_EUR_",pop,"_",sample1,"_",sample2,"_SDPRX_real_chr",chr,"_1.txt"))
    names(SDPRX_pop_chr) = c("rsID","A1",pop)

    SDPRX_all = rbind(SDPRX_all,SDPRX_pop_chr)
    
}

write.table(SDPRX_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog,"_",sample1,"_",sample2,"_SDPRX_real_EUR_",pop,"_beta_EUR.txt"),quote=F,sep='\t',row.names=F,col.names=T)


}
}
}
}

# Step3: Clean the previous result
h2=0.4
rhog=0.8

pop1=EUR
sample1=UKB

for sample2 in 25K 90K; do
for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for pop2 in EAS AFR SAS AMR; do
for chr in {1..22}; do

rm -rf /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/SDPRX/sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${pop1}_${pop2}_${sample1}_${sample2}_SDPRX_real_chr${chr}_2.txt

done
done
done
done
done