## Step0: Obtain PRScsx format GWAS
## subsample PRS
library(data.table)

h2 = 0.4
rhog = 0.8

type = "prune_snplist_1"
approx_list = c("TRUE")

for (sample_size in c("25K","90K")){
for (pop in c("EAS","AFR","SAS","AMR")){
for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

for (rpt in c(1:4)){

for (approx in approx_list){

sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/clean/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_", type, "_", pop,"_train_GWAS_approx",approx,"_ratio3.00_repeat",rpt,".txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","BETA","P")]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/subsample_data/PRScsx/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_", type, "_train_GWAS_approx",approx,"_PRScsx_repeat",rpt,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

}

}
}

}
}
}

## no subsample prune PRS
library(data.table)

h2 = 0.4
rhog = 0.8
pop = "EUR"
sample_size = "UKB"

for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/",pop,"_prune_pval1_r20.5_wc250_1.snplist"), header = FALSE)
sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/summary_data/",pop,"/discover/clean/",pop,"_sim",sim_i,"_p",p,"_rho",rhog,"_",sample_size,"_clean_real.txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","BETA","P")]
sumstat_data = sumstat_data[which(sumstat_data$SNP %in% prune_snplist$V1),]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/PRScsx/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_PRScsx_real_prune.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}

library(data.table)

h2 = 0.4
rhog = 0.8

for (sample_size in c("25K","90K")){
for (pop in c("EAS","AFR","SAS","AMR")){
for (sim_i in c(1:5)){
for (p in c(0.1, 0.01, 0.001, 5e-04)){

prune_snplist = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/prune_clump/snplist/",pop,"_prune_pval1_r20.5_wc250_1.snplist"), header = FALSE)
sumstat_data = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/summary_data/",pop,"/discover_validate/clean/",pop,"_sim",sim_i,"_p",p,"_rho",rhog,"_",sample_size,"_clean_real.txt"))
sumstat_data = sumstat_data[,c("SNP","A1","A2","BETA","P")]
sumstat_data = sumstat_data[which(sumstat_data$SNP %in% prune_snplist$V1),]

write.table(sumstat_data,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/data/sim_data/summary_data/discover_validate/PRScsx/", pop,"_sim",sim_i, "_h2", h2, "_p", p, "_rhog", rhog, "_", sample_size, "_PRScsx_real_prune.txt"),quote=F,sep='\t',row.names=F,col.names=T)

}
}
}
}

# Step1: Estimate beta
for subpop in EAS AFR SAS AMR; do

job_file="/gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/JointPRS_scenario2_${subpop}.txt"
> $job_file  # Empty the job file if it already exists

h2=0.4
rhog=0.8

sample1=UKB

type=prune_snplist_1
type_original="prune"
approx_list="TRUE"

for sample2 in 25K 90K; do
for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for rpt in {1..4}; do

for approx in ${approx_list}; do

# sample size
if [[ "${subpop}" == "EAS" && "${sample2}" == "25K" ]]; then
  sample_size1=311600; sample_size2=18750;  sample_size3=25000; sample_size4=25000; sample_size5=25000; out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_subEAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "EAS" && "${sample2}" == "90K" ]]; then
  sample_size1=311600; sample_size2=67500;  sample_size3=90000; sample_size4=90000; sample_size5=90000; out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_subEAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "AFR" && "${sample2}" == "25K" ]]; then
  sample_size1=311600; sample_size2=25000; sample_size3=18750;  sample_size4=25000; sample_size5=25000; out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_subAFR_SAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "AFR" && "${sample2}" == "90K" ]]; then
  sample_size1=311600; sample_size2=90000; sample_size3=67500;  sample_size4=90000; sample_size5=90000; out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_subAFR_SAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "SAS" && "${sample2}" == "25K" ]]; then
  sample_size1=311600; sample_size2=25000; sample_size3=25000; sample_size4=18750;  sample_size5=25000; out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_subSAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "SAS" && "${sample2}" == "90K" ]]; then
  sample_size1=311600; sample_size2=90000; sample_size3=90000; sample_size4=67500;  sample_size5=90000; out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_subSAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "AMR" && "${sample2}" == "25K" ]]; then
  sample_size1=311600; sample_size2=25000; sample_size3=25000; sample_size4=25000; sample_size5=18750;  out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_subAMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "AMR" && "${sample2}" == "90K" ]]; then
  sample_size1=311600; sample_size2=90000; sample_size3=90000; sample_size4=90000; sample_size5=67500;  out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_subAMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
else
  echo "Please provide a valid subpop: EUR, EAS, AFR, SAS, or AMR, and sample2 as 25K or 90K."
fi

for pop in EUR EAS AFR SAS AMR; do
  if [[ "$pop" == "$subpop" ]]; then
    if [[ "$pop" == "EUR" ]]; then
        file="data/sim_data/summary_data/subsample_data/PRScsx/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_${type}_train_GWAS_approx${approx}_PRScsx_repeat${rpt}.txt"
    else
        file="data/sim_data/summary_data/subsample_data/PRScsx/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_${type}_train_GWAS_approx${approx}_PRScsx_repeat${rpt}.txt"
    fi
  else
    if [[ "$pop" == "EUR" ]]; then
        file="data/sim_data/summary_data/discover_validate/PRScsx/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample1}_PRScsx_real_${type_original}.txt"
    else
        file="data/sim_data/summary_data/discover_validate/PRScsx/${pop}_sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_${sample2}_PRScsx_real_${type_original}.txt"
    fi
  fi
  declare "${pop}_sumstat=${file}"
done

sst_file="${EUR_sumstat},${EAS_sumstat},${AFR_sumstat},${SAS_sumstat},${AMR_sumstat}"

for chr in {1..22}; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/subsample/JointPRS/${out_name}_AMR_pst_eff_a1_b0.5_phiauto_chr${chr}.txt" ]]; then

echo "module load miniconda; conda activate py_env; cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/; python /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/JointPRS/JointPRS.py --ref_dir=/gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ref_data/PRScsx/1kg --bim_prefix=/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/All/All_test --sst_file=${sst_file} --rho_cons=1,1,1,1,1 --n_gwas=${sample_size1},${sample_size2},${sample_size3},${sample_size4},${sample_size5} --chrom=${chr} --pop=EUR,EAS,AFR,SAS,AMR --out_dir=result/sim_result/subsample/JointPRS --out_name=${out_name}" >> $job_file

fi

done
done

done
done
done
done

module load dSQ
cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/

dsq --job-file /gpfs/gibbs/pi/zhao/lx94/SWIFT/code/simulation/PRS/JointPRS_scenario2_${subpop}.txt --partition=scavenge,day --requeue --mem=10G --cpus-per-task=1 --ntasks=1 --nodes=1 --time=24:00:00 --mail-type=ALL
sbatch dsq-JointPRS_scenario2_${subpop}-$(date +%Y-%m-%d).sh

done

# Step2: Organize beta by chr pop for each param in each scenario
library(data.table)

h2=0.4
rhog=0.8 

sample1="UKB"

type = "prune_snplist_1"
approx_list = c("TRUE")

for (sample2 in c("25K","90K")){
for (sim_i in c(1:5)){
for (p in c(0.001,0.01,5e-04,0.1)){
for (pop in c("EUR","EAS","AFR","SAS","AMR")){

for (rpt in c(1:4)){

for (approx in approx_list){

for (subpop in c("EAS","AFR","SAS","AMR")){

if (subpop == "EAS") {
    out_name = paste0("sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_subEAS_AFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_subsample_",type,"_approx",approx,"_repeat",rpt)
} else if (subpop == "AFR") {
    out_name = paste0("sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_subAFR_SAS_AMR_",sample1,"_",sample2,"_JointPRS_subsample_",type,"_approx",approx,"_repeat",rpt)
} else if (subpop == "SAS") {
    out_name = paste0("sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_subSAS_AMR_",sample1,"_",sample2,"_JointPRS_subsample_",type,"_approx",approx,"_repeat",rpt)
} else if (subpop == "AMR") {
    out_name = paste0("sim",sim_i,"_h2",h2,"_p",p,"_rhog",rhog,"_EUR_EAS_AFR_SAS_subAMR_",sample1,"_",sample2,"_JointPRS_subsample_",type,"_approx",approx,"_repeat",rpt)
} else {
    next
}

JointPRS_all <- data.table()
for(chr in 1:22){
    
    JointPRS_pop_chr <- fread(paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/subsample/JointPRS/",out_name,"_",pop,"_pst_eff_a1_b0.5_phiauto_chr",chr,".txt"))
    JointPRS_pop_chr <- JointPRS_pop_chr[,c(2,4,6)]
    names(JointPRS_pop_chr) = c("rsID","A1",pop)

    JointPRS_all = rbind(JointPRS_all,JointPRS_pop_chr)
    
}

write.table(JointPRS_all,paste0("/gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/subsample/JointPRS/",out_name,"_beta_",pop,".txt"),quote=F,sep='\t',row.names=F,col.names=T)

}

}
}
}

}
}
}

# Step3: Clean the previous result
h2=0.4
rhog=0.8

sample1=UKB

type=prune_snplist_1
approx_list="TRUE"

for sample2 in 25K 90K; do
for subpop in EUR EAS AFR SAS AMR; do
for sim_i in {1..5}; do
for p in 0.001 0.01 5e-04 0.1; do
for rpt in {1..4}; do

for approx in ${approx_list}; do

# sample size
if [[ "${subpop}" == "EUR" ]]; then
  out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_subEUR_EAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "EAS" ]]; then
  out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_subEAS_AFR_SAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "AFR" ]]; then
  out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_subAFR_SAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "SAS" ]]; then
  out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_subSAS_AMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
elif [[ "${subpop}" == "AMR" ]]; then
  out_name="sim${sim_i}_h2${h2}_p${p}_rhog${rhog}_EUR_EAS_AFR_SAS_subAMR_${sample1}_${sample2}_JointPRS_subsample_${type}_approx${approx}_repeat${rpt}"
else
  echo "Please provide a valid subpop: EUR, EAS, AFR, SAS, or AMR."
fi

for chr in {1..22}; do

rm -rf /gpfs/gibbs/pi/zhao/lx94/SWIFT/result/sim_result/subsample/JointPRS/${out_name}_AMR_pst_eff_a1_b0.5_phiauto_chr${chr}.txt

done
done

done
done
done
done
done