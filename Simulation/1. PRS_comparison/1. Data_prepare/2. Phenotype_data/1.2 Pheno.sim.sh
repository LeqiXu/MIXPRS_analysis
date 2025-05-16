## 1. Phenotype Simulation
# simulate phenotype
# only discover participants
for pop in EUR EAS AFR SAS AMR
do
  mkdir -p /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/discover
  mkdir -p /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/validate
  mkdir -p /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/test
done

# discover participants
## non-EUR: 15K
rho=0.8 

for i in {1..5}; do 
for p in 0.001 0.01 5e-04 0.1; do
for pop in EAS AFR SAS AMR; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/discover/${pop}_sim${i}_p${p}_rho${rho}_15K.phen" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=phenotype_GWAS_simulation_${pop}_p${p}_rho${rho}_i${i}_d
#SBATCH --output=out_phenotype_GWAS_simulation_${pop}_p${p}_rho${rho}_i${i}_d.txt

  # discover phenotype
  if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/discover/${pop}_sim${i}_p${p}_rho${rho}_15K.phen" ]]; then
    /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/gcta64/gcta-1.94.1 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop}/discover/${pop} \
    --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/effect_data/${pop}_sim${i}_p${p}_rho${rho}.snplist \
    --keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/info/discover_15K_id.tsv \
    --simu-qt \
    --simu-causal-loci /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/effect_data/${pop}_sim${i}_p${p}_rho${rho}.txt \
    --simu-hsq 0.4 \
    --simu-rep 1 \
    --out /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/discover/${pop}_sim${i}_p${p}_rho${rho}_15K
  fi

EOT
fi
done
done
done

## non-EUR: 80K
rho=0.8 

for i in {1..5}; do 
for p in 0.001 0.01 5e-04 0.1; do
for pop in EAS AFR SAS AMR; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/discover/${pop}_sim${i}_p${p}_rho${rho}_80K.phen" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=phenotype_GWAS_simulation_${pop}_p${p}_rho${rho}_i${i}_l
#SBATCH --output=out_phenotype_GWAS_simulation_${pop}_p${p}_rho${rho}_i${i}_l.txt

  if [[ ${pop} != "EUR" ]] && [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/discover/${pop}_sim${i}_p${p}_rho${rho}_80K.phen" ]]; then
  /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/gcta64/gcta-1.94.1 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop}/discover/${pop} \
  --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/effect_data/${pop}_sim${i}_p${p}_rho${rho}.snplist \
  --keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/info/discover_80K_id.tsv \
  --simu-qt \
  --simu-causal-loci /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/effect_data/${pop}_sim${i}_p${p}_rho${rho}.txt \
  --simu-hsq 0.4 \
  --simu-rep 1 \
  --out /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/discover/${pop}_sim${i}_p${p}_rho${rho}_80K
  fi
EOT
fi
done
done
done

## EUR: UKB
rho=0.8 
pop=EUR

for i in {1..5}; do 
for p in 0.001 0.01 5e-04 0.1; do
if [[ ${pop} == "EUR" ]] && [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/discover/${pop}_sim${i}_p${p}_rho${rho}_UKB.phen" ]]; then
sbatch <<EOT
#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --requeue
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=phenotype_GWAS_simulation_${pop}_p${p}_rho${rho}_i${i}_l
#SBATCH --output=out_phenotype_GWAS_simulation_${pop}_p${p}_rho${rho}_i${i}_l.txt

    /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/gcta64/gcta-1.94.1 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/data/ukbb_data/geno_data/${pop} \
    --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/effect_data/${pop}_sim${i}_p${p}_rho${rho}.snplist \
    --simu-qt \
    --simu-causal-loci /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/effect_data/${pop}_sim${i}_p${p}_rho${rho}.txt \
    --simu-hsq 0.4 \
    --simu-rep 1 \
    --out /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/discover/${pop}_sim${i}_p${p}_rho${rho}_UKB
EOT
fi
done
done


# validate participants and test participants
for pop in EUR EAS AFR SAS AMR
do
  mkdir -p /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/validate
  mkdir -p /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/test
done

rho=0.8 

for i in {1..5}; do 
for p in 0.001 0.01 5e-04 0.1; do
for pop in EUR EAS AFR SAS AMR; do
if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/validate/${pop}_sim${i}_p${p}_rho${rho}_10K.phen" ]] || [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/test/${pop}_sim${i}_p${p}_rho${rho}_10K.phen" ]]; then
        sbatch <<EOT
#!/bin/bash
#SBATCH --partition=scavenge,day,week
#SBATCH --requeue
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=phenotype_GWAS_simulation_${pop}_p${p}_rho${rho}_i${i}
#SBATCH --output=out_phenotype_GWAS_simulation_${pop}_p${p}_rho${rho}_i${i}.txt

  # validate phenotype
  if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/validate/${pop}_sim${i}_p${p}_rho${rho}_10K.phen" ]]; then
  /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/gcta64/gcta-1.94.1 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop}/validate/${pop} \
  --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/effect_data/${pop}_sim${i}_p${p}_rho${rho}.snplist \
  --keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/info/validate_id.tsv \
  --simu-qt \
  --simu-causal-loci /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/effect_data/${pop}_sim${i}_p${p}_rho${rho}.txt \
  --simu-hsq 0.4 \
  --simu-rep 1 \
  --out /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/validate/${pop}_sim${i}_p${p}_rho${rho}_10K
  fi

  # test phenotype
  if [[ ! -e "/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/test/${pop}_sim${i}_p${p}_rho${rho}_10K.phen" ]]; then
  /gpfs/gibbs/pi/zhao/lx94/JointPRS/method/gcta64/gcta-1.94.1 --bfile /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/${pop}/test/${pop} \
  --extract /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/effect_data/${pop}_sim${i}_p${p}_rho${rho}.snplist \
  --keep /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/geno_data/info/test_id.tsv \
  --simu-qt \
  --simu-causal-loci /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/effect_data/${pop}_sim${i}_p${p}_rho${rho}.txt \
  --simu-hsq 0.4 \
  --simu-rep 1 \
  --out /gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/${pop}/test/${pop}_sim${i}_p${p}_rho${rho}_10K
  fi
EOT
fi
done
done
done


## 2. Organize pheno data for validation and testing
# validate and test participants pheno data from simulation genotype data
library(data.table)

rho=0.8 

for (i in c(1:5)){ 
for (p in c(0.001,0.01,0.0005,0.1)){
for (pop in c("EUR","EAS","AFR","SAS","AMR")){
  val_data_all = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/",pop,"/validate/",pop,"_sim",i,"_p",p,"_rho",rho,"_10K.phen"))
  test_data_all = fread(paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/",pop,"/test/",pop,"_sim",i,"_p",p,"_rho",rho,"_10K.phen"))
        
  val_data_10K = val_data_all[1:10000,]
  test_data_10K = test_data_all[1:10000,]

  val_data_10K[,c(3)] = scale(val_data_10K[,c(3)])
  test_data_10K[,c(3)] = scale(test_data_10K[,c(3)])

  colnames(val_data_10K) = c("FID","IID","pheno")
  colnames(test_data_10K) = c("FID","IID","pheno")

  write.table(val_data_10K, paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/",pop,"/validate/",pop,"_sim",i,"_p",p,"_rho",rho,"_10K_doubleid.tsv"), row.names=F, col.names=F, quote = F, sep = "\t")
  write.table(test_data_10K, paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/",pop,"/test/",pop,"_sim",i,"_p",p,"_rho",rho,"_10K_doubleid.tsv"), row.names=F, col.names=F, quote = F, sep = "\t")

  write.table(val_data_10K, paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/",pop,"/validate/",pop,"_sim",i,"_p",p,"_rho",rho,"_10K_doubleidname.tsv"), row.names=F, col.names=T, quote = F, sep = "\t")
  write.table(test_data_10K, paste0("/gpfs/gibbs/pi/zhao/lx94/JointPRS/revision/data/sim_data/pheno_data/",pop,"/test/",pop,"_sim",i,"_p",p,"_rho",rho,"_10K_doubleidname.tsv"), row.names=F, col.names=T, quote = F, sep = "\t")  
}
}
}
