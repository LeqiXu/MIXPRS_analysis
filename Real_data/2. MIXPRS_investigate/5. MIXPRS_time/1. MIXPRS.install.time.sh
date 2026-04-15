sbatch <<EOT
#!/bin/bash
#SBATCH --partition=day
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=mixprs_conda_test
#SBATCH --output=out_mixprs_conda_test.txt

module load miniconda

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/github/MIXPRS
/usr/bin/time -v conda env create -f environment.yml -n mixprs_conda_test
EOT

sbatch <<EOT
#!/bin/bash
#SBATCH --partition=day
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name=mixprs_pip_test
#SBATCH --output=out_mixprs_pip_test.txt

module load miniconda

cd /gpfs/gibbs/pi/zhao/lx94/SWIFT/github/MIXPRS
/usr/bin/time -v bash -lc '
  source $(conda info --base)/etc/profile.d/conda.sh
  conda create -y -n mixprs_pip_test python=3.8
  conda activate mixprs_pip_test
  pip install -r requirements.txt
'
EOT
