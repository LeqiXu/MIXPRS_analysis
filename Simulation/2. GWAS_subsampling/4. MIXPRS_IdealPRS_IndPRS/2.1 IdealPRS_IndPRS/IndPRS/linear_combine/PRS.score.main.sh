# Generate params list
# SDPRX_params="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/SDPRX/params.txt"
# SDPRX_out="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/SDPRX.params.txt"
# JointPRS_params="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/JointPRS/params.txt"
# JointPRS_out="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/JointPRS.params.txt"

# awk '{
#     NF-=2; line=$0; pop1=$1; pop2=$3
#     print line, pop1;
#     print line, pop2;
# }' $SDPRX_params | sort -u > $SDPRX_out

# awk 'BEGIN{
#     split("EUR AFR EAS AMR SAS", pops, " ")
# }
# {
#     NF-=2; line=$0
#     for(i in pops) print line, pops[i]
# }' $JointPRS_params | sort -u > $JointPRS_out


# PRS score
SDPRX_out="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/SDPRX.params.txt"
JointPRS_out="/gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/JointPRS.params.txt"

SDPRX_job=$(wc -l < ${SDPRX_out})
JointPRS_job=$(wc -l < ${JointPRS_out})

sbatch --array=1-${SDPRX_job} /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/PRS.SDPRX.inner.sh
sbatch --array=1-${JointPRS_job} /gpfs/gibbs/pi/zhao/xz674/simulation/b.method_calculate/IndPRS/linear_combine/PRS.JointPRS.inner.sh

# 1-${SDPRX_job}
# 1-${JointPRS_job}
