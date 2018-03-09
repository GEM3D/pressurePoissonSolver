#!/bin/bash
filebase="r2_strong"
cmdbase="./steady --mesh ./mesh_4levels.txt -l 10 -n 32 --gauss --neumann --divide 2"
petscopts="
-ksp_type gmres -ksp_pc_side right -pc_type hypre -pc_hypre_boomeramg_relax_type_all Jacobi
-pc_hypre_boomeramg_relax_type_coarse Jacobi -ksp_monitor_true_residual ascii
-pc_hypre_boomeramg_print_statistics 1 -ksp_converged_reason ascii "

PROCS=1
sbatchbase=" --ntasks-per-node 24 "
while [ $PROCS -le 64 ]; do
    sbatch $sbatchbase --job-name="${filebase}_${PROCS}" -n $PROCS \
        -e ${filebase}_${PROCS}.err -o ${filebase}_${PROCS}.out \
        slurmrun.sh $cmdbase :: $petscopts
    let PROCS*=2
done
