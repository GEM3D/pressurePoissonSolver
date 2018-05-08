#!/bin/bash
filebase="amgxstrong"
cmdbase="./steady --mesh ./mesh_4levels.txt -l 10 -n 32 --gauss --amgx ./amgx.json --neumann --divide 2"
sbatchbase="--gres=gpu:2 --ntasks-per-node 2 --ntasks-per-socket 2 --partition tesla -m block "
petscopts=""

PROCS=1
while [ $PROCS -le 16 ]; do
    sbatch $sbatchbase --job-name="${filebase}_${PROCS}" -n $PROCS \
        -e ${filebase}_${PROCS}.err -o ${filebase}_${PROCS}.out \
        slurmrun.sh $cmdbase :: $petscopts
    let PROCS*=2
done
