#! /bin/bash -login

## Exit if any command fails
set -e

## Load required modules
module load python/3.6.6

## Activate virtual environment
source env/bin/activate

case $1 in

    '-h' | '--help' | ' ' | '')
            echo -e '\e[31mSpecificy which workflow to unlock (i.e. genoProc or postImputation)'
            exit 2
            ;;
    
    'genoProc' | 'runGenopipe')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/genoProc --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status ./scripts/status.py
            ;;

    'postImputation' | 'runGenopipeImputation')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/postImputation --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status ./scripts/status.py
            ;;
esac

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun with sbatch runGenopipe sbatch runGenopipeImputation"