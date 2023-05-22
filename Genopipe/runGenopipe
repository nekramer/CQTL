#!/bin/bash
#SBATCH -J Genopipe
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p general
#SBATCH -o "%x-%j.out"

## Load python
module load python/3.6.6

## Load plink
module load plink

## Create and activate virtual environment with requirements
python3 -m venv env && source env/bin/activate && pip3 install -r config/requirements.txt

## Make directory for slurm logs
mkdir -p output/logs_slurm

## Execute workflow
snakemake -j 1 --latency-wait 500 -s snakefiles/genoProc --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --constraint={cluster.constraint} --parsable" --cluster-status ./scripts/status.py