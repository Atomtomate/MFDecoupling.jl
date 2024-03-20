#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --ntasks 96
#SBATCH -p large96
##SBATCH -A hhp00048

#SBATCH -J U0.1V0.25

export SLURM_CPU_BIND=none





declare -i  i=0


U=1.0
V=0.5
g=-0.161844
f=0.197858

echo "$U  $V   $g  $f"
time ./main "$U" "$V" "$g" "$f" 0 300 

wait
