#!/bin/bash
#
#SBATCH --job-name=single_test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=8096
#SBATCH --partition=hmem

module load OpenMPI
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

while getopts U:p:g: flag
do
    case "${flag}" in
        U) U_s=${OPTARG};;
        p) phi=${OPTARG};;
        g) gamma=${OPTARG};;
    esac
done

./solver $U_s $phi $gamma
