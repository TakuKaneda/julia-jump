#!/bin/bash
#SBATCH --job-name=cempc_test
#SBATCH --output=cempc_test_res.txt
#SBATCH --ntasks=4
#SBATCH --time=20:00
#SBATCH --mem-per-cpu=2000

julia parallel-run-cempc.jl
