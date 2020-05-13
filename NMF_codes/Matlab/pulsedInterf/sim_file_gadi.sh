#!/bin/bash
#PBS -q normal
#PBS -l walltime=3:00:00
#PBS -l ncpus=16
#PBS -l mem=32GB
#PBS -l jobfs=1GB
#PBS -l software=matlab_mq
#PBS -l wd
 
module load matlab/R2019b
nohup matlab -nodisplay -nosplash -nojvm <test_chirp_mitigation_nmf_semi.m>& test_chirp_mitigation_nmf_semi&
