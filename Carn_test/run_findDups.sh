#!/bin/bash
#SBATCH -n 2
#SBATCH -p fast
#SBATCH -t 02:00:00
#SBATCH --mem=60GB
#SBATCH -e Carn_acyltransf.err
#SBATCH -o Carn_acyltransf.out



python3.6 ../findDups_countLoss.py ../Carn_acyltransf.faa.nw ../Carn_acyltransf.faa ./23072021_1/ midpoint
