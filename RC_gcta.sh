#!/bin/bash

#PBS -N z.path
#PBS -l walltime=00:00:60:00
#PBS -l nodes=2:ppn=10
#PBS -l naccesspolicy=singletask
#PBS -q janus-small
#PBS -d /projects/KellerLab/simonson/
#PBS -o /projects/KellerLab/simonson/
#PBS -e /projects/KellerLab/simonson/


. /curc/tools/utils/dkinit

gcta --reml --mgrm giant.mgrm --pheno BMI.res.phe --reml-no-lrt --reml-maxit 20 --out giant.genes --thread-num 10
