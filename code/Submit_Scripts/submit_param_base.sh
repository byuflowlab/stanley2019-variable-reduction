#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1G
#SBATCH -J "P base"

module load python/2.7

spacing=4.0
boundary='amalia'
rose='northIslandRose'

python optParamRevision.py $spacing $boundary $rose 1


