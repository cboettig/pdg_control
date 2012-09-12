## Use the current working directory
#$ -cwd
## use bash commands
#$ -S /bin/bash
## combine error and output files
#$ -j y
## Parallel for MPI
##$ -pe mpi 20 
## Parallel for openmp:
#$ -pe threaded 4
## Output file name
#$ -o det_figure3.log

module load gcc openmpi R 
Rscript -e "source('~/.Rprofile'); library('knitr'); knit('figure3.Rmd')" 


