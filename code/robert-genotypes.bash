#!/bin/bash

# Set an environment variable for the output file path
export OUTPUT_FILE="/net/mulan/disk2/fredboe/research/ml/cluster_outputs/robert-genotypes.out"
export ERR_FILE="/net/mulan/disk2/fredboe/research/ml/cluster_outputs/robert-genotypes.err"

# Launch a Slurm job that uses the environment variable for the output file
sbatch <<EOT
#!/bin/bash
#SBATCH --job-name=quarto-robert-genotypes
#SBATCH --output=$OUTPUT_FILE
#SBATCH --err=$ERR_FILE
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1-00:00:00

quarto render robert-genotypes.qmd --to gfm

EOT