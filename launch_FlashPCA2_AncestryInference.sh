#!/bin/bash
#SBATCH --account=ctb-hussinju
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --job-name=FlashPCA2
#SBATCH --output=launch_FlashPCA2_AncestryInference.out


module load nextflow


nextflow run flashPCA.nf -c nextflow.config -resume


echo "DONE"
