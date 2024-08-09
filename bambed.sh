#!/bin/bash

#SBATCH --job-name=extract_regions
#SBATCH --output=extract_regions_%j.out
#SBATCH --error=extract_regions_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Load your conda environment
source /users/hinch/vxy098/miniconda3/etc/profile.d/conda.sh
conda activate denovo

# Define paths to input and intermediate output
original_bam="/well/hinch/projects/Brca2/PacBio/Trinity/results/aligned/pup1.sorted.bam"
input_bed="/well/hinch/projects/Brca2/PacBio/Trinity/results/sv/pup1/delta_0/haplotype.bed"
intermediate_bam="/well/hinch/projects/Brca2/PacBio/Trinity/results/sv/pup1/delta_0/pup1_intermidiate.bam"


# Filter BAM based on BED file
samtools view -b -L "$input_bed" -o "$intermediate_bam" "$original_bam"

echo "Intermediate BAM file created at $intermediate_bam"