#!/bin/bash

#SBATCH --job-name=sort_index_bam
#SBATCH --output=sort_index_bam_%j.out
#SBATCH --error=sort_index_bam_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=short

# Load your conda environment
source /users/hinch/vxy098/miniconda3/etc/profile.d/conda.sh
conda activate denovo

# Define paths to intermediate and final output
intermediate_bam="/well/hinch/projects/Brca2/PacBio/Trinity/results/sv/pup1/delta_0/intermediate.bam"
rnames_txt="/well/hinch/projects/Brca2/PacBio/Trinity/results/sv/pup1/delta_0/pup1_corrected_rnames.txt"
sorted_output_bam="/well/hinch/projects/Brca2/PacBio/Trinity/results/sv/pup1/delta_0/pup1_denovo_sorted.bam"

# Filter intermediate BAM for specific read names
# Filter intermediate BAM for specific read names
samtools view -bh -N "$rnames_txt" "$intermediate_bam" > "$sorted_output_bam"

# Sort and index the filtered BAM
samtools sort "$sorted_output_bam" -o "$sorted_output_bam"
samtools index "$sorted_output_bam"



echo "Sorted and indexed BAM created at $sorted_output_bam"
