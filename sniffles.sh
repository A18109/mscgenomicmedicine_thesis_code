#!/bin/bash
#SBATCH --job-name=sniffles_variants
#SBATCH --output=sniffles_variant_%j.out
#SBATCH --error=sniffles_variant_%j.err
#SBATCH --time=30:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
# Load conda environment
source /users/hinch/vxy098/miniconda3/etc/profile.d/conda.sh
conda activate denovo
##
INPUT_BAM="/well/hinch/projects/Brca2/PacBio/Trinity/results/aligned/sample.sorted.bam"
OUTPUT_VCF="/well/hinch/projects/Brca2/PacBio/Trinity/results/variants/sniffles/sample_sniffles.vcf"
###
sniffles --input $INPUT_BAM  --vcf $OUTPUT_VCF
echo "Sniffles SV calling for father complete."
