#!/bin/bash
#SBATCH --job-name=freebayes
#SBATCH --output=freebayes_%j.out
#SBATCH --error=freebayes_%j.err
#SBATCH --time=30:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

# Activate environment
source /users/hinch/vxy098/miniconda3/etc/profile.d/conda.sh
conda activate denovo

# Define input files and parameters
reference_genome="/well/hinch/projects/Brca2/PacBio/Trinity/reference/mm10_gatkorder.fa"
input_vcf="/well/hinch/projects/Brca2/PacBio/Trinity/results/sv/pup1/delta_0/mgpreduced.vcf.gz"
sorted_bam="/well/hinch/projects/Brca2/PacBio/Trinity/results/sv/pup1/delta_0/pup1_denovo_sorted.bam"
output_vcf="/well/hinch/projects/Brca2/PacBio/Trinity/results/sv/pup1/delta_0/freebayes/freebayes_haplotype.vcf"

# Ensure all necessary files are present
if [ ! -f "$reference_genome" ] || [ ! -f "$input_vcf" ] || [ ! -f "$sorted_bam" ]; then
    echo "Error: One or more input files do not exist."
    exit 1
fi

# Run FreeBayes for variant calling
echo "Running FreeBayes for variant calling..."
freebayes -f "$reference_genome" -@ "$input_vcf" "$sorted_bam" > "$output_vcf" 2>&1

# Check if FreeBayes was successful
if [ $? -ne 0 ]; then
    echo "Error: FreeBayes variant calling failed."
    echo "Review the 'freebayes_$SLURM_JOBID.err' file for details."
    exit 1
fi

echo "Variant calling completed successfully."
