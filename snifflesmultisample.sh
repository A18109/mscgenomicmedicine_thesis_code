#!/bin/bash
#SBATCH --job-name=sniffles_variants
#SBATCH --output=sniffles_variant_%j.out
#SBATCH --error=sniffles_variant_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4


# Load conda environment
source /users/hinch/vxy098/miniconda3/etc/profile.d/conda.sh
# This command sources the script that configures the environment to use conda.

conda activate denovo
# This command activates a specific conda environment named 'denovo'.

# Define input and output directories
INPUT_DIR="/well/hinch/projects/Brca2/PacBio/Trinity/results/aligned"
# Path to the directory containing input files.

OUTPUT_DIR="/well/hinch/projects/Brca2/PacBio/Trinity/results/variants/sniffles"
# Path to the directory where output files will be saved.

# Loop over each BAM file in the input directory
for bam_file in ${INPUT_DIR}/*.bam; do
    # This loop iterates over all .bam files in the INPUT_DIR.

    # Extract the base filename without extension to use as sample name
    sample_name=$(basename $bam_file .sorted.bam)
    # Extracts the sample name from the filename, removing the path and the '.sorted.bam' suffix.

    # Run sniffles to generate SNF file for each sample
    sniffles -i $bam_file --snf ${OUTPUT_DIR}/${sample_name}.snf --minsvlen 100 --output-rnames #SVs of length >= 100bp
    # Runs the sniffles command on each .bam file, generating a .snf file in the OUTPUT_DIR with the same sample name.

    # Echo completion of current sample
    echo "Sniffles SV calling for $sample_name complete."
    # Prints a message indicating completion of structural variant calling for the current sample.
done

echo "Sniffles SV calling for all samples complete."
# This line prints a message indicating that the script has processed all samples.
