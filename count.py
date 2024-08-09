import sys

def read_bed_file(bed_file):
    """
    Reads a BED file and returns a list of regions.

    Each region is represented as a tuple (chromosome, start, end).

    Args:
        bed_file (str): Path to the input BED file.

    Returns:
        List[Tuple[str, int, int]]: List of regions from the BED file.
    """
    bed_regions = []
    with open(bed_file, 'r') as bed:
        for line in bed:
            chrom, start, end = line.strip().split()
            bed_regions.append((chrom, int(start), int(end)))
    return bed_regions

def count_snps_in_region(vcf_file, bed_regions):
    """
    Counts the number of SNPs with genotypes 0/0, 0/1, and 1/1 within each region specified in the BED file.

    Args:
        vcf_file (str): Path to the input VCF file.
        bed_regions (List[Tuple[str, int, int]]): List of regions from the BED file.

    Returns:
        List[Tuple[str, int, int, int, int, int]]: List of regions with counts of each genotype.
    """
    region_snp_counts = []

    for chrom, start, end in bed_regions:
        # Initialize genotype counts
        counts = {'0/0': 0, '0/1': 0, '1/1': 0}

        with open(vcf_file, 'r') as vcf:
            for line in vcf:
                if line.startswith('#'):
                    continue

                columns = line.strip().split('\t')
                vcf_chrom = columns[0]
                pos = int(columns[1])
                gt = columns[9].split(':')[0]

                # Check if the SNP is within the current region
                if vcf_chrom == chrom and start <= pos <= end:
                    if gt in counts:
                        counts[gt] += 1

        # Append the counts for the current region
        region_snp_counts.append((chrom, start, end, counts['0/0'], counts['0/1'], counts['1/1']))

    return region_snp_counts

def write_tsv(output_tsv, region_snp_counts):
    """
    Writes the SNP counts within each region to a TSV file.

    Args:
        output_tsv (str): Path to the output TSV file.
        region_snp_counts (List[Tuple[str, int, int, int, int, int]]): List of regions with counts of each genotype.
    """
    with open(output_tsv, 'w') as tsv:
        # Write the header
        tsv.write("Chromosome\tStart\tEnd\tGT_0/0\tGT_0/1\tGT_1/1\n")

        # Write the SNP counts for each region
        for chrom, start, end, gt_0_0, gt_0_1, gt_1_1 in region_snp_counts:
            tsv.write(f"{chrom}\t{start}\t{end}\t{gt_0_0}\t{gt_0_1}\t{gt_1_1}\n")

if __name__ == '__main__':
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python vcf_bed_to_tsv.py <input_vcf> <input_bed> <output_tsv>")
        sys.exit(1)

    # Parse command line arguments
    input_vcf = sys.argv[1]
    input_bed = sys.argv[2]
    output_tsv = sys.argv[3]

    # Read BED file and get regions
    bed_regions = read_bed_file(input_bed)

    # Count SNPs in each region
    region_snp_counts = count_snps_in_region(input_vcf, bed_regions)

    # Write results to TSV file
    write_tsv(output_tsv, region_snp_counts)
