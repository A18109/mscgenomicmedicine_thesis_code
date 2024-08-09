import sys

def read_vcf(path):
    """
    Read a VCF file and extract the data lines, skipping the header lines.
    This function assumes that VCF header lines start with '#'.
    Returns a list of list of strings where each sublist represents a VCF record.
    """
    with open(path, 'r') as file:
        data = [line.strip().split('\t') for line in file if not line.startswith('#')]
    return data

def create_bed(vcf_data, output_path):
    """
    Create a BED file for specific regions around each SV based on its type.

    Args:
    vcf_data: List of VCF records obtained from read_vcf.
    output_path: Path to write the output BED file.

    For 'DEL' type:
    - Extracts two regions: (SVstart - 14.28kb, SVstart) and (SVend, SVend + 14.28kb).

    For 'INS' type:
    - Extracts two regions: (SVstart - 14.28kb, SVstart) and (SVstart + SVlen, SVstart + SVlen + 14.28kb).
    """
    with open(output_path, 'w') as file:
        for record in vcf_data:
            chrom = record[0]
            pos = int(record[1])
            info = record[7]
            svtype = [item for item in info.split(';') if item.startswith('SVTYPE=')][0].split('=')[1]
            svlen = int([item for item in info.split(';') if item.startswith('SVLEN=')][0].split('=')[1])
            end = int([item for item in info.split(';') if item.startswith('END=')][0].split('=')[1])

            padding = 14280  # Define the size of the region to extract around the SV points 

            if svtype == "DEL":
                # For deletions, handle regions before the start and after the end
                start1 = max(0, pos - padding)
                end1 = pos
                start2 = end
                end2 = end + padding
            elif svtype == "INS":
                # For insertions, handle regions before the start and the adjusted end by the length of the insertion
                start1 = max(0, pos - padding)
                end1 = pos
                start2 = pos + svlen
                end2 = start2 + padding

            # Write each region to the BED file
            file.write(f"{chrom}\t{start1}\t{end1}\n")
            file.write(f"{chrom}\t{start2}\t{end2}\n")

def main():
    """
    Main function to execute the script logic. Requires two command line arguments:
    - Input VCF file path
    - Output BED file path
    """
    if len(sys.argv) != 3:
        print("Usage: python create_bed.py <input_vcf> <output_bed>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_bed = sys.argv[2]

    vcf_data = read_vcf(input_vcf)
    create_bed(vcf_data, output_bed)

    print(f"Bed file created at {output_bed}")

if __name__ == "__main__":
    main()