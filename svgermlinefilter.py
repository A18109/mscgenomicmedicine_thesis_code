import sys

def read_bed(file):
    """
    Reads the BED file and returns a list of dictionaries, each representing a BED entry.

    Parameters:
    file (str): Path to the BED file

    Returns:
    list: List of dictionaries with keys 'chr', 'start', 'end', 'svtype', and 'svlen'
    """
    bed_data = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('chr'):
                parts = line.strip().split()
                bed_data.append({
                    "chr": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "svtype": parts[3],
                    "svlen": int(parts[4])
                })
    return bed_data

def read_vcf(file):
    """
    Reads the VCF file and returns the header lines and variant entries separately.

    Parameters:
    file (str): Path to the VCF file

    Returns:
    tuple: (list of header lines, list of variant entries)
    """
    vcf_header = []
    vcf_data = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                vcf_header.append(line.strip())
            else:
                vcf_data.append(line.strip())
    return vcf_header, vcf_data

def filter_vcf(bed_data, vcf_data, delta):
    """
    Filters the VCF data based on the BED entries and specified delta.

    Parameters:
    bed_data (list): List of BED entries
    vcf_data (list): List of VCF variant entries
    delta (int): Delta value for SVLEN comparison

    Returns:
    list: Filtered list of VCF variant entries
    """
    filtered_vcf_data = []
    for variant in vcf_data:
        parts = variant.split('\t')
        chr = parts[0]
        pos = int(parts[1])
        info = parts[7]

        # Extract SVTYPE and SVLEN from the INFO field
        svtype = None
        svlen = None
        for field in info.split(';'):
            if field.startswith('SVTYPE='):
                svtype = field.split('=')[1]
            if field.startswith('SVLEN='):
                svlen = int(field.split('=')[1])

        # If SVTYPE or SVLEN is not found, skip this variant
        if svtype is None or svlen is None:
            continue

        exclude = False
        for bed_entry in bed_data:
            # Check if the variant matches the criteria for exclusion
            if (chr == bed_entry['chr'] and
                svtype == bed_entry['svtype'] and
                abs(svlen) >= abs(bed_entry['svlen']) - 10 and
                abs(svlen) <= abs(bed_entry['svlen']) + 10 and
                bed_entry['start'] <= pos <= bed_entry['end']):
                exclude = True
                break

        # If the variant does not match any exclusion criteria, add it to the filtered list
        if not exclude:
            filtered_vcf_data.append(variant)

    return filtered_vcf_data

def main(bed_file, vcf_file, output_file):
    """
    Main function to read the BED and VCF files, filter the VCF data, and write the output.

    Parameters:
    bed_file (str): Path to the BED file
    vcf_file (str): Path to the VCF file
    output_file (str): Path to the output VCF file
    """
    # Read the BED and VCF files
    bed_data = read_bed(bed_file)
    vcf_header, vcf_data = read_vcf(vcf_file)
    # Filter the VCF data
    filtered_vcf_data = filter_vcf(bed_data, vcf_data, 10)

    # Write the filtered VCF data to the output file
    with open(output_file, 'w') as f:
        for line in vcf_header:
            f.write(line + '\n')
        for variant in filtered_vcf_data:
            f.write(variant + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python filter_vcf.py <bed_file> <vcf_file> <output_file>")
        sys.exit(1)

    bed_file = sys.argv[1]
    vcf_file = sys.argv[2]
    output_file = sys.argv[3]

    main(bed_file, vcf_file, output_file)


