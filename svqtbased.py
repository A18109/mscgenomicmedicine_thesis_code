import sys

def read_vcf(path):
    """
    Read a VCF file and extract data, including metadata lines starting with '##'.
    Args:
    path (str): Path to the VCF file.

    Returns:
    tuple: (list of metadata lines, list of lists containing the parsed VCF data split into fields)
    """
    metadata = []
    data = []
    with open(path, 'r') as file:
        for line in file:
            if line.startswith('##'):
                metadata.append(line.strip())
            elif not line.startswith('#'):
                data.append(line.strip().split('\t'))
    return (metadata, data)

def filter_vcf(vcf_data):
    """
    Filters the VCF data based on specified criteria.
    """
    filtered_data = []
    for record in vcf_data:
        chrom, pos, id, ref, alt, qual, filter_, info = record[:8]

        # Filter by chromosome
        if not chrom.startswith('chr') or not chrom[3:].isdigit() or not (1 <= int(chrom[3:]) <= 19):
            continue

        # Convert QUAL to float and check filter and quality
        try:
            qual = float(qual)
            if filter_ != 'PASS' or qual < 60:
                continue
        except ValueError:
            continue  # Skip records with non-numeric QUAL values

        # Parsing the INFO field into a dictionary
        info_dict = dict(item.split('=') for item in info.split(';') if '=' in item)

        # Check for required SVTYPE, STRAND, and SUPPORT
        svtype = info_dict.get('SVTYPE')
        strands = info_dict.get('STRAND')
        support = int(info_dict.get('SUPPORT', '0'))

        if svtype not in ['INS', 'DEL'] or strands != '+-':
            continue

        # If all conditions are met, add to filtered data
        filtered_data.append(record)

    return filtered_data

def save_filtered_vcf(metadata, filtered_data, output_path):
    """
    Saves the filtered VCF data to a new file, including metadata headers.
    """
    with open(output_path, 'w') as file:
        for line in metadata:
            file.write(line + '\n')
        file.write('#' + '\t'.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'father', 'pup1', 'pup2']) + '\n')  # Write the column headers
        for record in filtered_data:
            file.write('\t'.join(record) + '\n')

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_vcf> <output_vcf>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]

    # Read the VCF file, including metadata
    metadata, vcf_data = read_vcf(input_vcf)

    # Filter the VCF data
    filtered_data = filter_vcf(vcf_data)

    # Save the filtered data, including metadata, to a new VCF file
    save_filtered_vcf(metadata, filtered_data, output_vcf)

    print(f"Filtered VCF saved to {output_vcf}")

if __name__ == "__main__":
    main()