
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

def parse_info_field(info_field):
    """
    Parse the INFO field of a VCF record into a dictionary.
    Args:
    info_field (str): The INFO field string.

    Returns:
    dict: Dictionary with INFO field keys and values.
    """
    info_dict = {}
    for entry in info_field.split(';'):
        if '=' in entry:
            key, value = entry.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[entry] = True
    return info_dict

def filter_de_novo_mutations(vcf_data):
    """
    Filters the VCF data to identify potential de novo mutations in pup1.
    """
    de_novo_mutations = []
    for record in vcf_data:
        # Extract the genotype information from the FORMAT and respective samples
        format_index = record[8].split(':').index('GT')  # Find the index of 'GT' in FORMAT
        father_gt = record[9].split(':')[format_index]
        pup1_gt = record[10].split(':')[format_index]
        pup2_gt = record[11].split(':')[format_index]

        # Parse the INFO field
        info_dict = parse_info_field(record[7])

        # Get SVLEN and SUPPORT values
        svlen = abs(int(info_dict.get('SVLEN', 0)))
        support = int(info_dict.get('SUPPORT', 0))

        # Filter based on the given criteria
        if svlen > 14280:
            if pup2_gt in ['1/1', '0/1'] and father_gt in ['0/0', './.'] and pup1_gt in ['0/0', './.'] and support >= 10:
                de_novo_mutations.append(record)
        else:
            if pup2_gt in ['0/1'] and pup1_gt in ['./.', '0/0'] and father_gt in ['./.', '0/0'] and support >= 10:
                de_novo_mutations.append(record)

    return de_novo_mutations

def save_filtered_vcf(metadata, filtered_data, output_path):
    """
    Saves the filtered VCF data to a new file, including metadata headers.
    """
    with open(output_path, 'w') as file:
        for line in metadata:
            file.write(line + '\n')
        file.write('#' + '\t'.join(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'father', 'pup1', 'pup2']) + '\n')
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

    # Filter the VCF data for de novo mutations
    de_novo_mutations = filter_de_novo_mutations(vcf_data)

    # Save the filtered data, including metadata, to a new VCF file
    save_filtered_vcf(metadata, de_novo_mutations, output_vcf)

    print(f"Filtered VCF saved to {output_vcf}")

if __name__ == "__main__":
    main()