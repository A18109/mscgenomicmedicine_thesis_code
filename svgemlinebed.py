
import sys
delta = 50

def read_vcf(path):
    """
    Read a VCF file and extract data.
    Args:
    path (str): Path to the VCF file.

    Returns:
    list: List containing the parsed VCF data split into fields.
    """
    data = []
    with open(path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                data.append(line.strip().split('\t'))
    return data

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

def vcf_to_bed(vcf_data, output_path):
    """
    Convert VCF data to BED format and save to a file.
    Args:
    vcf_data (list): List of VCF records.
    output_path (str): Path to the output BED file.
    """
    with open(output_path, 'w') as file:
        for record in vcf_data:
            chrom = record[0]
            pos = int(record[1])
            info_dict = parse_info_field(record[7])
            svtype = info_dict.get('SVTYPE', 'NA')
            svlen = info_dict.get('SVLEN', 'NA')
            start = pos - delta
            end = pos + delta
            file.write(f"{chrom}\t{start}\t{end}\t{svtype}\t{svlen}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_vcf> <output_bed>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_bed = sys.argv[2]

    # Read the VCF file
    vcf_data = read_vcf(input_vcf)

    # Convert the VCF data to BED format
    vcf_to_bed(vcf_data, output_bed)

    print(f"BED file saved to {output_bed}")

if __name__ == "__main__":
    main()