import sys

def read_vcf(path):
    """
    Read a VCF file and extract the data lines, skipping header lines.
    """
    with open(path, 'r') as file:
        data = [line.strip().split('\t') for line in file if not line.startswith('#')]
    return data

def extract_read_names(vcf_data, output_path):
    """
    Extract read names supporting each SV and save to a text file.
    """
    with open(output_path, 'w') as file:
        for record in vcf_data:
            info = record[7]
            rnames = [item for item in info.split(';') if item.startswith('RNAMES=')][0].split('=')[1]
            file.write(rnames + '\n')

def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_rnames.py <input_vcf> <output_txt>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_txt = sys.argv[2]

    vcf_data = read_vcf(input_vcf)
    extract_read_names(vcf_data, output_txt)

    print(f"Read names file created at {output_txt}")

if __name__ == "__main__":
    main()