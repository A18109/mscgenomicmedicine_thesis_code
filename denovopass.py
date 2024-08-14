import sys
import gzip

def is_gzipped(file_path):
    """
    Check if a file is gzipped.

    Parameters:
    file_path (str): Path to the file.

    Returns:
    bool: True if the file is gzipped, False otherwise.
    """
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def filter_pass_variants(input_vcf, output_vcf):
    """
    Filter variants with 'PASS' in the FILTER column and QUAL value >= 500 from a VCF file.

    Parameters:
    input_vcf (str): Path to the input VCF file.
    output_vcf (str): Path to the output VCF file.
    """
    if is_gzipped(input_vcf):
        open_func = gzip.open
    else:
        open_func = open

    with open_func(input_vcf, 'rt') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                # Write the header lines to the output file
                outfile.write(line)
                continue

            fields = line.strip().split("\t")
            filter_value = fields[6]
            qual_value = float(fields[5])

            # Only write the line if the FILTER value is 'PASS' and QUAL value >= 500
            if filter_value == 'PASS' and qual_value >= 500:
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_vcf_pass.py <input_vcf> <output_vcf>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]

    filter_pass_variants(input_vcf, output_vcf)
