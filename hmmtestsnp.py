import sys
import gzip
# this script takes in a bed file and a vcf file and outputs a new vcf file with positions defined in the bed file
def read_bed_file(bed_file):
    positions = set()
    with open(bed_file, 'r') as bf:
        for line in bf:
            chrom, pos = line.strip().split()
            positions.add((chrom, int(pos)))
    return positions

def extract_variants(vcf_file, positions, output_file):
    with gzip.open(vcf_file, 'rt') as vf, gzip.open(output_file, 'wt') as of:
        for line in vf:
            if line.startswith('#'):
                of.write(line)
            else:
                parts = line.strip().split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                if (chrom, pos) in positions:
                    of.write(line)

def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_variants.py <vcf_file> <bed_file> <output_file>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    bed_file = sys.argv[2]
    output_file = sys.argv[3]

    positions = read_bed_file(bed_file)
    extract_variants(vcf_file, positions, output_file)

if __name__ == "__main__":
    main()
