import sys

def parse_info(info):
    info_dict = {}
    for item in info.split(';'):
        key_value = item.split('=')
        if len(key_value) == 2:
            key, value = key_value
            info_dict[key] = value
    return info_dict

def process_vcf(input_vcf, output_bed):
    with open(input_vcf, 'r') as vcf, open(output_bed, 'w') as bed:
        for line in vcf:
            if line.startswith('#'):
                continue

            columns = line.strip().split('\t')
            chrom = columns[0]
            pos = int(columns[1])
            info = columns[7]

            info_dict = parse_info(info)
            svtype = info_dict.get('SVTYPE')
            end = int(info_dict.get('END', pos))
            svlen = int(info_dict.get('SVLEN', 0))

            if svtype == 'DEL':
                start = max(0, pos - 14280)
                stop = end + 14280
            elif svtype == 'INS':
                start = max(0, pos - 14280)
                stop = pos + abs(svlen) + 14280
            else:
                continue

            bed.write(f'{chrom}\t{start}\t{stop}\n')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python vcf_to_bed.py <input_vcf> <output_bed>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_bed = sys.argv[2]

    process_vcf(input_vcf, output_bed)
