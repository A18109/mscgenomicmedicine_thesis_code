import sys
import gzip
import os
import math

def load_vcf(file_path, chromosome):
    if not os.path.exists(file_path):
        print(f"Error: The file {file_path} does not exist.", file=sys.stderr)
        sys.exit(1)
    
    snps = {}
    header = []
    with gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                header.append(line)
                continue
            parts = line.strip().split('\t')
            if parts[0] != chromosome:
                continue
            pos = int(parts[1])
            format_fields = parts[8].split(':')
            sample_fields = parts[9].split(':')

            if 'UG' in format_fields:
                ug_index = format_fields.index('UG')
                unphased_genotype = sample_fields[ug_index]
                snps[pos] = {
                    'unphased_genotype': unphased_genotype,
                    'line': line.strip()
                }
            elif 'GT' in format_fields:
                gt_index = format_fields.index('GT')
                unphased_genotype = sample_fields[gt_index]
                snps[pos] = {
                    'unphased_genotype': unphased_genotype,
                    'line': line.strip()
                }
            else:
                print(f"Error: 'UG' or 'GT' field not found in format for position {pos} in {file_path}", file=sys.stderr)
    
    print(f"Loaded {len(snps)} positions from {file_path}", file=sys.stderr)
    sample_positions = list(snps.keys())[:10] + list(snps.keys())[-10:]
    for pos in sample_positions:
        print(f"Sample SNP position from {file_path}: {pos} -> {snps[pos]['unphased_genotype']}", file=sys.stderr)
    
    return snps, header

def initialize_hmm_parameters():
    states = ['B6', 'CAST']

    transition_prob = {
        'B6': {'B6': 0.995, 'CAST': 0.005},
        'CAST': {'B6': 0.005, 'CAST': 0.995}
    }

    emission_prob = {
        'B6': {'equal': 0.98, 'not_equal': 0.02},
        'CAST': {'equal': 0.02, 'not_equal': 0.98}
    }

    log_transition_prob = {s: {s2: math.log(p) for s2, p in sp.items()} for s, sp in transition_prob.items()}
    log_emission_prob = {s: {o: math.log(p) for o, p in sp.items()} for s, sp in emission_prob.items()}

    return states, log_transition_prob, log_emission_prob

def get_emission_log_prob(state, observation, emit_prob):
    return emit_prob[state].get(observation, float('-inf'))

def viterbi(observed_sequence, states, start_prob, trans_prob, emit_prob):
    V = [{}]
    path = {}

    for state in states:
        V[0][state] = math.log(start_prob[state]) + get_emission_log_prob(state, observed_sequence[0], emit_prob)
        path[state] = [state]
    print(f"Initial log probabilities: {V[0]}", file=sys.stderr)

    for t in range(1, len(observed_sequence)):
        V.append({})
        new_path = {}

        for state in states:
            max_log_prob, prev_state = max(
                (V[t-1][prev_state] + trans_prob[prev_state][state] + get_emission_log_prob(state, observed_sequence[t], emit_prob), prev_state)
                for prev_state in states
            )
            V[t][state] = max_log_prob
            new_path[state] = path[prev_state] + [state]
        path = new_path
        if t % 100000 == 0:
            print(f"Step {t}: {V[t]}", file=sys.stderr)

    n = len(observed_sequence) - 1
    max_log_prob, state = max((V[n][state], state) for state in states)
    print(f"Final log probabilities: {V[n]}", file=sys.stderr)
    return path[state]

def identify_b6_positions(pup1_snps, most_likely_states):
    b6_positions = []

    for pos, state in zip(sorted(pup1_snps.keys()), most_likely_states):
        if state == 'B6':
            b6_positions.append(pos)

    return b6_positions

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python3 identify_b6_positions.py <B6_father_vcf> <pup1_vcf> <pup2_vcf> <chromosome> <known_cast_vcf>", file=sys.stderr)
        sys.exit(1)
    
    b6_father_vcf = sys.argv[1]
    pup1_vcf = sys.argv[2]
    pup2_vcf = sys.argv[3]
    chromosome = sys.argv[4]
    known_cast_vcf = sys.argv[5]

    print(f"Processing chromosome: {chromosome}", file=sys.stderr)

    b6_snps, _ = load_vcf(b6_father_vcf, chromosome)
    pup1_snps, _ = load_vcf(pup1_vcf, chromosome)
    pup2_snps, _ = load_vcf(pup2_vcf, chromosome)
    known_cast_snps, _ = load_vcf(known_cast_vcf, chromosome)

    print(f"Loaded {len(b6_snps)} SNPs from B6 father VCF", file=sys.stderr)
    print(f"Loaded {len(pup1_snps)} SNPs from pup1 VCF", file=sys.stderr)
    print(f"Loaded {len(pup2_snps)} SNPs from pup2 VCF", file=sys.stderr)
    print(f"Loaded {len(known_cast_snps)} SNPs from known CAST VCF", file=sys.stderr)

    states, trans_prob, emit_prob = initialize_hmm_parameters()

    observed_sequence = []
    for pos in sorted(pup1_snps.keys()):
        pup1_ug = pup1_snps[pos]['unphased_genotype']
        b6_ug = b6_snps[pos]['unphased_genotype'] if pos in b6_snps else None
        pup2_ug = pup2_snps[pos]['unphased_genotype'] if pos in pup2_snps else None
        known_cast = known_cast_snps.get(pos, {}).get('unphased_genotype')

        print(f"Position: {pos}, pup1_ug: {pup1_ug}, b6_ug: {b6_ug}, pup2_ug: {pup2_ug}, known_cast: {known_cast}", file=sys.stderr)

        if pup1_ug == '0/1':
            if (b6_ug is None) and (pup2_ug == '0/1'):
                if known_cast == '1/1':
                    observed_sequence.append('not_equal')
                else:
                    observed_sequence.append('equal')
            elif (b6_ug is None or b6_ug == '0/0') and (pup2_ug is None or pup2_ug == '0/0'):
                if known_cast == '1/1':
                    observed_sequence.append('not_equal')
                else:
                    observed_sequence.append('equal')
            else:
                observed_sequence.append('equal')
        elif pup1_ug == '1/1':
            if pup2_ug in {'0/1', '1/1'} and known_cast == '1/1':
                observed_sequence.append('not_equal')
            elif pup2_ug in {'0/1', '1/1'} and known_cast != '1/1':
                observed_sequence.append('equal')
            else:
                observed_sequence.append('equal')
        elif pup1_ug == '0/0':
            observed_sequence.append('equal')

        if pos in {4073346, 3194487, 3175966}:
            print(f"Observation appended for position {pos}: {observed_sequence[-1]}", file=sys.stderr)
    
    print(f"Total positions in pup1: {len(pup1_snps)}", file=sys.stderr)
    print(f"Observed sequence: {observed_sequence[:1000]}", file=sys.stderr)

    start_prob = {state: 1/len(states) for state in states}

    most_likely_states = viterbi(observed_sequence, states, start_prob, trans_prob, emit_prob)
    
    print(f"Most likely states: {most_likely_states[:1000]}", file=sys.stderr)
    print(f"Total most likely states: {len(most_likely_states)}", file=sys.stderr)

    b6_positions = identify_b6_positions(pup1_snps, most_likely_states)
    
    print(f"B6 positions: {b6_positions[:1000]}", file=sys.stderr)
    print(f"Total B6 positions: {len(b6_positions)}", file=sys.stderr)

    sample_positions = b6_positions[:10] + b6_positions[-10:]
    for pos in sample_positions:
        print(f"Sample B6 position: {chromosome}\t{pos}", file=sys.stderr)

    for pos in b6_positions:
        print(f"{chromosome}\t{pos}")

    if not b6_positions:
        print("No positions with 'B6' state identified.", file=sys.stderr)
    else:
        print(f"Identified {len(b6_positions)} positions with 'B6' state.", file=sys.stderr)
