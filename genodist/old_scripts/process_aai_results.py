import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process BLASTP outputs and calculate AAI (Average Amino acid Identity).')
    parser.add_argument('--blast1', required=True, help='Forward BLASTP output file')
    parser.add_argument('--blast2', required=True, help='Reverse BLASTP output file')
    parser.add_argument('--output', required=True, help='Output file path')
    parser.add_argument('--genome1', required=True, help='Name of the first genome')
    parser.add_argument('--genome2', required=True, help='Name of the second genome')
    parser.add_argument('--identity', type=float, default=30.0, help='Minimum percent identity threshold (default: 30.0)')
    parser.add_argument('--coverage', type=float, default=70.0, help='Minimum alignment coverage threshold (default: 70.0)')
    return parser.parse_args()

def parse_blastp_output(blast_output, min_identity=30, min_coverage=70):
    """
    Parse BLASTP output file and return a list of hits that pass the identity and coverage thresholds.
    Input: BLASTP output file path. Format: "-outfmt "6 qseqid sseqid pident length qlen slen"
    Output: list of hits
    """
    hits = []
    with open(blast_output) as f:
        for line in f:
            parts = line.strip().split("\t")
            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])  # percent identity
            align_length = int(parts[3])
            qlen = int(parts[4])
            slen = int(parts[5])

            # Calculate coverage
            q_coverage = align_length / qlen * 100
            s_coverage = align_length / slen * 100

            if q_coverage >= min_coverage and s_coverage >= min_coverage:  # apply coverage filter
                if pident >= min_identity:  # apply identity filter
                    hits.append({
                        "query": qseqid,
                        "subject": sseqid,
                        "identity": pident,
                        "alignment_length": align_length,
                    })
    
    return hits

def find_bbh(bbh1, bbh2):
    best_hits1 = {a['query']: a for a in bbh1}  # dictionary of best hits from bbh1
    best_hits2 = {a['query']: a for a in bbh2}  # dictionary of best hits from bbh2
    bbh = [best_hits1[query] for query in best_hits1 if best_hits1[query]['subject'] in best_hits2 and best_hits2[best_hits1[query]['subject']]['subject'] == query]
    return bbh

def calculate_aai(bbh, query_count):
    """
    Calculate AAI with more comprehensive metrics
    
    Args:
        bbh: List of bidirectional best hits
        query_count: Total number of query genes (for normalization)
    
    Returns:
        Dictionary with AAI metrics
    """
    num_hits = len(bbh)
    if num_hits == 0:
        return {
            'aai': 0,
            'std_dev': 0,
            'hit_count': 0,
            'query_count': query_count,
            'hit_percentage': 0
        }
    
    total_identity = sum([h['identity'] for h in bbh])
    mean_identity = total_identity / num_hits
    variance = sum([(h['identity'] - mean_identity) ** 2 for h in bbh]) / num_hits
    std_dev = variance ** 0.5
    hit_percentage = (num_hits / query_count) * 100 if query_count > 0 else 0
    
    return {
        'aai': mean_identity,
        'std_dev': std_dev,
        'hit_count': num_hits,
        'query_count': query_count,
        'hit_percentage': hit_percentage
    }

def main():
    # Parse command-line arguments
    args = parse_arguments()
    
    # Get parameters from command-line arguments
    forward_blast = args.blast1
    reverse_blast = args.blast2
    output_file = args.output
    genome1 = args.genome1
    genome2 = args.genome2
    min_identity = args.identity
    min_coverage = args.coverage
    
    # Parse BLAST outputs
    bbh1 = parse_blastp_output(forward_blast, min_identity, min_coverage)
    bbh2 = parse_blastp_output(reverse_blast, min_identity, min_coverage)
    
    # Find bidirectional best hits
    bbh = find_bbh(bbh1, bbh2)
    
    # Count total genes in each genome
    gene_count1 = len(set(h['query'] for h in bbh1))
    gene_count2 = len(set(h['query'] for h in bbh2))

    # Calculate AAI
    results = calculate_aai(bbh, min(gene_count1, gene_count2))
    
    # Write AAI value to output file
    # Write header and AAI results to output file
    with open(output_file, 'w') as f:
        f.write(f'{genome1}\t{genome2}\t'
            f'{results["aai"]:.2f}\t{results["std_dev"]:.2f}\t'
            f'{results["hit_count"]}\t{results["query_count"]}\t'
            f'{results["hit_percentage"]:.2f}\n')

if __name__ == "__main__":
    main()
