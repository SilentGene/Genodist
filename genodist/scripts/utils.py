import sys
import hashlib

def split_genome(fasta_file, chunk_file, chop_len=1000, min_len=700, step=200):
    """
    Split the query genome into smaller segments and write them to a file.
    
    Args:
        query_file (str): Path to the input FASTA file
        chunk_file (str): Path to the output file for genome chunks
        chop_len (int): Length of each chunk (default: 1000)
        min_len (int): Minimum length for valid chunks (default: 700)
        step (int): Step size for sliding window (default: 200)
    
    Returns:
        int: Number of chunks generated
    """
    # Initialize chunk counter
    chunk_count = 0

    # Generate a unique hash for the output chunks
    file_hash = hashlib.md5(fasta_file.encode("utf-8")).hexdigest()
    
    # Open output file
    with open(chunk_file, 'w') as out_f:
        # Process each sequence in the FASTA file
        with open(fasta_file) as in_f:
            n = 0
            seq = []
            
            # Read through the FASTA file
            for line in in_f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    # Process the previous sequence if it exists
                    if seq:
                        sequence = ''.join(seq)
                        # Generate chunks for the sequence
                        for i in range(0, len(sequence)-chop_len+1, step):
                            chunk = sequence[i:i+chop_len]
                            # Only write chunks that meet the minimum length
                            if len(chunk) >= min_len:
                                chunk_count += 1
                                # Write chunk to output file
                                out_f.write(f'>{file_hash}_{n+1}\n')
                                out_f.write(f'{chunk}\n')
                                n += 1
                    
                    # Start new sequence
                    seq = []
                else:
                    seq.append(line)
            
            # Process the last sequence
            if seq:
                sequence = ''.join(seq)
                for i in range(0, len(sequence)-chop_len+1, step):
                    chunk = sequence[i:i+chop_len]
                    if len(chunk) >= min_len:
                        chunk_count += 1
                        out_f.write(f'>{file_hash}_{n+1}\n')
                        out_f.write(f'{chunk}\n')
                        n += 1


def cal_ani(blast_output, len_cut=100):
    """
    Parse BLAST output and calculate ANI.
    """

    # Set Identity and Alignment Percentage cut-off values
    cvg_cut = 70
    sum_id = 0
    count = 0
    qr_best = {}
    ani = 0

   # Process BLAST output to calculate ANI
    with open(blast_output, 'r') as bl:
        for line in bl:
            cols = line.strip().split('\t')
            qseqid, pident, length = cols[0], float(cols[2]), int(cols[3])
            if length < len_cut:
                continue
            if qseqid in qr_best:
                continue  # only use best hit for every query segment
            qr_best[qseqid] = True
            if (length * 100 / 1020) < cvg_cut:
                continue
            sum_id += pident
            count += 1
    
    if count == 0:
        print(f"No valid hits found in BLAST result {blast_output}. ANI value is set to 0", file=sys.stderr)
        ani = 0
    else:
        ani = sum_id / count
    
    return ani


def write_value(outfile, genome1, genome2, value):
    with open(outfile, 'a') as f:
        f.write(f'{genome1}\t{genome2}\t{value:2f}\n')


def rename_fasta_headers(fasta_file, output_file):
    """
    Rename headers in a FASTA file.
    """
    file_hash = hashlib.md5(fasta_file.encode("utf-8")).hexdigest()
    with open(fasta_file) as f, open(output_file, 'w') as out:
        n = 0
        for line in f:
            if line.startswith('>'):
                out.write(f'>{file_hash}_{n+1}\n')
                n += 1
            else:
                out.write(line)

#== For AAI calculation==
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
#== End of AAI calculation ==
