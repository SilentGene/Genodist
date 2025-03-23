import os
import sys
import argparse
from subprocess import check_call, CalledProcessError
import hashlib  # for generating unique sequence IDs

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate Average Amino Acid Identity (AAI) between two genomes')
    parser.add_argument('-1', '--genome1', required=True, help='Query genome protein sequence in FASTA format')
    parser.add_argument('-2', '--genome2', required=True, help='Subject genome protein sequence in FASTA format')
    parser.add_argument('-o', '--output', required=True, help='Output tsv file')
    parser.add_argument('-tmp', '--tmpdir', default='.', help='Temporary directory to store intermediate files')
    parser.add_argument('-n', '--threads', type=int, default=1, help='Number of threads to use for BLAST')
    parser.add_argument('-p', '--program', choices=['blast', 'diamond'], default='blast', help='Program to use for alignment')
    parser.add_argument('-id', '--identity', type=float, default=30, help='Minimum identity percentage to consider a hit')
    parser.add_argument('-cov', '--coverage', type=float, default=70, help='Minimum coverage percentage to consider a hit')
    return parser.parse_args()

def rename_fasta_headers(fasta_file, output_file):
    # file name hash 32 characters long
    file_hash = hashlib.md5(fasta_file.encode("utf-8")).hexdigest()
    
    with open(fasta_file) as f, open(output_file, 'w') as out:
        n = 0
        for i, line in enumerate(f):
            if line.startswith('>'):
                out.write(f'>{file_hash}_{n+1}\n')
                n += 1
            else:
                out.write(line)
    return output_file

def run_makeblastdb_prot(faa_file):

    makeblastdb_cmd = f"makeblastdb -dbtype prot -in {faa_file}"

    try:
        check_call(makeblastdb_cmd, shell=True)
    except CalledProcessError as e:
        print(f"Error executing command: {makeblastdb_cmd}", file=sys.stderr)
        sys.exit(e.returncode)
    

def run_blastp(query, db, output, threads):
    blast_cmd = (
        f"blastp -query {query} -db {db} "
        f"-evalue 1e-5 -outfmt '6 qseqid sseqid pident length qlen slen' -out {output} -num_threads {threads} -max_target_seqs 1"
    )
    try:
        check_call(blast_cmd, shell=True)
    except CalledProcessError as e:
        print(f"Error executing command: {blast_cmd}", file=sys.stderr)
        sys.exit(e.returncode)

def parse_blastp_output(blast_output, min_identity=30, min_coverage=70):
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

            if q_coverage >= min_identity and s_coverage >= min_identity:  # apply coverage filter
                if pident >= min_identity:  # apply identity filter
                    hits.append({
                            "query": qseqid,
                            "subject": sseqid,
                            "identity": pident,
                            "alignment_length": align_length,
                        })
    
    return hits

def find_bbh(bbh1, bbh2):
    best_hits1 = {a['query']: a for a in bbh1}
    best_hits2 = {a['query']: a for a in bbh2}
    bbh = [a for a in bbh1 if a['subject'] in best_hits2 and best_hits2[a['subject']]['subject'] == a['query']]
    return bbh

def calculate_aai(bbh):
    num_hits = len(bbh)
    print(f'Number of BBH: {num_hits}')
    if num_hits == 0:
        return 0, 0
    
    total_identity = sum([h['identity'] for h in bbh])
    mean_identity = total_identity / num_hits
    variance = sum([(h['identity'] - mean_identity) ** 2 for h in bbh]) / num_hits
    std_dev = variance ** 0.5
    return mean_identity, std_dev
    

def main():
    args = parse_arguments()

    # Create output directory if it doesn't exist
    os.makedirs(args.tmpdir, exist_ok=True)

    # Rename headers in FASTA files
    genome1_renamed = rename_fasta_headers(args.genome1, os.path.join(args.tmpdir, os.path.basename(args.genome1) + '.renamed'))
    genome2_renamed = rename_fasta_headers(args.genome2, os.path.join(args.tmpdir, os.path.basename(args.genome2) + '.renamed'))

    # Create BLAST database for both genomes
    run_makeblastdb_prot(genome1_renamed)
    run_makeblastdb_prot(genome2_renamed)

    # Run BLAST+ for both directions
    blast_output1 = os.path.join(args.tmpdir, f'{os.path.basename(args.genome1)}_vs_{os.path.basename(args.genome2)}.blast')
    blast_output2 = os.path.join(args.tmpdir, f'{os.path.basename(args.genome2)}_vs_{os.path.basename(args.genome1)}.blast')

    run_blastp(genome1_renamed, genome2_renamed, blast_output1, args.threads)
    run_blastp(genome2_renamed, genome1_renamed, blast_output2, args.threads)

    # Parse BLAST outputs
    bbh1 = parse_blastp_output(blast_output1)  # g1 vs g2
    bbh2 = parse_blastp_output(blast_output2)  # g2 vs g1

    # find bbh
    bbh = find_bbh(bbh1, bbh2)

    # Calculate AAI using BBH results
    aai, aai_dev = calculate_aai(bbh)

    # Print results
    print(f'Average Amino Acid Identity (AAI): {aai:.2f} ± {aai_dev:.2f}')

    # Write AAI value
    with open(args.output, 'w') as afo:
        afo.write(f'{os.path.basename(args.genome1)}\t{os.path.basename(args.genome2)}\t{aai:.2f}\t{aai_dev:.2f}\n')

if __name__ == "__main__":
    main()
