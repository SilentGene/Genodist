#!/usr/bin/env python3


import os
import sys
import argparse
from subprocess import check_call, CalledProcessError
import hashlib  # for generating unique sequence IDs


def parse_arguments():
    """
    Parse command-line arguments for the script.
    """
    parser = argparse.ArgumentParser(description='Calculate Average Nucleotide Identity (ANI) between two genomes')
    parser.add_argument('-1', '--genome1', required=True, help='Query genome DNA sequence in FASTA format')
    parser.add_argument('-2', '--genome2', required=True, help='Subject genome DNA sequence in FASTA format')
    parser.add_argument('-tmp', '--tmpdir', required=True, help='tempory directory to store intermediate files')
    parser.add_argument('-n', '--threads', type=int, default=1, help='Number of threads to use for BLAST')
    parser.add_argument('-p', '--program', choices=['blast', 'diamond'], default='blast', help='Program to use for alignment')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument("-w", "--win", type=int, default=1000, help="Window size in the ANI calculation (in bp).")
    parser.add_argument("-s", "--step", type=int, default=200, help="Step size in the ANI calculation (in bp).")
    parser.add_argument("-l", "--len", type=int, default=700, help="Minimum contig length (in bp).")
    parser.add_argument("-a", "--aln", type=int, default=50, help="Minimum alignment length (in bp).")
    parser.add_argument("-i", "--id", type=float, default=70, help="Minimum alignment identity (in %).")
    return parser.parse_args()

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
    
    return chunk_count
    

def run_makeblastdb_nucl(fasta_file, output_dir):
    """
    Create a BLAST database from a nucleotide FASTA file.
    """
    db_name = os.path.basename(fasta_file)
    db_path = os.path.join(output_dir, db_name)

    makeblastdb_cmd = f'makeblastdb -dbtype nucl -in {fasta_file} -out {db_path} -parse_seqids -blastdb_version 5'

    try:
        check_call(makeblastdb_cmd, shell=True)
    except CalledProcessError as e:
        print(f"Error executing command: {makeblastdb_cmd}", file=sys.stderr)
        sys.exit(e.returncode)
    
    return db_path

def run_blastn(query, db, output, threads, cut_id=70):
    """
    Run BLAST+ alignment.
    """
    blast_cmd = (
        f'blastn -task blastn -query {query} -db {db} '
        f'-xdrop_gap 150 -evalue 1e-15 -perc_identity {cut_id} '
        f'-dust no -outfmt 6 '
        f'-max_target_seqs 1 -max_hsps 1 '
        f'-out {output} -num_threads {threads}'
    )
    print(f"Running BLAST command: {blast_cmd}")
    try:
        check_call(blast_cmd, shell=True)
    except CalledProcessError as e:
        print(f"Error executing command: {blast_cmd}", file=sys.stderr)
        sys.exit(e.returncode)

def parser_blast(blast_output, len_cut=100):
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

def main():
    # Parse command-line arguments
    args = parse_arguments()
    
    query_genome_path = args.genome1
    subject_genome_path = args.genome2
    tmp_dir = args.tmpdir

    # Create output directory if it doesn't exist
    os.makedirs(tmp_dir, exist_ok=True)

    # Split query genome and write segments to file
    query_chunks = os.path.join(tmp_dir, os.path.basename(query_genome_path)+'.chunks')
    split_genome(query_genome_path, query_chunks, 
                       chop_len=args.win, min_len=args.len, step=args.step)
    
    # Create BLAST database for the subject genome  
    subject_db = run_makeblastdb_nucl(subject_genome_path, tmp_dir)

    # Run BLAST+ alignment
    blast_output = os.path.join(tmp_dir, f'{os.path.basename(query_genome_path)}_vs_{os.path.basename(subject_genome_path)}.blast')
    run_blastn(query_chunks, subject_db, blast_output, args.threads, args.id)

    # Calculate ANI from BLAST output
    ani = parser_blast(blast_output, args.aln)

    # write ANI value
    with open(args.output, 'w') as afo:
        afo.write(f'{os.path.basename(query_genome_path)}\t{os.path.basename(subject_genome_path)}\t{ani:.2f}\n')
    

if __name__ == '__main__':
    main()
