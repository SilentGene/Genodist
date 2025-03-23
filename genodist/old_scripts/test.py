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


split_genome("../../genodist_test/14AMGs+MnitroUQ/Mnitro_retentate.fa", "test.chunk.txt")
