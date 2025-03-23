import os
import glob
from scripts.utils import parse_blastp_output
from scripts.utils import find_bbh
from scripts.utils import calculate_aai

# Load config
configfile: "config.yaml"

# Define parameters
genome_dir = config["genome_dna_dir"]
faa_dir = config["genome_faa_dir"]
file_ext = config["file_ext"]
output_dir = config["output_dir"]
program = config.get("aai_program", "blast") # blast or diamond
tmp_dir = os.path.join(output_dir, "tmp")
dir_blast_output = os.path.join(output_dir, "blast_output")
dir_aais = os.path.join(output_dir, "AAIs")
dir_renamed = os.path.join(output_dir, "renamed")
dir_log = os.path.join(output_dir, "log")

###==== START: Define parameters according to the input type  ====###
# if genome_dir is defined, a protein_dir will be created
if genome_dir:
    faa_dir = os.path.join(output_dir, "protiens")
    os.makedirs(faa_dir, exist_ok=True)

# Create output directories
for d in [output_dir, tmp_dir, dir_blast_output, dir_aais, dir_renamed, dir_log]:
    os.makedirs(d, exist_ok=True)

# Define input files
if genome_dir:
    (GENOME,) = glob_wildcards(f"{genome_dir}/{{GENOME}}.{file_ext}")
elif faa_dir:
    (GENOME,) = glob_wildcards(f"{faa_dir}/{{GENOME}}.{file_ext}")
else:
    raise ValueError("Either genome_dir or faa_dir must be defined")

# Define suffix for different types of input
if genome_dir:
    genome_suffix = file_ext
    faa_suffix = "faa"
elif faa_dir:
    faa_suffix = file_ext

###==== END: Define parameters according to the input type  ====###

GENOME_PAIRS = [(g1, g2) for i, g1 in enumerate(GENOME) for g2 in GENOME if g1 != g2]

# Rule all
rule all:
    input:
        os.path.join(output_dir, "AAI_heatmap.pdf")

# run prodigal only when genome_dir is defined
rule run_prodigal:
    input:
        os.path.join(genome_dir, "{genome}" + "." + genome_suffix) if genome_dir else ""
    output:
        os.path.join(faa_dir, "{genome}" + "." + faa_suffix)
    params:
        genome="{genome}",
        gff_out=os.path.join(faa_dir, "{genome}" + ".gff"),
        ffn_out=os.path.join(faa_dir, "{genome}" + ".ffn")
    threads: 1
    log: os.path.join(dir_log, "{genome}" + "_prodigal.log")
    conda: "../envs/prodigal.yaml"
    shell:
        "prodigal -i {input} -o {params.gff_out} -d {params.ffn_out} -a {output} -q > {log} 2>&1"

# Rule for renaming FASTA headers
rule rename_headers:
    input:
        os.path.join(faa_dir, "{genome}" + faa_suffix)
    output:
        os.path.join(dir_renamed, "{genome}" + ".renamed.faa")
    log: os.path.join(dir_log, "rename_headers_" + "{genome}" + ".log")
    shell:
        """
        python -c 'from scripts.utils import rename_fasta_headers; \
        rename_fasta_headers("{input}", "{output}")' \
        > {log} 2>&1
        """

# Rule for creating BLAST protein database
rule blast_makedb_prot:
    input:
        os.path.join(dir_renamed, "{genome}" + ".renamed.faa")
    output:
        multiext(os.path.join(dir_renamed, "{genome}"),
            ".pdb", ".phr", ".pin", ".pot",
            ".psq", ".ptf", ".pto"
        )
    params:
        db_path=os.path.join(dir_renamed, "{genome}")
    conda: "../envs/ani_aai.yaml"
    log: os.path.join(dir_log, "makeblastdb_prot_" + "{genome}" + ".log")
    shell:
        "makeblastdb -dbtype prot -in {input} -out {params.db_path} > {log} 2>&1"

# Rule for running BLASTP from genome1 to genome2
rule blastp:
    input:
        query=os.path.join(dir_renamed, "{genome1}" + ".renamed.faa"),
        db=multiext(os.path.join(dir_renamed, "{genome2}"),
            ".pdb", ".phr", ".pin", ".pot",
            ".psq", ".ptf", ".pto"
        )
    output:
        os.path.join(dir_blast_output, "{genome1}_vs_{genome2}.blastp.tsv")
    params:
        db=os.path.join(dir_renamed, "{genome2}"),
        evalue=1e-5
    conda: "../envs/ani_aai.yaml"
    log: os.path.join(dir_log, "{genome1}_vs_{genome2}.blastp.log")
    threads: 4
    shell:
        """
        blastp -query {input.query} -db {params.db} \
            -evalue {params.evalue} \
            -outfmt '6 qseqid sseqid pident length qlen slen' \
            -out {output} -num_threads {threads} -max_target_seqs 1 > {log} 2>&1
        """

# Rule for calculating AAI
rule calculate_aai:
    input:
        blast1=os.path.join(dir_blast_output, "{genome1}_vs_{genome2}.blastp.tsv"),
        blast2=os.path.join(dir_blast_output, "{genome2}_vs_{genome1}.blastp.tsv")
    output:
        os.path.join(dir_aais, "{genome1}_vs_{genome2}.aai.tsv")
    params:
        genome1="{genome1}",
        genome2="{genome2}",
        identity=30,
        coverage=70
    log: os.path.join(dir_log, "{genome1}_vs_{genome2}_aai.log")
    run:
        import sys
        
        # Redirect stdout and stderr to log file
        with open(log[0], "w") as log_file:
            sys.stdout = log_file
            sys.stderr = log_file
            
            try:
                # Parse BLAST outputs
                bbh1 = parse_blastp_output(input.blast1, params.identity, params.coverage)
                bbh2 = parse_blastp_output(input.blast2, params.identity, params.coverage)
                
                # Find bidirectional best hits
                bbh = find_bbh(bbh1, bbh2)
                
                # Count total genes in each genome
                gene_count1 = len(set(h['query'] for h in bbh1))
                gene_count2 = len(set(h['query'] for h in bbh2))
                
                # Calculate AAI
                results = calculate_aai(bbh, min(gene_count1, gene_count2))
                
                # Write AAI results to output file
                with open(output[0], 'w') as f:
                    f.write(f'{params.genome1}\t{params.genome2}\t'
                           f'{results["aai"]:.2f}\t{results["std_dev"]:.2f}\t'
                           f'{results["hit_count"]}\t{results["query_count"]}\t'
                           f'{results["hit_percentage"]:.2f}\n')
            except Exception as e:
                print(f"Error processing AAI: {str(e)}")
                sys.exit(1)

# Rule to aggregate results
rule aggregate_results:
    input:
        expand(os.path.join(dir_aais, "{genome1}_vs_{genome2}.aai.tsv"), 
            zip, genome1=[genome1 for genome1, genome2 in GENOME_PAIRS], 
                 genome2=[genome2 for genome1, genome2 in GENOME_PAIRS])
    output:
        os.path.join(output_dir, "AAI_results.tsv")
    shell:
        """
        echo -e "Genome1\tGenome2\tAAI\tStdev\tBidirectionalHits\tQuery_Count\tHit_Percentage" > {output}
        cat {input} >> {output}
        """

# Rule to draw heatmap
rule draw_heatmap:
    input:
        os.path.join(output_dir, "AAI_results.tsv")
    output:
        os.path.join(output_dir, "AAI_heatmap.pdf")
    conda: "../envs/viz.yaml"
    log: os.path.join(dir_log, "AAI_heatmap.log")
    shell:
        "python scripts/draw_heatmap.py -i {input} -o {output} > {log} 2>&1"
