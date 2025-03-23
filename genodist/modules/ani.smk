import os
import glob
from scripts.utils import split_genome
from scripts.utils import cal_ani
from scripts.utils import write_value

# Load config
configfile: "config.yaml"

# Define parameters
input_dir = config["genome_dna_dir"]
file_ext = config["file_ext"]
output_dir = config["output_dir"]
program = config["program"]
mode = config["mode"]
dir_blast_db = os.path.join(output_dir, "blast_db")
dir_chunked = os.path.join(output_dir, "chunked")
dir_blast_output = os.path.join(output_dir, "blast_output")
dir_anis = os.path.join(output_dir, "ANIs")
dir_log = os.path.join(output_dir, "log")

# Create output directory if not exists
os.makedirs(output_dir, exist_ok=True)

# Define input files
(GENOME,) = glob_wildcards(f"{input_dir}/{{GENOME}}.{file_ext}")  # f-string is used, so the wildcard is needed to be double-braced
GENOME_PAIRS = [(g1, g2) for i, g1 in enumerate(GENOME) for g2 in GENOME if g1 != g2]

# print(GENOME_PAIRS)


# Rule all
rule all:
    input:
        os.path.join(output_dir, "ANI_heatmap.pdf")


# Rule for blast+ makeblastdb
rule blast_makedb_nuc:
    input:
        os.path.join(input_dir, "{genome}" + "." + file_ext)
    output:
        multiext(os.path.join(dir_blast_db, "{genome}"),
            ".ndb", ".nhr", ".nin", ".not",
            ".nsq", ".ntf", ".nto"
        )
    params:
        db_path=os.path.join(dir_blast_db, "{genome}")
    conda: "../envs/ani_aai.yaml"
    log: os.path.join(dir_log, "makeblastdb_" + "{genome}" + ".log")
    shell:
        "makeblastdb -dbtype nucl -in {input} -out {params.db_path} -blastdb_version 5 -parse_seqids > {log} 2>&1"

# Rule for chunking the input files
rule chunk_input:
    input:
        os.path.join(input_dir, "{genome}" + "." + file_ext)
    output:
        os.path.join(dir_chunked, "{genome}" + ".chunks")
    params:
        chop_len=1000,
        min_len=700,
        step=200
    run:
        split_genome(input[0], output[0], params.chop_len, params.min_len, params.step)

# Rule for running blastn
rule blastn:
    input:
        query=os.path.join(dir_chunked, "{genome1}" + ".chunks"),
        db=multiext(os.path.join(dir_blast_db, "{genome2}"),
            ".ndb", ".nhr", ".nin", ".not",
            ".nsq", ".ntf", ".nto"
        )
    output:
        os.path.join(dir_blast_output, "{genome1}_vs_{genome2}.blastn.tsv"),
    params:
        db=os.path.join(dir_blast_db, "{genome2}"),
        cut_id=70,
        evalue=1e-15
    conda: "../envs/ani_aai.yaml"
    log: os.path.join(dir_log, "{genome1}_vs_{genome2}.blastn.log")
    threads: 2
    shell:
        """
        blastn -task blastn -query {input.query} -db {params.db} \
            -xdrop_gap 150 -evalue {params.evalue} -perc_identity {params.cut_id} \
            -dust no -outfmt 6 \
            -max_target_seqs 1 -max_hsps 1 \
            -out {output} -num_threads {threads} > {log} 2>&1
        """

# Rule for calculating ANI
rule calculate_ani:
    input:
        os.path.join(dir_blast_output, "{genome1}_vs_{genome2}.blastn.tsv")
    output:
        os.path.join(dir_anis, "{genome1}_vs_{genome2}.ani.tsv")
    params:
        len_cut=100,
        genome1="{genome1}",
        genome2="{genome2}"
    log: os.path.join(output_dir, "log", "{genome1}_vs_{genome2}_ani.log")
    run:
        ani = cal_ani(input[0], params.len_cut)
        write_value(output[0], params.genome1, params.genome2, ani)

# Rule to aggregate results
rule aggregate_results:
    input:
        expand(os.path.join(dir_anis, "{genome1}_vs_{genome2}.ani.tsv"), 
            zip, genome1=[genome1 for genome1, genome2 in GENOME_PAIRS], 
                 genome2=[genome2 for genome1, genome2 in GENOME_PAIRS])
    output:
        os.path.join(output_dir, "ANI_results.tsv")
    shell:
        """
        echo -e "genome1\tgenome2\tANI" > {output}
        cat {input} >> {output}
        """

# Rule to draw heatmap
rule draw_heatmap:
    input:
        os.path.join(output_dir, "ANI_results.tsv")
    output:
        os.path.join(output_dir, "ANI_heatmap.pdf"),
    conda: "../envs/viz.yaml"
    log: os.path.join(dir_log, "ANI_heatmap.log")
    shell:
        "python scripts/draw_heatmap.py -i {input} -o {output} > {log} 2>&1"
