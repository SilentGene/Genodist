import os
import glob

# Load config
configfile: "config.yaml"

# Define parameters
input_dir = config["genome_dna_dir"]
file_ext = config["file_ext"]
output_dir = config["output_dir"]

# Create output directory if not exists
os.makedirs(output_dir, exist_ok=True)

# Define input files
(GENOME,) = glob_wildcards(f"{input_dir}/{{GENOME}}.{file_ext}")  # f-string is used, so the wildcard is needed to be double-braced

# print(GENOME_PAIRS)


# Rule all
rule all:
    input:
        expand(os.path.join(output_dir, "prodigal_output", "{genome}.faa"), genome=GENOME)
        

# Rule to run Prodigal
rule run_prodigal:
    input:
        os.path.join(input_dir, "{genome}" + "." + file_ext)
    output:
        os.path.join(output_dir, "prodigal_output", "{genome}.faa")
    threads: 2
    params:
        genome="{genome}"
    conda: "../envs/prodigal.yaml"
    log: os.path.join(output_dir, "log", "{genome}_prodigal.log")
    shell:
        "prodigal -i {input} -o {params.genome}.gff -d {params.genome}.ffn -a {output} -q > {log} 2>&1"