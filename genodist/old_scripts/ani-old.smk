import os
import glob

# Load config
configfile: "config.yaml"

# Define parameters
input_dir = config["genome_dna_dir"]
file_ext = config["file_ext"]
output_dir = config["output_dir"]
program = config["program"]
mode = config["mode"]

# Create output directory if not exists
os.makedirs(output_dir, exist_ok=True)

# Define input files
(GENOME,) = glob_wildcards(f"{input_dir}/{{GENOME}}.{file_ext}")  # f-string is used, so the wildcard is needed to be double-braced
GENOME_PAIRS = [(g1, g2) for i, g1 in enumerate(GENOME) for g2 in GENOME[i+1:]]

# print(GENOME_PAIRS)


# Rule all
rule all:
    input:
        os.path.join(output_dir, "ani_heatmap.pdf")

# Rule for calculating ANI
rule calculate_ani:
    input:
        g1=os.path.join(input_dir, "{g1}" + "." + file_ext),
        g2=os.path.join(input_dir, "{g2}" + "." + file_ext)
    output:
        os.path.join(output_dir, "blast_output" ,"{g1}_vs_{g2}_ani.tsv")
    threads: 2
    params:
        program=program,
        tmp=os.path.join(output_dir, "blast_tmp")
    conda: "../envs/ani.yaml"
    log: os.path.join(output_dir, "log", "{g1}_vs_{g2}_ani.log")
    shell:
        "python scripts/calculate_ani.py -1 {input.g1} -2 {input.g2} -tmp {params.tmp} -o {output} -n {threads} -p {params.program} > {log} 2>&1"

# Rule to aggregate results
rule aggregate_results:
    input:
        expand(os.path.join(output_dir, "blast_output", "{g1}_vs_{g2}_ani.tsv"), zip, g1=[g1 for g1, g2 in GENOME_PAIRS], g2=[g2 for g1, g2 in GENOME_PAIRS])
    output:
        os.path.join(output_dir, "all_ani_results.tsv")
    shell:
        "cat {input} > {output}"

# Rule to draw heatmap
rule draw_heatmap:
    input:
        os.path.join(output_dir, "all_ani_results.tsv")
    output:
        os.path.join(output_dir, "ani_heatmap.pdf"),
    conda: "../envs/viz.yaml"
    log: os.path.join(output_dir, "log", "ani_heatmap.log")
    shell:
        "python scripts/draw_heatmap.py -i {input} -o {output} > {log} 2>&1"
