#!/usr/bin/env python3

# ~~~[   Genodist - Genomic Distance Analysis Tool ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    ||  Copyright (C) 2025 Heyu Lin     ||
# ~~~||~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~||~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#    ||  This program is free software: you can redistribute it and/or modify ||
#    ||  it under the terms of the GNU General Public License as published by ||
#    ||  the Free Software Foundation, either version 3 of the License, or    ||
#    ||  (at your option) any later version.                                  ||
#    ||                                                                       ||
#    ||  This program is distributed in the hope that it will be useful,      ||
#    ||  but WITHOUT ANY WARRANTY; without even the implied warranty of       ||
#    ||  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        ||
#    ||  GNU General Public License for more details.                         ||
#    ||                                                                       ||
#    ||  You should have received a copy of the GNU General Public License    ||
#    ||  along with this program. If not, see <https://www.gnu.org/licenses/>.||
# ~~~||~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~||~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

__author__ = "Heyu Lin"
__copyright__ = "Copyright 20225"
__credits__ = ["Heyu Lin"]
__license__ = "GPL3"
__maintainer__ = "Heyu Lin"
__email__ = "heyu.lin AT qut.edu.au"
__status__ = "Development"

"""
A command-line interface for running the genodist Snakemake pipeline
"""

import os
import sys
import argparse
import subprocess
import logging
from ruamel.yaml import YAML
import multiprocessing
import time  # Import time module for execution timing

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def setup_logging(verbosity):
    """Configure logging based on verbosity level"""
    log_levels = {
        0: logging.WARNING,  # Default
        1: logging.INFO,     # -v
        2: logging.DEBUG     # -vv
    }
    
    # Get the appropriate level or default to DEBUG for high verbosity values
    level = log_levels.get(min(verbosity, 2), logging.DEBUG)
    
    # Configure logging
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Log the verbosity setting
    logging.debug(f"Logging level set to: {logging.getLevelName(level)}")

def create_config(mode, genome_dir, protein_dir, output_dir, program="blast", file_ext="fa"):
    """Create the configuration file in the output directory with the provided parameters"""
    try:
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        logging.debug(f"Ensured output directory exists: {output_dir}")
        
        # Define config file path in the output directory
        config_file = os.path.join(output_dir, "config.yaml")
        logging.debug(f"Config file will be created at: {config_file}")
        
        yaml = YAML()
        yaml.preserve_quotes = True
        yaml.indent(mapping=2, sequence=4, offset=2)
        
        # Create config dictionary
        config = {
            "genome_dna_dir": genome_dir,
            "genome_faa_dir": protein_dir,
            "file_ext": file_ext,
            "output_dir": output_dir,
            "program": program,
            "mode": mode
        }
        
        
        # Create a CommentedMap for proper comment handling
        from ruamel.yaml.comments import CommentedMap
        cm = CommentedMap(config)
        
        # Add comments above each item with empty lines between them
        cm.yaml_set_start_comment("Configuration file for genodist pipeline\n")
        cm.yaml_set_comment_before_after_key("genome_dna_dir", before="\nInput directory containing genome fna files")
        cm.yaml_set_comment_before_after_key("genome_faa_dir", before="\nInput directory containing genome faa files")
        cm.yaml_set_comment_before_after_key("file_ext", before="\nFile extension for genome files (e.g., fna, fasta)")
        cm.yaml_set_comment_before_after_key("output_dir", before="\nOutput directory for results")
        cm.yaml_set_comment_before_after_key("program", before="\nProgram for ANI and AAI calculation: \"blast\" or \"diamond\"")
        cm.yaml_set_comment_before_after_key("mode", before="\nMode to use: \"ani\" or \"aai\"")
        
        # Write the configuration to file
        with open(config_file, 'w') as f:
            yaml.dump(cm, f)
            
        logging.info(f"Configuration file created: {config_file}")
        return config_file
    except Exception as e:
        logging.error(f"Error creating configuration file: {e}")
        return False

def run_snakemake(threads, config_file):
    """Run the Snakemake workflow with the specified config file"""
    try:
        cmd = [
            "snakemake",
            "--cores", str(threads),
            "--sdm",
            "conda",
            "--configfile", config_file
        ]
        
        # Change to the script directory
        os.chdir(SCRIPT_DIR)
        logging.debug(f"Changed working directory to: {SCRIPT_DIR}")
        
        logging.debug(f"Running Snakemake command: {' '.join(cmd)}")
        
        # Run Snakemake
        logging.info("Starting Snakemake workflow...")
        result = subprocess.run(cmd, check=True)
        if result.returncode == 0:
            logging.info("Snakemake workflow completed successfully")
            return True
        else:
            logging.warning(f"Snakemake returned non-zero exit code: {result.returncode}")
            return False
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running Snakemake: {e}")
        return False
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        return False

def main():
    # Track the start time
    start_time = time.time()
    
    # Create the main parser
    parser = argparse.ArgumentParser(
        description="Genodist - Genomic Distance Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Add common arguments
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=multiprocessing.cpu_count(),
        help=f"Number of threads to use (default: all CPUs available: {multiprocessing.cpu_count()})"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="count",
        default=1,  # Default to INFO level
        help="Increase verbosity (can be used multiple times, e.g. -v, -vv)"
    )
    
    # Create subparsers
    subparsers = parser.add_subparsers(
        dest="mode",
        title="analysis modes",
        help="Analysis mode"
    )
    subparsers.required = True  # Make mode selection required
    
    # ANI subparser
    ani_parser = subparsers.add_parser(
        "ani", 
        help="ANI analysis",
        description="Calculate Average Nucleotide Identity between genomes"
    )
    
    # ANI-specific arguments
    ani_parser.add_argument(
        "-i", "--genome",
        required=True,
        help="Input directory containing genomic DNA files (FASTA format)"
    )
    
    ani_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for results"
    )
    
    ani_parser.add_argument(
        "-p", "--program",
        choices=["blast", "diamond"],
        default="blast",
        help="Program to use for sequence comparison (default: blast)"
    )
    
    ani_parser.add_argument(
        "-e", "--ext", "--suffix",
        default="fa",
        help="File extension for genome files (default: fa)"
    )
    
    # AAI subparser
    aai_parser = subparsers.add_parser(
        "aai", 
        help="AAI analysis",
        description="Calculate Average Amino Acid Identity between genomes"
    )
    
    # AAI-specific arguments
    # add_mutually_exclusive_group() to allow input of either genome DNA or protein files
    aai_group = aai_parser.add_mutually_exclusive_group(required=True)
    aai_group.add_argument(
        "-i", "-g" ,"--genome",
        help="Input directory containing genomic DNA files (FASTA format)"
    )

    aai_group.add_argument(
        "-f", "-faa", "--protein",
        help="Input directory containing protein files (FASTA format)"
    )
    
    aai_parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output directory for results"
    )
    
    aai_parser.add_argument(
        "-p", "--program",
        choices=["blast", "diamond"],
        default="blast",
        help="Program to use for sequence comparison (default: blast)"
    )
    
    aai_parser.add_argument(
        "-e", "--ext", "--suffix",
        help="File extension for protein files (default: fa when --genome is used, faa when --protein is used)"
    )
    
    # Parse arguments
    args = parser.parse_args()

    # Set the file extension based on the mode and input type
    if args.mode == "aai" and args.ext is None:
        if args.genome:
            args.ext = "fa"
        elif args.protein:
            args.ext = "faa"
    
    # Set up logging based on verbosity
    setup_logging(args.verbose)
    
    # Display the selected parameters
    logging.info(f"Running {args.mode.upper()} analysis")
    logging.info(f"Genome directory: {args.genome}")
    logging.info(f"Protein directory: {args.protein}")
    logging.info(f"Output directory: {args.output}")
    logging.info(f"Using program: {args.program}")
    logging.info(f"File extension: {args.ext}")
    logging.info(f"Using {args.threads} threads")
    
    # Create the configuration file in the output directory
    logging.info("Creating configuration...")
    config_file = create_config(args.mode, args.genome, args.protein, args.output, args.program, args.ext)
    if not config_file:
        logging.error("Failed to create configuration. Aborting.")
        return 1
    
    # Run Snakemake with the new config file
    logging.info(f"Running {args.mode.upper()} analysis with {args.program}...")
    if not run_snakemake(args.threads, config_file):
        logging.error("Analysis failed.")
        return 1
    
    logging.info(f"Analysis complete! Results are available in {args.output}")
    logging.info(f"Configuration saved to {config_file}")
    
    # Calculate and display the total execution time
    elapsed_time = time.time() - start_time
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    if hours > 0:
        time_str = f"{int(hours)}h {int(minutes)}m {seconds:.2f}s"
    elif minutes > 0:
        time_str = f"{int(minutes)}m {seconds:.2f}s"
    else:
        time_str = f"{seconds:.2f}s"
    
    logging.info(f"Total execution time: {time_str}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
