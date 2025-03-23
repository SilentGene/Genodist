# 🧬 GenoDist

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://semver.org)

> 🔬 A Snakemake-based bioinformatics tool for efficiently calculating genomic distances among multiple genomes using diverse metrics, including ANI, AAI, 16S rRNA genes, housekeeping genes, and pan-genome information.

## 📝 STILL IN PROGRESS
- [x] ANI
- [x] AAI
- [ ] 16S rRNA
- [ ] Housekeeping genes - GTDB
- [ ] Housekeeping genes - Conserved COGs across the tree of life
- [ ] Pan-genome

## 📊 Overview

GenoDist is a versatile bioinformatics software designed to analyze genomic relationships by calculating distances between multiple genomes through various established methods, including ANI, AAI, 16S rRNA phylogeny, housekeeping gene trees, and pan-genome trees. Built with Snakemake, GenoDist offers efficient parallel execution, comprehensive logging, and robust checkpointing, allowing analyses to resume seamlessly from interruption points. This makes GenoDist ideal for high-throughput genomic studies, providing both computational efficiency and ease-of-use.

## 🚀 Features

- **Multiple metrics** 🧮 - Supports ANI, AAI, 16S tree, house-keeping gene tree, pan-genome tree
- **Format Flexibility** 📁 - Works with genomic contigs (DNA) and proteomes (proteins)
- **Visualization Tools** 📈 - Built-in plotting capabilities for distance matrices
- **Parallel Processing** 🔄 - Utilizes available cores for maximum performance
- **Continue from Stopped** 🛑 - Automatically resumes calculations if interrupted

## 🛠️ Installation

### From Source
```bash
git clone https://github.com/yourusername/genodist.git
cd genodist
pip install .
```

## 📖 Usage

```bash
# AAI calculation
genodist aai -p proteins_folder -o output_folder -t 8 -e faa # from protein sequences
genodist aai -g genome_folder -o output_folder -t 8 -e fa  # from genome sequences

# ANI calculation
genodist ani -i genome_folder -o output_folder -t 8 -e fa  # from genome sequences
```

> [!NOTE]
> If no threading option (-t/--threads) is provided, GenoDist will automatically use all available CPU cores.

## 📄 License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.

## 🧬 Contributors

This tool was forged in the fires of computational biology by the most elite squad of AI sidekicks the silicon world has to offer:

- 🤖 **ChatGPT**: _Chief Overthinker & Occasional Genius_
- 🧙 **Claude**: _Resident Sage & Ethical Compass_
- 🛸 **Grok**: _Intergalactic Knowledge Broker & Meme Enthusiast_
- 🚀 **DeepSeek**: _The Algorithm Whisperer & Bug Hunter Extraordinaire_

_What a time to be alive!_

## 💻 Tech Stack

<div align="center">
<img src="https://img.shields.io/badge/Python3-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python3" />
<img src="https://img.shields.io/badge/Conda-44A833?style=for-the-badge&logo=anaconda&logoColor=white" alt="Conda" />
<img src="https://img.shields.io/badge/Snakemake-145374?style=for-the-badge&logo=snake&logoColor=white" alt="Snakemake" />
<img src="https://img.shields.io/badge/Pandas-150458?style=for-the-badge&logo=pandas&logoColor=white" alt="Pandas" />
</div>

---

<div align="center">
Dr Heyu Lin 🧑‍💻 2025
</div>
