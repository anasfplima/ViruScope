# ViruScope (CLI version) - Short Documentation

## Overview
ViruScope is a Python command line tool for viral primer retrieval from literature, in-silico primer generation, scoring, annotation and combinations.

## Features

- **AROLit**: Literature mining of PubMed/PDFs for published primers.  
- **ViruScope Core**: Automated in silico primer generation, validation, and scoring.
- **Primer Combinations**: Automated generation of Forward-Reverse primer pairs based on customizable parameters
- **Primer Scoring**:  
  - GC content  
  - Melting temperature (Tm)  
  - Conservation Score (PPI / PPI3)  
  - Structural stability (No-Fold%) via primer3  
- **Database Integration**: Store and query primers with SQLite.  
- **Visualization**: Primer overlap (Venn diagrams), distributions, and locus annotations.  

---

## Installing
Clone and install dependencies:

```bash
git clone https://github.com/anasfplima/ViruScope.git
cd ViruScope
python -m venv .venv
.venv\Scripts\activate   # Windows
# or source .venv/bin/activate   # Mac/Linux

pip install -r requirements.txt
```

## External tools
- NCBI BLAST+ (makeblastdb, blastn) must be installed and on PATH for BLAST-related methods.
- If using PubMed fetch via WSL, ensure WSL has Entrez-Direct (esearch/efetch) installed.

## Important notes
- Paths with special characters may break BLAST command-line calls — avoid underscores and spaces as well for safety.
- Suitable for larger datasets. If that is not the case, consider downloading the GUI version for ease of use.

## Usage

The main functionalities can be accessed through the classes defined in `sourcecode.py`. Below are some short usage examples:

### Example 1: Fetching MEDLINE Records and scrape PDFs with AROLit

```python
from sourcecode import arolit

# Initialize the arolit class
arolit_instance = arolit()

# Fetch MEDLINE records for a specific query and save as .nbib file
arolit_instance.fetch_pubmed_medline("COVID-19", "output_file.nbib")

# After retrieving the PDFs through a reference manager (e.g.: Zotero), scrape articles for oligonucleotide data
arolit_instance.oligos_to_csv_mult("examples/papers/", "literature_primers.csv")
```

### Example 2: Primer Search

```python
from sourcecode import viruscope

# Initialize the viruscope class
viruscope_instance = viruscope()

# Load primers from a CSV file
primers = viruscope_instance.load_primers_from_csv("sample_primers.csv")

# Generate primers from a genome file
viruscope_instance.generate_primers("sample_genome.fasta", "SampleVirus", "output_directory", "primers.fasta", 18, 25)
```

### Example 3: Running BLAST

```python
# Make BLAST database
viruscope_instance.make_blastdb_ref("blastdirectory", "reference_genome.fasta")

# Run BLAST against the reference genome
viruscope_instance.run_blast("primers.fasta", "blastdb_path", "output_directory")
```
