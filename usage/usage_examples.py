"""
Quickstart example for ViruScope

This script shows how to:
1. Load example data
2. Run some AROLit and iSOP key functions
3. Export results

You can find the outputs of this workflow inside the 'output' folder.
"""
#%%
import viruscopeCLI # if you want the extra utilities
#or for individual imports:
from viruscopeCLI import viruscope
from viruscopeCLI import arolit

# Initialize tool
vs = viruscope()
al = arolit()

#%%

# IMPORTANT!!! Due to BLAST+ and SeqIO handling relative paths inconsistently,
# run the functions with full paths instead to avoid raising errors.

# Fetch articles related to a virus
al.fetch_pubmed_medline(
    query="'Ebola' AND 'PCR'",
    output_file=r"\examples\output\sample_articles.nbib",
)

# You can upload the file to a reference manager like Zotero to get the PDF files automatically and export them to a folder.
# NOTE:Zotero exports PDFs to individual folders, to pull them all from their folders to the main folder in order to use the function:
# 1. Open PowerShell in the destination directory (make sure the subfolder is in this directory)
# 2. Run:
#   Get-ChildItem -Recurse -Filter *.pdf | ForEach-Object {
#       Move-Item $_.FullName -Destination . -Force
#   }

#%%
# Extract sequences from articles
al.oligos_to_csv_mult(
    folder_path=r"examples\sample_data\sample_papers",
    output_csv=r"examples\output\sample_oligos.csv",
)

#%%
# Export to SQL database
al.csv_to_sql(
    csv_path=r"examples\output\sample_oligos.csv",
    database=r"examples\output\sample_oligos.db",
)

#%%
# Make BLAST database (use full path because BLAST+ is finnicky with relative paths)
vs.make_blastdb_alig(
    save_to=r"\examples\output\blastdb",
    fasta_file=r"\examples\sample_data\reference_aligned.fasta"
)

# Generate primers with BLAST validation for the alignment simultaneously
vs.generate_primers(
    genome_file=r"\examples\sample_data\reference_genome.fasta",
    name="Demo",
    save_to=r"\examples\output",
    range_start=17,
    range_end=18,
    blast=True,
    db_path=r'\examples\output\blastdb\blastalign'
)

# Parse primers - GC, Tm, PPI, PPI3, Conservation Scores
vs.parse_silico_primers(
    align_file=r"\examples\sample_data\reference_aligned.fasta",
    save_to=r"\examples\output"
)

# Extract loci and map locations of primers - Using a GenBank file option
loci_list = vs.extract_loci_from_genbank(
    genbank_file=r"\examples\sample_data\reference.gb"
)

mapped_loci =vs.map_ref_coords_to_alignment(
    ref_seq=r"\examples\sample_data\reference_genome.fasta",
    aligned_ref_seq=r"\examples\sample_data\reference_aligned.fasta",
    loci=loci_list
)

vs.annotate_primers_with_locus(
    dict_name='filtered_silico_primers',
    loci=mapped_loci
)

# Calculate combinations
vs.calculate_primer_combinations(
    dict_name='filtered_silico_primers',
    conservation_score_threshold=99
)

# Calculate fold scores
vs.fold_scores_for_combinations(
    dict_name='primer_combinations',
)

# Export results to CSV
vs.save_to_csv(
    data='primer_combinations',
    path=r'\examples\output\combinations.csv'
)