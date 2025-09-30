
#%%# Imports
from Bio import Entrez
from Bio import SeqIO
from Bio import AlignIO
import os
import subprocess
import re
import pandas as pd
import numpy as np
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv
import sqlite3
from itertools import product
from concurrent.futures import ThreadPoolExecutor
import primer3
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import subprocess
import pymupdf4llm
import sqlite3
import functools
from threading import Thread


#%% wrapper for timeout function

# takes the timeout parameter and returns the decorator
def timeout(sec):
    # takes target function and adds the extra code to it without altering it
    def deco(func):
        @functools.wraps(func)
        
        # replaces func when it's called
        def wrapper(*args, **kwargs):
            # creates list that stores the result of running the function (output or error)
            res = [Exception(f'Function {func.__name__} timeout {sec} seconds exceeded.')]
            
            # helper function that calls func
            def newFunc():
                try:
                    # replaces res[0] with the successful output if it exists
                    res[0] = func(*args, **kwargs)
                except Exception as e:
                    # returns the exception if function fails
                    res[0] = e
            
            # initiates thread that runs newFunc
            t = Thread(target=newFunc)
            # thread ends automatically if main program is exited
            t.daemon = True
            
            try:
                # starts the thread, calls newFunc
                t.start()
                # waits for the thread to finish in {timeout} seconds, otherwise stops waiting
                t.join(sec)
            except Exception as e:
                print('Error starting thread.')
                raise e
                
            # save the result    
            ret = res[0]
            
            # if the result is an Exception, raise it
            if isinstance(ret, BaseException):
                raise ret
            return ret
        return wrapper
    return deco

#%% AROLit
class arolit:

    # Initialize class
    def __init__(self):
        self.doi_dict = {}
        self.PMID_dict = {}
        self.primers_dict = {}
        pass

    def fetch_pubmed_medline(self, query: str, output_file: str):
        """
        Fetch MEDLINE records from PubMed using Entrez Direct (via WSL).

        Simulates:
            wsl bash -c "esearch -db pubmed -query 'query' | efetch -format medline"

        Args:
            query (str): PubMed search query.
            output_file (str): Path to save the MEDLINE result.
        """
        try:
            # Build full command as a string for bash
            full_cmd = f"esearch -db pubmed -query \"{query}\" | efetch -format medline"

            # Run it inside WSL bash
            result = subprocess.run(
                ["wsl", "bash", "-c", full_cmd],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )

            if result.returncode != 0:
                print("Command failed:", result.stderr.decode())
                return

            # Write output to file
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(result.stdout.decode())

            return(f"Saved results to: {output_file}")

        except Exception as e:
            return(f"Error running command in WSL: {e}")
        
            
    def fetch_doi_dict(self, nbib_file):
        """
        Parses an NBIB file to extract DOIs and PMIDs from PubMed entries.
        Args:
            nbib_file (str): Path to the NBIB file to be parsed.
        Returns:
            tuple: Two dictionaries:
                - doi_dict: Keys are DOIs (formatted as full URLs),
                - PMID_dict: Keys are PMIDs (formatted as full URLs)
        Notes:
            - Only articles with a DOI will be added to doi_dict.
            - Articles without a DOI will have their PMID added to PMID_dict.
            - Assumes each article entry in the NBIB file is separated by an empty line.
        """

        self.doi_dict = {}
        self.PMID_dict = {} # some articles don't have DOIs so we can try to save PMID instead (always has PMID bc it was fetched from PubMed db)
        pmid = None # so it exists outside of the loop for safety
        has_doi = False # initialize condition
        
        with open(nbib_file, 'r', encoding='utf-8') as file:
            for line in file:
                if line.startswith('PMID- '): 
                    pmid = line.strip().split('- ')[1]
                    has_doi = False # reset because new article
                    
                elif line.startswith('AID - ') and '[doi]' in line:
                    match = re.search(r'10\.\d{4,9}/\S+', line)
                    if match:
                        doi = match.group(0)
                        self.doi_dict.setdefault(f"https://doi.org/{doi}", [])
                        has_doi = True # mark article as having DOI
                        
                elif line.strip() == '': # spots the empty line inbetween entries
                    if pmid and not has_doi:
                        self.PMID_dict.setdefault(f"https://pubmed.ncbi.nlm.nih.gov/{pmid}", [])
                    pmid = None
                    has_doi = False # reset back to default second time for safe measures


        return self.doi_dict, self.PMID_dict
    
    @timeout(300)
    def extract_text_from_pdf(self, file_path): 
        '''Extracts text from a PDF file using pymupdf4llm, with a timeout of 5 minutes to avoid hanging on large files.
        Only use if you want to extract text from a single file, otherwise use oligos_to_csv_mult()'''
        text = pymupdf4llm.to_markdown(file_path)
        return text
    
    def oligo_search_new(self,article):
        pattern = re.compile(
            r'''(?ix)
            (?<![A-Za-z0-9])                                  # left boundary (no letter/number just before)
            (?:                                                # start sequence
            [ACGTURKSMYWBHNDV]                               # 1st base
            (?:[\s\-\u2010\u2011\u2012\u2013\u2014\u2212]*  # separators: space/tab/newline or PDF dash variants
            [ACGTURKSMYWBHNDV]                            # next base
            ){16,}                                           # => total ≥17 bases
            )
            (?![A-Za-z0-9])                                    # right boundary (no letter/number just after)
        ''')
        clean = lambda s: re.sub(r'[\s\-\u2010\u2011\u2012\u2013\u2014\u2212]', '', s).upper()

        hits = [clean(m.group(0)) for m in pattern.finditer(article)]
        # keep only unique, length-validated primers (safety net)
        seen = set()
        hits = [h for h in hits if len(h) >= 17 and not (h in seen or seen.add(h))]
        return hits

    def oligo_search_new2(self, article):
        pattern = re.compile(
            r'''(?ix)
            (?<![A-Za-z0-9])                                  # left boundary (no letter/number just before)
            (?:                                                # start sequence
                [ACGTURKSMYWBHNDV/()]                               # 1st base
                (?:
                    [\s\-\u2010\u2011\u2012\u2013\u2014\u2212]*  # separators: space/tab/newline or PDF dash variants
                    [ACGTURKSMYWBHNDV/()]*                            # next base
                )                                           # => total ≥17 bases
            )
            (?![A-Za-z0-9])                                    # right boundary (no letter/number just after)
        ''')

        clean = lambda s: re.sub(r'[\s\-\u2010\u2011\u2012\u2013\u2014\u2212]', '', s).upper()
        hits = [clean(m.group(0)) for m in pattern.finditer(article)]
        hits = [hit.replace(' ','') for hit in hits]
        hits = [hit.strip('()') for hit in hits]
        # keep only unique, length-validated primers (safety net)
        seen = set()
        hits = [h for h in hits if len(h) >= 17 and not (h in seen or seen.add(h))]
        return hits

    def oligos_to_csv_mult(self, folder_path, output_csv):
        """
        Extracts oligonucleotide sequences from all PDF files in a specified folder and writes them to a CSV file.
        For each PDF file in the given folder, this method:
            - Extracts the text content.
            - Searches for oligonucleotide sequences using two different search methods.
            - Collects unique sequences and associates them with the source file's name (used as a DOI placeholder).
        The resulting CSV file contains two columns: 'Sequence' and 'Source'.
        Args:
            folder_path (str): Path to the folder containing PDF files to process.
            output_csv (str): Path to the output CSV file where results will be saved.
        Returns:
            str: Status message indicating success or the reason for skipping a file.
        """

        files = [f for f in os.listdir(folder_path) if f.endswith('.pdf')]
        data = []
        
        for file in files:
            path = os.path.join(folder_path, file)
            print(f'Processing {file}...')
            
            doi = os.path.splitext(file)[0] # Placeholder, objetivo é que o titulo seja o DOI. Esta função divide o nome numa tupla com antes e depois da extensão daí o [0]
            
            try:
                # Extrair texto do PDF
                text = self.extract_text_from_pdf(path)

                oligos = self.oligo_search_new(text)
                
                oligos_2 = self.oligo_search_new2(text)

                for sequence in oligos:
                    data.append([sequence, doi])

                seen_sequences = {row[0] for row in data}

                for sequence in oligos_2:
                    if sequence not in seen_sequences:
                        data.append([sequence, doi])
                        seen_sequences.add(sequence)

            except Exception as e:
                return(f'Skipping {file} due to error: {e}')
        
        with open(output_csv, mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerow(['Sequence', 'Source'])
            writer.writerows(data)
        
        return 'Primers successfully parsed to CSV.'

    def csv_to_sql(self, csv_path, database):
        """
        Imports data from a CSV file into a SQLite database, resetting relevant tables and establishing relationships.
        The method performs the following steps:
        1. Connects to the specified SQLite database.
        2. Resets the 'sequences', 'references', and 'sequence_references' tables, including their autoincrement counters.
        3. Loads the CSV file into a pandas DataFrame.
        4. Inserts unique sequences from the 'Sequence' column into the 'sequences' table.
        5. Inserts unique references from the 'Source' column into the 'references' table.
        6. Creates associations between sequences and references in the 'sequence_references' table.
        7. Commits the changes and closes the database connection.
        Args:
            csv_path (str): Path to the CSV file containing the data to import.
            database (str): Path to the SQLite database file.
        Notes:
            - The CSV file must contain at least 'Sequence' and 'Source' columns.
            - The method uses 'INSERT OR IGNORE' to avoid duplicate entries.
            - Table names 'references' and 'sequence_references' are used; 'references' is quoted due to being a SQL reserved keyword.
        """

        conn = sqlite3.connect(database)
        cursor = conn.cursor()
        
        # Reset das tabelas para evitar problemas com o autoincrement
        cursor.execute("DELETE FROM sequence_references")  
        cursor.execute("DELETE FROM sequences")            
        cursor.execute("DELETE FROM \"references\"")
        
        cursor.execute("DELETE FROM sqlite_sequence WHERE name='sequences'")
        cursor.execute("DELETE FROM sqlite_sequence WHERE name='references'")

        # Carrega o ficheiro CSV em forma de dataframe, onde cada coluna do CSV é uma coluna da dataframe
        df = pd.read_csv(csv_path)

        # Insere sequências na base de dados
        for sequence in df['Sequence'].unique(): # Extrai todas as sequências únicas da coluna Sequences
            cursor.execute("INSERT OR IGNORE INTO sequences (sequence) VALUES (?)", (sequence,)) # Ignora o Insert se a sequência já existir

        # Insere referências na base de dados
        for reference in df['Source'].unique(): # Igual ao de cima
            cursor.execute("INSERT OR IGNORE INTO \"references\" (reference) VALUES (?)", (reference,)) # Usar \"references\" pois é uma keyword reservada no SQL

        # Junta sequências com as respetivas referências
        for _, row in df.iterrows(): # Itera sobre cada linha. Como dá return a um par Index-Conteúdo, ignoramos o index com '_'
            sequence = row['Sequence']
            reference = row['Source']
            cursor.execute("""
                INSERT OR IGNORE INTO sequence_references (sequence_id, reference_id)
                VALUES (
                    (SELECT id FROM sequences WHERE sequence = ?),
                    (SELECT id FROM "references" WHERE reference = ?)
                )
            """, (sequence, reference))

        # Commit changes and close connection
        conn.commit()
        conn.close()

    def export_for_score(self, input_db, output_csv, prefix='Primer'):
        # Creates CSV file with indexed sequences and the articles they show up in

        # Connect to the database
        ''' Database used for input should have 3 tables: Sequences (contains an
        ID column and the corresponding sequence), References (contains an ID
        column and the corresponding reference), and Sequences_references 
        (contains a column with the sequence ID from Sequences and the reference
        ID from References) '''
        conn = sqlite3.connect(input_db)
        cursor = conn.cursor()
            
        # Run the query
        query = """
        SELECT s.ID, s.Sequence, GROUP_CONCAT(r.Reference, ' | ') as Ref_list
        FROM sequences s
        LEFT JOIN sequence_references sr ON s.ID = sr.sequence_id
        LEFT JOIN "references" r ON sr.reference_id = r.ID
        GROUP BY s.ID, s.Sequence;
        """
        data = cursor.execute(query).fetchall()
            
        # Close connection
        conn.close()
            
        # Convert to DataFrame and format Sequence ID
        df = pd.DataFrame(data, columns=["Sequence ID", "Sequence", "References"])
        df["Sequence ID"] = df["Sequence ID"].apply(lambda x: f"{prefix}{x}")
            
        # Convert DOIs into full links
        def format_doi(ref_string):
            if ref_string is None:
                return ""
            doi_list = ref_string.split(" | ")
            # Change this as needed depending on how your DOI was saved as
            formatted_dois = [f"https://doi.org/{doi[0:7]}/{doi[7:]}" for doi in doi_list]
            return " | ".join(formatted_dois)
            
        # Apply formatting to the References column
        df["References"] = df["References"].apply(format_doi)
            
        # Calculate primer length
        df.insert(2, "Primer length", [len(seq) for seq in df["Sequence"].tolist()])
            
        # Save to CSV
        df.to_csv(output_csv, index=False)
            
        # Save to primers_dict
        self.primers_dict = df.to_dict(orient='records')
            
        return(f"CSV file saved to {output_csv}")
    
#%% iSOP/Viruscope
class viruscope:
    
    # Initialize class
    def __init__(self):
        self.primers_dict = {}
        self.filtered_primers = {}
        self.silico_primers = []
        self.filtered_silico_primers = []
        self.IUPAC_CODES = {
            "A": ["A"], "T": ["T"], "C": ["C"], "G": ["G"],
            "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"],
            "K": ["G", "T"], "M": ["A", "C"], "B": ["C", "G", "T"],
            "D": ["A", "G", "T"], "H": ["A", "C", "T"], "V": ["A", "C", "G"],
            "N": ["A", "C", "G", "T"], "U": ["T"], "/": [""], "(": [""], ")": [""], "-": ["A", "C", "G", "T"]
        }
        self.AMBIGUITY_MAP = {
            b"A": {b"A": 1.0, b"R": 0.5, b"W": 0.5, b"M": 0.5, b"D": 0.33, b"H": 0.33, b"V": 0.33, b"N": 0.25},
            b"C": {b"C": 1.0, b"Y": 0.5, b"M": 0.5, b"S": 0.5, b"B": 0.33, b"H": 0.33, b"V": 0.33, b"N": 0.25},
            b"G": {b"G": 1.0, b"R": 0.5, b"K": 0.5, b"S": 0.5, b"B": 0.33, b"D": 0.33, b"V": 0.33, b"N": 0.25},
            b"T": {b"T": 1.0, b"Y": 0.5, b"K": 0.5, b"W": 0.5, b"B": 0.33, b"D": 0.33, b"H": 0.33, b"N": 0.25},
            b"R": {b"A": 0.5, b"G": 0.5, b"R": 1.0, b"N": 0.25},
            b"Y": {b"C": 0.5, b"T": 0.5, b"Y": 1.0, b"N": 0.25},
            b"S": {b"C": 0.5, b"G": 0.5, b"S": 1.0, b"N": 0.25},
            b"W": {b"A": 0.5, b"T": 0.5, b"W": 1.0, b"N": 0.25},
            b"K": {b"G": 0.5, b"T": 0.5, b"K": 1.0, b"N": 0.25},
            b"M": {b"A": 0.5, b"C": 0.5, b"M": 1.0, b"N": 0.25},
            b"B": {b"C": 0.33, b"G": 0.33, b"T": 0.33, b"B": 1.0, b"N": 0.25},
            b"D": {b"A": 0.33, b"G": 0.33, b"T": 0.33, b"D": 1.0, b"N": 0.25},
            b"H": {b"A": 0.33, b"C": 0.33, b"T": 0.33, b"H": 1.0, b"N": 0.25},
            b"V": {b"A": 0.33, b"C": 0.33, b"G": 0.33, b"V": 1.0, b"N": 0.25},
            b"N": {b"A": 0.25, b"C": 0.25, b"G": 0.25, b"T": 0.25, b"N": 1.0},
        }
        pass
    
    # Function for faster dict to csv saving
    def save_to_csv(self, data, path):
        """
        Saves the specified attribute data of the class to a CSV file.
        Parameters:
            data (str): The name of the attribute containing the data to be saved. The attribute is the name of the ViruScope list containing the primer information.
            path (str): The file path where the CSV file will be saved.
        Returns:
            None
        Notes:
            - If the specified attribute does not exist or is None, a message is printed and no file is saved.
        """
        
        data = getattr(self, data, None)
        if data is not None:
            df = pd.DataFrame(data)
            df.to_csv(path, index=False)
            print(f"Data saved to {path}.")
        else:
            print(f"No attribute named {data}.")
        
    # Load primers from a CSV file
    def load_primers_from_csv(self, csv_path):
        """
        Loads primer information from a CSV file and converts numeric values to appropriate types.
        Args:
            csv_path (str): The file path to the CSV file containing primer data.
        Returns:
            list[dict]: A list of dictionaries, each representing a primer with keys as column headers.
                        Numeric values are converted to int or float where possible; otherwise, they remain as strings.
        Raises:
            FileNotFoundError: If the specified CSV file does not exist.
            csv.Error: If there is an error reading the CSV file.
        """

        def convert_value(val):
            try:
                return int(val)
            except ValueError:
                try:
                    return float(val)
                except ValueError:
                    return val  # Leave as string if not numeric

        with open(csv_path, newline='') as f:
            reader = csv.DictReader(f)
            primers = []
            for row in reader:
                converted_row = {key: convert_value(value) for key, value in row.items()}
                primers.append(converted_row)
        return primers  
    
    # Creates a list of all possible primers from the reference sequence
    def generate_primers(self, genome_file, name: str, save_to, range_start: int, range_end: int, step=1, blast=False, db_path='',  output_file=''):
        '''
        
        Parameters
        ----------
        genome_file: input file with the genome the primers are gonna be created for.
        name : name to give the primers for identification (ex.: 'Ebola' > Ebola1, Ebola2).
        range_start : minimum primer length.
        range_end : maximum primer length.
        step : how much the sliding window moves forward each iteration. Default is 1.

        Returns
        -------
        None.

        '''
        os.chdir(save_to)
        
        # Clear previous runs if necessary
        self.silico_primers = []
        self.filtered_silico_primers = []
        
        # Open and store sequence
        sequence = (SeqIO.read(genome_file, 'fasta')).seq
        
        seq_len = len(sequence)
        primer_id = 0 
        
        # Loops over every possible length for the primer
        for length in range(range_start, range_end+1):
            # Sliding window of size 'length' and step = value over the sequence
            # Obtains the start and end positions of the primer
            for start in range(0, seq_len-length+1, step): # second argument prevents index out of range
                end = start + length
                
                # Forward primer
                forward = sequence[start:end]
                
                # Reverse primer
                reverse = forward.reverse_complement()
                
                # Add the entries to the dictionary of in silico primers initialized with the class
                for primer_type, seq in [("Forward", forward), ("Reverse", reverse)]:
                    primer_id += 1
                    primer_dict = {
                        'Primer ID': f"{name}Silico{primer_id}",
                        'Sequence': str(seq),
                        'Primer length': len(seq),
                        'Ref Start': start if primer_type == "Forward" else end,
                        'Ref End': end if primer_type == "Forward" else start,
                        'Primer Start': 1,
                        'Primer End': len(seq),
                        'Alignment Sequence': str(seq),  # Placeholder: actual alignment not done here
                        'Length in Alignment': len(seq),  # Same as above
                        'GC Content': '',
                        'Melting Temperature': '',
                        'Type': primer_type,
                        'PPI': 0.0,   # Placeholder
                        'PPI3': 0.0,  # Placeholder
                        'Score': 0.0  # Placeholder
                    }
                    self.silico_primers.append(primer_dict)
        
        # Save before blast            
        records = []
        
        for primer in self.silico_primers:
            seq_id = primer.get("Primer ID", "unknown_id")
            seq_str = primer.get("Sequence", "")
            if seq_str:  # Avoid empty sequences
                record = SeqRecord(Seq(seq_str), id=seq_id, description="")
                records.append(record)
        
        SeqIO.write(records, 'primersbeforeblast.fasta', "fasta")

        # BLAST aligned ref genome
        if blast:
            
            self.run_blast('primersbeforeblast.fasta', db_path, save_to)
            # Opens and filters the BLAST output file
            blast_data = pd.read_csv('primer_hits.txt', sep='\t', names=["Primer ID", "genome seqid", "Alignment Sequence", "pident", "Length in Alignment", "mismatch",
               "Primer Start", "Primer End", "Ref Start", "Ref End", "evalue"])
            # Drop partial primers
            blast_data = blast_data[blast_data['Length in Alignment'] >= 17]
            # Clean up the BLAST data
            blast_data = blast_data[['Primer ID', 'Ref Start', 'Ref End']]
            blast_data = blast_data.drop_duplicates(subset='Primer ID', keep='first')
            
            # Convert BLAST dataframe to a lookup dict
            ref_lookup = blast_data.set_index('Primer ID')[['Ref Start', 'Ref End']].to_dict('index')

            # Update the original primer list
            for primer in self.silico_primers:
                primer_id = primer.get('Primer ID')
                if primer_id in ref_lookup:
                    primer['Ref Start'] = ref_lookup[primer_id]['Ref Start']
                    primer['Ref End'] = ref_lookup[primer_id]['Ref End']
                    self.filtered_silico_primers.append(primer)
                    
            # Export the generated primers as a FASTA file
            records = []
        
            for primer in self.filtered_silico_primers:
                seq_id = primer.get("Primer ID", "unknown_id")
                seq_str = primer.get("Sequence", "")
                if seq_str:  # Avoid empty sequences
                    record = SeqRecord(Seq(seq_str), id=seq_id, description="")
                    records.append(record)
        
            SeqIO.write(records, 'primersafterblast.fasta', "fasta")
        

        return print('In Silico primer database created successfully.')
    
    # Generates a Venn diagram comparing two sets of primers
    def primer_venn_diagram(self, file1, file2, label1="File1", label2="File2"):
        """
        Generates and saves a Venn diagram illustrating the overlap between two sets of primer sequences from FASTA files.
        Args:
            file1 (str): Path to the first FASTA file containing primer sequences.
            file2 (str): Path to the second FASTA file containing primer sequences.
            label1 (str, optional): Label for the first set in the Venn diagram. Defaults to "File1".
            label2 (str, optional): Label for the second set in the Venn diagram. Defaults to "File2".
        Returns:
            None
        Side Effects:
            - Displays the Venn diagram using matplotlib.
            - Saves the Venn diagram as 'primers_venn.png' in the current working directory.
            - Prints a confirmation message upon saving the diagram.
        """

        # Read sequences (use IDs or sequences for comparison)
        primers1 = {str(record.seq) for record in SeqIO.parse(file1, "fasta")}
        primers2 = {str(record.seq) for record in SeqIO.parse(file2, "fasta")}

        # Make the Venn diagram
        venn2([primers1, primers2], set_labels=(label1, label2))
        plt.title("Primer Overlap")
        plt.show()

        # Save as PNG
        plt.savefig("primers_venn.png", dpi=300)
        plt.close()

        return print("Venn diagram saved as primers_venn.png")    
    
    # Parses the silico primers, calculates GC content, melting temperature, PPI, PPI3, and conservation scores
    def parse_silico_primers(self, align_file, save_to):
        """
        Parses and processes in silico primers by performing a series of calculations and saving results.
        This method executes the following steps:
            1. Calculates GC and melting temperature (Tm) for in silico primers.
            2. Runs PPI and PPI3 calculations concurrently using a thread pool.
            3. Computes conservation scores for the filtered primers and saves the results to a CSV file.
        Args:
            align_file (str): Path to the alignment file used for PPI and PPI3 calculations.
            save_to (str): Directory or file path where the conservation scores CSV will be saved.
        Side Effects:
            - Prints progress messages to the console.
            - Saves conservation scores to 'silico_conservation_scores.csv' in the specified location.
        Note:
            - The method assumes that 'filtered_silico_primers' is a valid intermediate result used across calculations.
            - If filtered_silico_primers is not set, use load_primers_from_csv() and set it to (instance name).filtered_silico_primers before calling this method.
            - The function was built for silico primers but can be used for any primer list as long as it has the same structure.
        """

        self.GC_MT_calculations('filtered_silico_primers')
        print('GC MT calculations done.')

        # Use ThreadPoolExecutor for independent tasks
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(self.PPI_calculations, align_file, 'filtered_silico_primers'),
                executor.submit(self.PPI3_calculations, align_file, 'filtered_silico_primers')
            ]
            # Wait for both tasks to complete
            for future in futures:
                future.result()
        print('PPI and PPI3 calculations done.')

        # Run conservation scores after PPI and PPI3 are complete
        self.conservation_scores('filtered_silico_primers', save_to, 'silico_conservation_scores.csv')
        print('Conservation scores calculations done.')

    # Filters the silico primers based on conservation score threshold
    def filter_silico_primers(self, threshold=90):
        '''
        Filters the silico primers based on the conservation score threshold.
        
        Parameters
        ----------
        threshold : float, optional
            The minimum conservation score for a primer to be considered valid.
            Default is 90.
        
        Returns
        -------
        list
            A list of filtered primers that meet the conservation score criteria.
        '''
        self.filtered_silico_primers = [primer for primer in self.silico_primers if primer['Score'] >= threshold]
        return print('Silico primers filtered successfully.')
    
          
    # Create fasta file with all the primers for the BLAST
    def primers_to_fasta(self, primers_csv, fasta_output):
        '''
        Parameters
        ----------
        primers_csv : str
            path to the primers csv file created in export_for_score().
        fasta_output : str
            name/path of the resulting fasta file.
        
        '''
        data = pd.read_csv(primers_csv)
        with open(fasta_output, 'w') as file:
            for primer, row in data.iterrows():
                file.write(f">{row['Sequence ID']}\n{row['Sequence']}\n")
            file.close()
            
        return print(f"FASTA file {fasta_output} created successfully.")
        
    # Create a fasta file with the reference genome for the BLAST
    def get_ref_from_alig(self,fasta_path, ref_id, save_to):
        '''
        Parameters
        ----------
        fasta_path : str
            path to the alignment fasta file.
        ref_id : str
            ID assigned to the reference genome by NCBI.
        save_to : str
            path to where you want to save the resulting file.

        Returns
        -------
        None.

        '''
        os.chdir(save_to)
        
        #ATTENTION: couldn't use the MAFFT without dupes bc the reference genome got removed in those FIX THIS
        ref_fasta = open('reference_genome.fasta', 'w')
        for record in SeqIO.parse(fasta_path, 'fasta'):
            if record.id == ref_id:
                ref_fasta.write('>'+record.id+'\n')
                ref_fasta.write(str(record.seq.upper())+'\n')
        print(f'File "reference_genome.fasta" successfully saved to {save_to}.')
        ref_fasta.close()
                
        
    # Create BLAST database with reference genome to prep for nucleotide BLAST
    def make_blastdb_ref(self, save_to, fasta_file):
        '''

        Parameters
        ----------
        save_to : str
            path to where you want to save the resulting file.
        fasta_file : str
            path to the reference_genome.fasta generated with get_ref_from_alig().

        Returns
        -------
        Success or failure status.

        '''
        os.chdir(save_to)
        
        # Sends prompt to command line. ATTENTION: DO NOT HAVE SPECIAL CHARACTERS IN YOUR PATHS, NCBI BLAST CANNOT PARSE THOSE
        command_line="makeblastdb -in " +fasta_file+ " -dbtype nucl -out blastref -logfile log1.txt"
        child_process = subprocess.Popen(str(command_line),stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        blast_results, err = child_process.communicate()
        print(blast_results.decode(), err.decode())
        status = child_process.wait()
        
        return print("Database created successfully.") if status == 0 else print("Failed to create database. Check the log file for error handling.")
    
    
    # Create BLAST database with the MAFFT alignment to prep for nucleotide BLAST
    def make_blastdb_alig(self, save_to, fasta_file):
        '''
        
        Parameters
        ----------
        save_to : str
            path to where you want to save the resulting file.
        fasta_file : str
            path to the MAFFT alignment file.

        Returns
        -------
        Success or failure status

        '''
        os.chdir(save_to)
        
        # Sends prompt to command line. ATTENTION: DO NOT HAVE SPECIAL CHARACTERS IN YOUR PATHS, NCBI BLAST CANNOT PARSE THOSE
        command_line="makeblastdb -in " +fasta_file+ " -dbtype nucl -out blastalign -logfile log1.txt"
        child_process = subprocess.Popen(str(command_line),stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        blast_results, err = child_process.communicate()
        print(blast_results.decode(), err.decode())
        status = child_process.wait()
        
        return print("Database created successfully.") if status == 0 else print("Failed to create database. Check the log file for error handling.")
    
    
    # Run BLAST against reference genome/alignment
    def run_blast(self, primers_fasta, db_path, save_to):
        '''
        Parameters
        ----------
        primers_fasta : str
            FASTA file with all the primers. If you don't have one,
            create it using primers_to_fasta().
        db_path : str
            path to the database you want to run the BLAST against.
        save_to : str
            path where you want to save the results file.

        '''
        os.chdir(save_to)
        # Sends prompt to command line. ATTENTION: DO NOT HAVE SPECIAL CHARACTERS IN YOUR PATHS, NCBI BLAST CANNOT PARSE THOSE
        # tried gapopen 1 wordsize 10 and gapopen2 wordsize 18
        command_line='blastn -gapopen 1 -word_size 10 -query '+primers_fasta+' -db '+db_path+' -task blastn-short -outfmt "6 qseqid sseqid qseq pident length mismatch qstart qend sstart send evalue" -evalue 0.01 -out primer_hits.txt -logfile log.txt'
        child_process = subprocess.Popen(str(command_line),stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        blast_results, err = child_process.communicate()
        print(blast_results.decode(), err.decode())
        status = child_process.wait()
        
        return print("BLAST run successfully.") if status == 0 else print("Failed to run BLAST. Check the log file for error handling.")
    
    # Run BLAST exhaustively against reference genome/alignment for silico primers
    def run_blast_exhaustive(self, primers_fasta, db_path, save_to):
        '''
        Parameters
        ----------
        primers_fasta : str
            FASTA file with all the primers. If you don't have one,
            create it using primers_to_fasta().
        db_path : str
            path to the database you want to run the BLAST against.
        save_to : str
            path where you want to save the results file.

        '''
        os.chdir(save_to)
        # Sends prompt to command line. ATTENTION: DO NOT HAVE SPECIAL CHARACTERS IN YOUR PATHS, NCBI BLAST CANNOT PARSE THOSE
        # tried gapopen 1 wordsize 10 and gapopen2 wordsize 18
        command_line='blastn -gapopen 0 -gapextend 0 -word_size 4 -query '+primers_fasta+' -db '+db_path+' -task blastn -reward 1 -penalty 1 -soft_masking false -outfmt "6 qseqid sseqid qseq pident length mismatch qstart qend sstart send evalue" -evalue 1000 -out primer_hits_exhaustive.txt -logfile log.txt'
        child_process = subprocess.Popen(str(command_line),stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        blast_results, err = child_process.communicate()
        print(blast_results.decode(), err.decode())
        status = child_process.wait()
            
        return print("BLAST run successfully.") if status == 0 else print("Failed to run BLAST. Check the log file for error handling.")
        
    # Parse the results from the BLAST
    def parse_blast(self, results_file, input_csv, output_csv, virus='Primer'):
        '''
        
        Parameters
        ----------
        results_file : path to the txt file obtained from run_blast().
        input_csv : path to the csv created in export_for_score().
        output_csv : path to where you want to save your parsed primers.
        virus : format of the primer index. The default is 'Primer'.

        '''
        # Opens and filters the BLAST output file
        blast_data = pd.read_csv(results_file, sep='\t', names=["Primer ID", "genome seqid", "Alignment Sequence", "pident", "Length in Alignment", "mismatch",
           "Primer Start", "Primer End", "Ref Start", "Ref End", "evalue"])
        blast_data = blast_data[['Primer ID', 'Ref Start', 'Ref End','Primer Start', 'Primer End', 'Alignment Sequence', 'Length in Alignment']]
        
        # Opens the file with all the extracted primers
        og_data = pd.read_csv(input_csv)
        
        # Only keeps the primers that exist in both files
        merged_data = pd.merge(blast_data, og_data, left_on='Primer ID', right_on='Sequence ID')
        
        # Gets rid of duplicate ID column
        merged_data.drop(columns=['Sequence ID'], inplace=True)
        
        # For alignment: gets rid of duplicate entries (there's one for each sequence)
        # Keeps only 1 entry for each primer+position pair
        merged_data = merged_data.drop_duplicates(subset=['Primer ID', 'Ref Start', 'Ref End'], keep='first')
        
        # New IDs
        merged_data['Primer ID'] = [virus + str(i+1) for i in range(len(merged_data))]
        
        merged_data = merged_data[['Primer ID', 'Sequence', 'Primer length','Ref Start', 'Ref End', 'Primer Start', 'Primer End','Alignment Sequence', 'Length in Alignment']]
        
        merged_data.to_csv(output_csv, index=False)
        
        # Saves to class 
        self.filtered_primers = merged_data.to_dict(orient='records')
        
        return print(f'Primers parsed to {output_csv}.') # Ask if it's to only keep the primers with no mismatches
    
    # Generate all possible sequences if a primer has degenerate bases
    def expand_degenerate(self,sequence):
        possible_bases = [self.IUPAC_CODES[base] for base in sequence]
        sequences = ["".join(p) for p in product(*possible_bases)]
        return sequences
    
    # Calculate GC content
    def calculate_gc(self,sequence):
        gc_count = sequence.count("G") + sequence.count("C")
        gc_content = (gc_count / len(sequence)) * 100
        return gc_content
    
    # Calculate melting temperatures
    def calculate_tm(self,sequence):
        gc_count = sequence.count("G") + sequence.count("C")
        total_count = sequence.count('A')+sequence.count('T')+gc_count
            
        if len(sequence) <= 13:
            tm = (sequence.count('A')+sequence.count('T'))*2 + gc_count*4
        else:
            tm = 64.9 + 41*(gc_count-16.4)/total_count
        return  tm
    
    # Handle cases based on whether or not primer has degenerated bases (scrapped but kept as an extra)
    def calculate_gc_tm_range(self, sequence):
        expanded_sequences = self.expand_degenerate(sequence)
        gc_values = [self.calculate_gc(seq) for seq in expanded_sequences]
        tm_values = [self.calculate_tm(seq) for seq in expanded_sequences]
        if min(gc_values) != max(gc_values):
            gc_value = f'{min(gc_values):.2f}% - {max(gc_values):.2f}%'
        else:
            gc_value =  f'{min(gc_values):.2f}%'
        if min(tm_values) != max(tm_values):
            tm_value = f'{min(tm_values):.2f}°C - {max(tm_values):.2f}°C'
        else:
            tm_value = f'{min(tm_values):.2f}°C'
        return gc_value, tm_value
    
    # New function that keeps the expanded sequences and calculates GC and Tm for each
    def calculate_gc_tm_expanded(self, sequence, primer):
        """
        Calculates GC content and melting temperature (Tm) for all expanded sequences of a primer.

        Args:
            sequence (str): The primer sequence, potentially containing degenerate bases.
            primer (dict): The primer metadata, including its ID and other information.

        Returns:
            list: A list of dictionaries, each containing the expanded sequence and its GC/Tm values.
        Notes:
            - The function expands degenerate bases in the primer sequence to generate all possible sequences.
            - It calculates GC content and Tm for each expanded sequence.
            - Each resulting dictionary includes the original primer metadata with updated sequence, GC content, Tm, and a unique Primer ID.
            - This function is called by GC_MT_calculations() to process a list of primers. Use this function when you want to process a single primer, otherwise use GC_MT_calculations().
        """
        expanded_sequences = self.expand_degenerate(sequence)
        expanded_primers = []
        for i, seq in enumerate(expanded_sequences):
            gc = self.calculate_gc(seq)
            tm = self.calculate_tm(seq)
            new_primer = primer.copy()
            new_primer['Sequence'] = seq
            new_primer['GC Content'] = round(gc, 2)
            new_primer['Melting Temperature'] = round(tm, 2)
            new_primer['Primer ID'] = f"{primer['Primer ID']}_var{i+1}"
            expanded_primers.append(new_primer)
        return expanded_primers

    # Appends GC and Tm information to class dictionary
    def GC_MT_calculations(self, dict_name):
        """
        Appends GC and Tm information to the specified class dictionary.

        Args:
            dict_name (str): The name of the dictionary attribute to update.
        """
        primer_list = getattr(self, dict_name, None)
        expanded_list = []
        for primer in primer_list:
            expanded = self.calculate_gc_tm_expanded(primer['Sequence'], primer)
            expanded_list.extend(expanded)
        # Replace the original list with the expanded list
        setattr(self, dict_name, expanded_list)

    # Classifies primers as Forward or Reverse
    def primer_type(self):
        """ 
        Classifies primers as Forward or Reverse based on their reference start and end positions.
        Updates the 'Type' key in each primer dictionary within the filtered_primers list. If you don't have a filtered_primers list, use load_primers_from_csv() and
        set it to (instance name).filtered_primers before calling this method.

        In silico primers are already classified when generated.
        """
        for primer in self.filtered_primers:
            if primer['Ref Start'] < primer['Ref End']:
                primer['Type'] = 'Forward'
            else:
                primer['Type'] = 'Reverse'
    
    # Calculates pairwise identity for a given column in the alignment
    def pairwise_identity(self, column, ambiguity_map):
        column_nogaps = column[column != b'-']  # Remove gaps
        if len(column_nogaps) < (0.3*len(column)):
            return 0.0
        
        unique, counts = np.unique(column_nogaps, return_counts=True)
        total_pairs = (len(column_nogaps) * (len(column_nogaps) - 1)) / 2
        
        identical_pairs = sum(
            (counts[i] * (counts[i] - 1)) / 2 for i in range(len(unique))
        )
        
        for i in range(len(unique)):
            for j in range(i + 1, len(unique)):
                b1, b2 = unique[i], unique[j]
                if b1 in ambiguity_map and b2 in ambiguity_map[b1]:
                    identical_pairs += counts[i] * counts[j] * ambiguity_map[b1][b2]
        
        return (identical_pairs / total_pairs) if total_pairs > 0 else 0.0
            
    # Calculates PPI (Percentage of Pairwise Identity) for primers in the alignment
    def PPI_calculations(self, align_file, dict_name):
        primer_list = getattr(self, dict_name, None)
        alignment = AlignIO.read(align_file, 'fasta')
        
        def process_primer(primer):
            start = primer['Ref Start']
            end = primer['Ref End']
            window = alignment[:, min(start, end):max(start, end)]
            alignment_array = np.array([list(record.seq) for record in window], dtype='S1')
            num_columns = alignment_array.shape[1]
            
            pairwise_identity_sum = 0
            valid_columns = 0
            
            for col_idx in range(num_columns):
                column = alignment_array[:, col_idx]
                if np.count_nonzero(column == b'-') < len(column) - 1:
                    valid_columns += 1
                    pairwise_identity_sum += self.pairwise_identity(column, self.AMBIGUITY_MAP)
            
            pairwise_identity_avg = (pairwise_identity_sum / valid_columns) * 100 if valid_columns > 0 else 0
            primer['PPI'] = round(pairwise_identity_avg, 2)

        with ThreadPoolExecutor() as executor:
            list(executor.map(process_primer, primer_list))

    # Calculates PPI3 (Percentage of Pairwise Identity for the last 3 bases) for primers in the alignment
    def PPI3_calculations(self, align_file, dict_name):
        primer_list = getattr(self, dict_name, None)
        alignment = AlignIO.read(align_file, 'fasta')
        
        def process_primer(primer):
            start = primer['Ref Start']
            end = primer['Ref End']
            window = alignment[:, min(start, end):max(start, end)]
            
            if primer['Type'] == 'Forward':
                start3 = primer['Ref End'] - 3
                end3 = primer['Ref End']
                window = alignment[:, start3:end3]
            else:
                alignment_rev_comp = MultipleSeqAlignment(
                    SeqRecord(seq.seq.reverse_complement(), id=seq.id)
                    for seq in window
                    )
                window = alignment_rev_comp[:, 0:3]
           
            alignment_array = np.array([list(record.seq) for record in window], dtype='S1')
            num_columns = alignment_array.shape[1]
            
            pairwise_identity_sum = 0
            valid_columns = 0
            
            for col_idx in range(num_columns):
                column = alignment_array[:, col_idx]
                if np.count_nonzero(column == b'-') < len(column) - 1:
                    valid_columns += 1
                    pairwise_identity_sum += self.pairwise_identity(column, self.AMBIGUITY_MAP)
            
            pairwise_identity_avg = (pairwise_identity_sum / valid_columns) * 100 if valid_columns > 0 else 0
            primer['PPI3'] = round(pairwise_identity_avg, 2)

        with ThreadPoolExecutor() as executor:
            list(executor.map(process_primer, primer_list))

    # Calculates conservation scores based on PPI and PPI3
    def conservation_scores(self, dict_name, save_to, csv_name):
        """
        Calculates conservation scores for a list of primers and saves the results to a CSV file.
        The method retrieves a list of primers from an attribute of the instance specified by `dict_name`.
        For each primer, it computes a 'Score' as the average of the 'PPI' and 'PPI3' values, rounded to two decimal places.
        The updated list is then saved as a CSV file in the specified directory.
        Args:
            dict_name (str): The name of the attribute containing the list of primer dictionaries.
            save_to (str): The directory path where the CSV file will be saved.
            csv_name (str): The name of the CSV file to write the results to.
        Raises:
            AttributeError: If the specified attribute does not exist.
            KeyError: If 'PPI' or 'PPI3' keys are missing in any primer dictionary.
            OSError: If changing the directory or writing the file fails.
        """
        
        primer_list = getattr(self, dict_name, None)
        for primer in primer_list:
            primer['Score'] = round(((primer['PPI']+primer['PPI3'])/2),2)
        os.chdir(save_to)
        with open(csv_name, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=primer_list[0].keys())
            writer.writeheader()
            writer.writerows(primer_list)

    # Extracts gene/locus regions from a GenBank file or NCBI accession
    def extract_loci_from_genbank(self, genbank_file=None, ncbi_id=None, email=None):
        """
        Extracts gene/locus regions from a GenBank file or NCBI accession.
        Returns a list of dicts: [{'name': gene, 'start': X, 'end': Y}, ...]
        """
        loci = []
        if genbank_file:
            record = SeqIO.read(genbank_file, "genbank")
        elif ncbi_id and email:
            Entrez.email = email
            with Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="text") as handle:
                record = SeqIO.read(handle, "genbank")
        else:
            raise ValueError("Provide either a genbank_file or both ncbi_id and email.")
        # Check for exclusive/inclusive end positions
        for feature in record.features:
            if feature.type in ["gene", "CDS"]:
                name = feature.qualifiers.get("gene", ["unnamed"])[0]
                start = int(feature.location.start)
                end = int(feature.location.end)
                loci.append({'name': name, 'start': start, 'end': end})

        # Merging loci that appear multiple times with different start/end positions
        # Based on the graph at: https://www.ncbi.nlm.nih.gov/nuccore/NC_002549.1?report=graph
        merged = {}
        for locus in loci:
            name = locus['name']
            start = locus['start']
            end = locus['end']
            if name not in merged:
                merged[name] = {'name': name, 'start': start, 'end': end}
            else:
                merged[name]['start'] = min(merged[name]['start'], start)
                merged[name]['end'] = max(merged[name]['end'], end)
        loci = list(merged.values())
        return loci
    
    # Maps reference coordinates to alignment coordinates
    def map_ref_coords_to_alignment(self, ref_seq, aligned_ref_seq, loci):
        """
        Maps loci coordinates from unaligned reference to aligned reference sequence.
        ref_seq: str, unaligned reference sequence
        aligned_ref_seq: str, aligned reference sequence (with gaps)
        loci: list of dicts [{'name': ..., 'start': ..., 'end': ...}]
        Returns: list of dicts with mapped coordinates
        """
        # Load the reference sequence (no gaps)
        ref_record = SeqIO.read(ref_seq, "fasta")
        ref_seq = str(ref_record.seq)

        # Load the aligned reference sequence (with gaps)
        aligned_record = SeqIO.read(aligned_ref_seq, "fasta")
        aligned_ref_seq = str(aligned_record.seq)

        ref_pos = 0
        align_pos = 0
        ref_to_align = {}
        while ref_pos < len(ref_seq) and align_pos < len(aligned_ref_seq):
            if aligned_ref_seq[align_pos] != '-':
                ref_to_align[ref_pos] = align_pos
                ref_pos += 1
            align_pos += 1

        mapped_loci = []
        for locus in loci:
            start = locus['start']
            end = locus['end'] - 1  # GenBank end is exclusive
            align_start = ref_to_align.get(start)
            align_end = ref_to_align.get(end)
            if align_start is not None and align_end is not None:
                mapped_loci.append({
                    'name': locus['name'],
                    'aligned_start': align_start + 1, # Python uses 0-based indexing, GenBank uses 1-based
                    'aligned_end': align_end + 1  # keep exclusive
                })
            else:
                mapped_loci.append({
                    'name': locus['name'],
                    'aligned_start': None,
                    'aligned_end': None
                })
        return mapped_loci
    
    def annotate_primers_with_locus(self, dict_name, loci):
        """
        Annotate primers with the gene/locus they bind to.

        Args:
            dict_name (str): The name of the dictionary attribute containing the primers to annotate.
            loci (list): A list of loci to annotate the primers with. You can use extract_loci_from_genbank() and map_ref_coords_to_alignment() to get this list.
        Returns:
            None: The method updates the 'Locus' key in each primer.
        """
        primer_list = getattr(self, dict_name, None)
        if primer_list is None:
            print(f"No primer list named {dict_name}")
            return

        for primer in primer_list:
            primer_start = min(primer['Ref Start'], primer['Ref End'])
            primer_end = max(primer['Ref Start'], primer['Ref End'])
            found = False
            for locus in loci:
                if primer_start >= locus['aligned_start'] and primer_end <= locus['aligned_end']:
                    primer['Locus'] = locus['name']
                    found = True
                    break
            if not found:
                primer['Locus'] = 'Intergenic/Unknown'
    
        self.filtered_silico_primers = primer_list
    
    # Calculate primer3 parameters
    def calc_homodimer(self,seq):
        return str(primer3.calc_homodimer(seq)).split(",")[2]

    def calc_hairpin(self, seq):
        return str(primer3.calc_hairpin(seq)).split(",")[2]

    def calc_heterodimer(self,row):
        return str(primer3.calc_heterodimer(row['Sequence_fwd'], row['Sequence_rev'])).split(",")[2]
    
    # Calculate primer combinations
    def calculate_primer_combinations(self, dict_name, min_amplicon_length=70, max_amplicon_length=1000, conservation_score_threshold=95, min_GC=40, max_GC=60, min_tm=50, max_tm=70):
        """
        Generates all valid combinations of forward and reverse primers from a specified primer dictionary,
        applying a series of filters based on conservation score, GC content, melting temperature, and amplicon length.
        Parameters:
            dict_name (str): The attribute name of the primer dictionary to use for combinations.
            min_amplicon_length (int, optional): Minimum allowed length of the amplicon (default: 70).
            max_amplicon_length (int, optional): Maximum allowed length of the amplicon (default: 1000).
            conservation_score_threshold (float, optional): Minimum conservation score required for primers (default: 95).
            min_GC (float, optional): Minimum GC content percentage for primers (default: 40).
            max_GC (float, optional): Maximum GC content percentage for primers (default: 60).
            min_tm (float, optional): Minimum melting temperature (°C) for primers (default: 50).
            max_tm (float, optional): Maximum melting temperature (°C) for primers (default: 70).
        Returns:
            None: The resulting valid primer combinations are stored in the `self.primer_combinations` attribute as a list of dictionaries.
        Notes:
            - The method filters primers by conservation score, GC content, and melting temperature before generating combinations.
            - Only combinations where the forward primer is upstream of the reverse primer and the amplicon length is within the specified range are considered.
            - Additional properties such as hairpin, homodimer, heterodimer, GC content, melting temperature, amplicon length, and combination score are calculated for each primer pair.
            - Intermediate filtering steps and distributions are visualized using histograms.
        """

        primers_df = pd.DataFrame(getattr(self, dict_name, None))

        #Filters need to be applied as otherwise the combinations will be too large (over 91 billion) > memory error

        # Filter primers by conservation score > 95 first
        primers_df = primers_df[primers_df['Score'].astype(float) > conservation_score_threshold].copy()
        print(len(primers_df), " primers after conservation score filtering")

        # Filter primers by GC content between 40% and 60%
    
        # Make histogram of GC content percentages and identify peak then determine margins
        plt.figure(figsize=(8, 5))
        plt.hist(primers_df['GC Content'], bins=20, color='skyblue', edgecolor='black')
        plt.title('GC Content Distribution of All Primers')
        plt.xlabel('GC Content (%)')
        plt.ylabel('Number of Primers')
        plt.grid(True)
        plt.show()
        primers_df = primers_df[(primers_df['GC Content'].astype(float) > min_GC) & (primers_df['GC Content'].astype(float) < max_GC)].copy()
        print(len(primers_df), " primers after GC content filtering")

        # Filter primers by melting temperature between 50 and 70 degrees Celsius
        # Make histogram of melting temperatures percentages and identify peak then determine margins
        plt.figure(figsize=(8, 5))
        plt.hist(primers_df['Melting Temperature'], bins=20, color='skyblue', edgecolor='black')
        plt.title('Melting Temperature Distribution of All Primers')
        plt.xlabel('Melting Temperature (°C)')
        plt.ylabel('Number of Primers')
        plt.grid(True)
        plt.show()
        primers_df = primers_df[(primers_df['Melting Temperature'].astype(float) > min_tm) & (primers_df['Melting Temperature'].astype(float) < max_tm)].copy()
        print(len(primers_df), " primers after melting temperature filtering")

        # Separate primers by type
        forward_df = primers_df[primers_df['Type'] == 'Forward'].copy()
        reverse_df = primers_df[primers_df['Type'] == 'Reverse'].copy()
        print(f"Forward primers: {len(forward_df)}, Reverse primers: {len(reverse_df)}")

        # Create all combinations of forward and reverse primers (Cartesian product)
        combos_df = forward_df.assign(key=1).merge(
            reverse_df.assign(key=1),
            on='key',
            suffixes=('_fwd', '_rev'))
        print(len(combos_df), "combinations before filtering")

        # Filter combinations by amplicon length and ensure the Forward primer is upstream of the reverse primer
        combos_df = combos_df[(combos_df['Ref Start_fwd'] < combos_df['Ref Start_rev']) &
                              (combos_df['Ref End_fwd'] < combos_df['Ref Start_rev']) &
                              (combos_df['Ref End_rev'] - combos_df['Ref Start_fwd'] >= min_amplicon_length) &
                              (combos_df['Ref End_rev'] - combos_df['Ref Start_fwd'] <= max_amplicon_length) &
                              (abs(combos_df['Melting Temperature_fwd'] - combos_df['Melting Temperature_rev']) <= 5)]
        print(len(combos_df), "combinations after filtering by amplicon length and primer positions")
        # combos_df.to_excel('combos_df.xlsx', index=False)
        # display(combos_df)

        combos_df['fwd_hairpin'] = combos_df['Sequence_fwd'].apply(self.calc_hairpin)
        combos_df['rev_hairpin'] = combos_df['Sequence_rev'].apply(self.calc_hairpin)
        combos_df['fwd_homodimer'] = combos_df['Sequence_fwd'].apply(self.calc_homodimer)
        combos_df['rev_homodimer'] = combos_df['Sequence_rev'].apply(self.calc_homodimer)
        combos_df['heterodimer'] = combos_df.apply(self.calc_heterodimer, axis=1)

        # Calculate other properties of the combinations
        combos_df['GC_content'] = combos_df[['GC Content_fwd', 'GC Content_rev']].astype(float).mean(axis=1)
        combos_df['Melting Temperature'] = combos_df[['Melting Temperature_fwd', 'Melting Temperature_rev']].astype(float).mean(axis=1)
        combos_df['Amplicon_length'] = combos_df['Ref End_rev'] - combos_df['Ref Start_fwd']
        combos_df['Combination_score'] = combos_df[['Score_fwd', 'Score_rev']].astype(float).mean(axis=1)

        self.primer_combinations = combos_df.to_dict(orient='records')
    
    # Calculate fold scores
    def fold_score(self, delta_g, dg_threshold=-7000, slope=1.0, T=310):
        '''
        Convert ΔG (kcal/mol) into a 0-100 'goodness' score based on thermodynamics principles.

        Parameters:
        delta_g       : ΔG value (negative for binding)
        dg_threshold  : ΔG where penalty becomes significant (kcal/mol)
        slope         : controls steepness of penalty transition
        T             : temperature in Kelvin (default 310K = ~37ºC)

        Notes:
        Higher score = better primer (less unwanted folding)
        This function is called by fold_scores_for_combinations(). Use this function when you just want to calculate the fold score for a single ΔG value.

        '''
        R = 1.987  # cal/(mol·K)

        # Logistic penalty: probability that primer is "lost" to unwanted folding
        penalty = 1 / (1 + np.exp((delta_g - (dg_threshold/2)) / (slope * R * T))) # threshold is halved so its value approximates 0%

        # Convert to a "goodness" percentage: 100 = ideal, 0 = bad
        score = (1 - penalty) * 100
        return np.round(score, 2)
    
    # Calculate fold scores for primer combinations
    def fold_scores_for_combinations(self, dict_name):
        """
        Calculate fold scores for primer combinations based on ΔG values.

        Parameters:
            dict_name (str): The attribute name of the primer combinations dictionary.
        Returns:
            None: The method updates the primer combinations with fold scores and saves them back to the class
        """
        combos_df = pd.DataFrame(getattr(self, dict_name, None))

        # Format outputs from primer3 into floats
        for col in ['fwd_hairpin', 'rev_hairpin', 'fwd_homodimer', 'rev_homodimer', 'heterodimer']:
            combos_df[col] = combos_df[col].astype(str).str.replace(' dg=', '').astype(float)

        
        # Calculate self fold scores for each combination
        combos_df['fwd_self_score'] = self.fold_score(
            combos_df['fwd_hairpin'].values,
            -2000,
            1.0,
            combos_df['Melting Temperature_fwd'].values + 273.15
        )
        combos_df['rev_self_score'] = self.fold_score(
            combos_df['rev_hairpin'].values,
            -2000,
            1.0,
            combos_df['Melting Temperature_rev'].values + 273.15
        )
        combos_df['Self_Fold_Score'] = combos_df[['fwd_self_score','rev_self_score']].astype(float).mean(axis=1).round(2)

        # Calculate homo fold scores for each combination
        combos_df['fwd_homo_score'] = self.fold_score(
            combos_df['fwd_homodimer'].values,
            -5000,
            1.0,
            combos_df['Melting Temperature_fwd'].values + 273.15
        )
        combos_df['rev_homo_score'] = self.fold_score(
            combos_df['rev_homodimer'].values,
            -5000,
            1.0,
            combos_df['Melting Temperature_rev'].values + 273.15
        )
        combos_df['Homo_Fold_Score'] = combos_df[['fwd_homo_score', 'rev_homo_score']].astype(float).mean(axis=1).round(2)
        
        # Calculate dimer fold scores for each combination
        combos_df['Dimer_Fold_Score'] = self.fold_score(
            combos_df['heterodimer'].values,
            -5000,
            1.0,
            combos_df['Melting Temperature'].values + 273.15
        )

        # Calculate No-fold %
        combos_df['No-fold score'] = combos_df[['Self_Fold_Score', 'Homo_Fold_Score', 'Dimer_Fold_Score']].mean(axis=1).round(2)

        # Clean up dictionary by removing redundant/unnecessary columns
        combos_df = combos_df.drop(columns=['Alignment Sequence_fwd',
                                            'Length in Alignment_fwd',
                                            'key',
                                            'Alignment Sequence_rev',
                                            'Length in Alignment_rev',
                                            ])
        # Save back to class
        self.primer_combinations = combos_df.to_dict(orient='records')

#%% Utilities

# Counts the number of sequences in a fasta file
def count_sequences(fasta):
    """
    Counts the number of sequences in a FASTA file.
    Parameters:
        fasta (str or file-like object): Path to the FASTA file or a file-like object containing FASTA-formatted sequences.
    Returns:
        int: The number of sequences found in the FASTA file.
    Prints:
        The number of sequences to the standard output.
    """

    count = sum(1 for _ in SeqIO.parse(fasta, "fasta"))
    print("Number of sequences:", count)
    return count
    
# Change the order of the sequences in a FASTA file
def reorder_sequences(input_fasta, output_fasta, reference_id, save_to):
    """
    Reorders sequences in a FASTA file so that a specified reference sequence appears first, then writes the reordered sequences to a new FASTA file.
    Args:
        input_fasta (str): Path to the input FASTA file containing the sequences.
        output_fasta (str): Name of the output FASTA file to write the reordered sequences.
        reference_id (str): The sequence ID of the reference genome to be placed first.
        save_to (str): Directory path where the output FASTA file will be saved.
    Returns:
        Print: a message indicating success or if the reference genome was not found.
    Side Effects:
        - Writes a new FASTA file with the reordered sequences to the specified directory.
        - Use this to avoid having the reference genome lost inbetween calculations (e.g.: during dupe removal).
    """

    # Read all sequences
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    # Separate the reference genome
    ref_seq = None
    other_seqs = []

    for seq in sequences:
        if seq.id == reference_id:
            ref_seq = seq
        else:
            other_seqs.append(seq)

    # Write the sequences back with the reference genome first
    if ref_seq:
        os.chdir(save_to)
        with open(output_fasta, "w") as output_handle:
            SeqIO.write([ref_seq] + other_seqs, output_handle, "fasta")
        print(f"Reordered FASTA saved as {output_fasta}")
    else:
        print("Reference genome not found in the FASTA file.")
        

# Remove duplicates in a FASTA file (try to have refseq as first entry)

def remove_dupes(input_fasta, output_fasta, save_to):
    """
    Removes duplicate sequences from a FASTA file and writes the unique sequences to a new FASTA file.
    Args:
        input_fasta (str): Path to the input FASTA file containing sequences.
        output_fasta (str): Name of the output FASTA file to write unique sequences.
        save_to (str): Directory path where the output FASTA file will be saved.
    Notes:
        - Only the first occurrence of each unique sequence is retained.
        - Use reorder_sequences() beforehand to ensure a specific sequence (e.g., reference genome) is kept.
    """

    unique_sequences = {}
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequence = str(record.seq)
        if sequence not in unique_sequences:
            unique_sequences[sequence] = record.id
    os.chdir(save_to)      
    with open(output_fasta, 'w') as file:
        for sequence, seq_id in unique_sequences.items():
            file.write(f'>{seq_id}\n{sequence}\n')
            

# For proteins: remove sequences with ambiguous (X) nucl, stop codons (*) or excessive gaps (-)

def remove_ambiguous(input_fasta, output_fasta, save_to):
    os.chdir(save_to)
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if "X" not in record.seq and "*" not in record.seq and record.seq.count("-") < len(record.seq) * 0.05:
                SeqIO.write(record, out_f, "fasta")

    print(f"Ambiguous sequences removed. Cleaned sequences saved to {output_fasta}")
    

def extract_sequence(fasta_file, output, sequence_ID, save_to):
    """
    Extracts a specific sequence from a FASTA file by its sequence ID and writes it to an output file.
    Args:
        fasta_file (str): Path to the input FASTA file containing sequences.
        output (str): Name of the output file where the extracted sequence will be saved.
        sequence_ID (str): The ID of the sequence to extract from the FASTA file.
        save_to (str): Directory path where the output file will be saved.
    """
    
    os.chdir(save_to)
    extract = open(output, 'w')
    for record in SeqIO.parse(fasta_file, 'fasta'):
        if record.id == sequence_ID:
            extract.write('>'+record.id+'\n')
            extract.write(str(record.seq.upper())+'\n')
    extract.close()
    print(f'Sequence extracted to {output}.')
    

# Export Primer IDs and sequences from a dictionary

def export_primers_to_fasta(primer_list, output_file):
    """
    Exports a list of primer dictionaries to a FASTA file using Primer ID as header and Sequence as sequence.

    Args:
        primer_list (list): List of primer dictionaries.
        output_file (str): Path to output FASTA file.
    """
    records = []
    for primer in primer_list:
        seq_id = primer.get("Primer ID", "unknown_id")
        seq_str = primer.get("Sequence", "")
        if seq_str:  # Avoid empty sequences
            record = SeqRecord(Seq(seq_str), id=seq_id, description="")
            records.append(record)
    
    SeqIO.write(records, output_file, "fasta")

def save_to_csv(data, path):
    """
    Saves the given data to a CSV file at the specified path.
    Parameters:
        data (iterable or dict): The name of the list of dictionaries to be saved.
        path (str): The file path where the CSV will be saved.
    Returns:
        None
    Example:
        save_to_csv('viruscope.primer_combinations', 'output.csv')
    """

    df = pd.DataFrame(data)
    df.to_csv(path, index=False)
