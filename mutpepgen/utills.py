from typing import List, Tuple, Dict, Set
from datetime import datetime
import logging
import pandas as pd
import re
import os

logger = logging.getLogger(__name__)
# Constants
INIT_VARIABLES: Dict[str, int] = {
    'MAIN_counter': 0,
    'FILE_ERROR_counter': 0,
    'FILE_SUCCESS_counter': 0,
    'TRANSCRIPT_FOUND_counter': 0,
    'TRANSCRIPT_NOT_FOUND_counter': 0,
    'SUBSTITUTION_FOUND_counter': 0,
    'SUBSTITUTION_NOT_FOUND_counter': 0,
    'SUBSTITUTION_SUCCESS_counter': 0,
    'SUBSTITUTION_ERROR_counter': 0,
    'POSITION_FOUND_counter': 0,    
    'POSITION_2ND_ATTEMPT_FOUND_counter': 0,
    'POSITION_3RD_ATTEMPT_FOUND_counter': 0,
    'POSITION_NOT_FOUND_counter': 0,
    'UNIPROTtoGRch38_NOT_FOUND_counter': 0,
    'MULTI_SEQ_FOUND_counter': 0,
    'MULTI_SEQ_POSITION_FOUND_counter': 0,
    'SUBSTITUTION_FOUND_2ND_ATTEMPT_counter': 0,
    'SUBSTITUTION_FOUND_3RD_ATTEMPT_counter': 0,
}
class dTime:
    """
    Simple class to return a timestamp string for logging or filenames.
    Format: YYYYMMDD_HHMMSS
    """
    def __init__(self):
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    def __str__(self):
        return self.timestamp

    def now(self):
        """Return current timestamp string"""
        return datetime.now().strftime("%Y%m%d_%H%M%S")
    
class Files_Manager:
    def __init__(self):
        self.file = None
    def check_file(self):
        return self.file
    #https://stackoverflow.com/questions/40745686/python-process-file-using-multiple-cores
    @staticmethod
    def file_block(fp, number_of_blocks, block):
        '''
        A generator that splits a file into blocks and iterates
        over the lines of one of the blocks.

        '''

        assert 0 <= block and block < number_of_blocks
        assert 0 < number_of_blocks

        fp.seek(0,2)
        file_size = fp.tell()

        ini = file_size * block / number_of_blocks
        end = file_size * (1 + block) / number_of_blocks

        if ini <= 0:
            fp.seek(0)
        else:
            fp.seek(ini-1)
            fp.readline()

        while fp.tell() < end:
            yield fp.readline()
            
            
    def _find_fileType(self):
        if len(self.file.split('.')) > 1:
            return self.file.split('.')[-1]
        else:
            raise ValueError("File type not found")
    
    def check_permission(self):
        return self.file.mode
    
    def check_file_name(self):
        return self.file.name
    
    def check_file_encoding(self):
        return self.file.encoding
    
    def check_file_type(self):
        return self.file.mode
            
    def open_file(self, file_path):
        self.file = open(file_path, 'r')
        return self.file

    def read_file(self):
        return self.file.read()

    def close_file(self):
        self.file.close()
        self.file = None

class logs:
    def __init__(self):
        self.cwd = os.getcwd()
        self.logs = None
        self.dTime = dTime()
        
    def _init_logsFile(self):
        file_path = os.path.join(self.cwd, f'{self.dTime}_proccess.logs')
        logging.basicConfig(filename=file_path, filemode='w', format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)
        
    def _add_info(self, info):
        logging.info(info)
    
    def _add_error(self, error):
        logging.error(error)
        
    def _add_warning(self, warning):
        logging.warning(warning)
        
    def _add_debug(self, debug):
        logging.debug(debug)

class UniProtParser:
    """
    Parser for UniProt database format to extract Ensembl transcript IDs
    and map them to protein sequences.
    """
    
    def __init__(self, log_callback=None):
        """
        Initialize the parser
        
        Args:
            log_callback: Function to call for logging messages
        """
        self.sequences = {}
        self.log = log_callback if log_callback else print
    
    def parse_file(self, file_path):
        """
        Parse a UniProt format file and extract Ensembl IDs and sequences
        
        Args:
            file_path: Path to the UniProt format file
            
        Returns:
            Dictionary of ENST IDs to sequences
        """
        self.log(f"Parsing UniProt format file: {file_path}")
        
        try:
            # Determine delimiter based on file extension
            file_ext = os.path.splitext(file_path)[1].lower()
            if file_ext in ['.tsv', '.txt']:
                delimiter = '\t'
            else:  # .csv
                delimiter = ','
            
            # Read the file
            df = pd.read_csv(file_path, delimiter=delimiter, low_memory=False)
            
            # First, identify the relevant columns
            ensembl_col = self._find_ensembl_column(df)
            sequence_col = self._find_sequence_column(df)
            
            if not ensembl_col:
                # Look for a column containing Ensembl IDs in the data
                for col in df.columns:
                    sample = df[col].dropna().astype(str).head(10)
                    if any('ENST' in str(val) for val in sample):
                        ensembl_col = col
                        break
            
            if not ensembl_col:
                self.log("No column containing Ensembl IDs found")
                return self.sequences
            
            if not sequence_col:
                # Try to find the sequence column by checking for long string content
                for col in df.columns:
                    if col != ensembl_col:  # Skip the Ensembl column
                        sample = df[col].dropna().astype(str).head(5)
                        if any(len(str(val)) > 50 and self._is_likely_protein_sequence(str(val)) for val in sample):
                            sequence_col = col
                            break
            
            if not sequence_col:
                self.log("No sequence column found")
                return self.sequences
            
            self.log(f"Using column '{ensembl_col}' for Ensembl IDs and '{sequence_col}' for sequences")
            
            # Process each row
            count = 0
            for idx, row in df.iterrows():
                if pd.isna(row[ensembl_col]) or pd.isna(row[sequence_col]):
                    continue
                
                # Get the sequence
                sequence = str(row[sequence_col]).strip()
                
                # Parse Ensembl IDs from the Ensembl column
                ensembl_text = str(row[ensembl_col])
                
                # For UniProt format, Ensembl IDs are often in quotes and semicolon-separated
                # Example: "ENST00000436697.3; ENSP00000484893.1; ENSG00000225973.4.";"ENST00000567948.1; ...
                
                # Extract all ENST IDs
                enst_ids = self._extract_all_enst_ids(ensembl_text)
                
                if enst_ids:
                    for enst_id in enst_ids:
                        # Remove version number if present
                        if "." in enst_id:
                            enst_id = enst_id.split(".")[0]
                        self.sequences[enst_id] = sequence
                        count += 1
            
            self.log(f"Successfully mapped {count} ENST IDs to sequences from UniProt format")
            
        except Exception as e:
            self.log(f"Error parsing UniProt format file: {str(e)}")
            raise
        
        return self.sequences
    
    def _find_ensembl_column(self, df):
        """Find the column containing Ensembl IDs"""
        ensembl_keywords = ['ensembl', 'enst', 'transcript']
        
        for col in df.columns:
            col_lower = str(col).lower()
            if any(keyword in col_lower for keyword in ensembl_keywords):
                return col
        
        return None
    
    def _find_sequence_column(self, df):
        """Find the column containing protein sequences"""
        sequence_keywords = ['sequence', 'seq', 'protein']
        
        for col in df.columns:
            col_lower = str(col).lower()
            if any(keyword in col_lower for keyword in sequence_keywords):
                return col
        
        return None
    
    def _extract_all_enst_ids(self, text):
        """Extract all ENST IDs from text"""
        # Handle UniProt format with quotes and semicolons
        # First, find all quoted sections
        quoted_sections = re.findall(r'"([^"]*)"', text)
        
        # For each section, find ENST IDs
        enst_ids = []
        
        # Process quoted sections
        for section in quoted_sections:
            matches = re.findall(r'(ENST\d+(?:\.\d+)?)', section)
            enst_ids.extend(matches)
        
        # Also look for ENSTs not in quotes
        matches = re.findall(r'(ENST\d+(?:\.\d+)?)', text)
        enst_ids.extend(matches)
        
        # Remove duplicates
        return list(set(enst_ids))
    
    def _is_likely_protein_sequence(self, text):
        """Check if a string is likely to be a protein sequence"""
        # Protein sequences typically contain only amino acid letters
        amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
        
        # Clean the string and convert to uppercase
        text = re.sub(r'\s', '', text).upper()
        
        # Check if at least 80% of characters are amino acids
        if len(text) == 0:
            return False
            
        amino_acid_count = sum(1 for char in text if char in amino_acids)
        return amino_acid_count / len(text) >= 0.8
    
    def save_to_fasta(self, output_path):
        """Save the parsed sequences to a FASTA file"""
        with open(output_path, 'w') as f:
            for enst_id, sequence in self.sequences.items():
                f.write(f">{enst_id}\n{sequence}\n")
                
        self.log(f"Saved {len(self.sequences)} sequences to {output_path}")


# if __name__ == "__main__":
#     parser = UniProtParser()
#     sequences = parser.parse_file("uniprot_data.tsv")
#     parser.save_to_fasta("uniprot_sequences.fasta")



def parse_maf_file(file_path):
    """
    Parse a MAF (Mutation Annotation Format) file and extract ENST IDs and protein mutations.
    
    Args:
        file_path (str): Path to the MAF file.
        
    Returns:
        list: A list of dictionaries containing mutation information.
    """
    mutations = []
    header = None
    
    with open(file_path, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
                
            # Parse the header line
            if header is None:
                header = line.strip().split('\t')
                continue
                
            # Parse data lines
            fields = line.strip().split('\t')
            if len(fields) < len(header):
                # Pad with empty strings if the line has fewer fields than the header
                fields.extend([''] * (len(header) - len(fields)))
                
            # Create a dictionary mapping column names to values
            mutation_data = dict(zip(header, fields))
            mutations.append(mutation_data)
    
    return mutations

def extract_transcripts_and_mutations(mutations):
    """
    Extract ENST IDs and protein mutations from parsed MAF data.
    
    Args:
        mutations (list): List of mutation dictionaries from parse_maf_file.
        
    Returns:
        dict: Dictionary mapping ENST IDs to their protein mutations.
    """
    transcript_to_mutations = {}
    
    for mutation in mutations:
        # Get the primary transcript ID
        primary_transcript = mutation.get('Transcript_ID', '')
        
        # Extract information from the all_effects field
        all_effects = mutation.get('all_effects', '')
        
        # Process the primary transcript if it's an ENST ID
        if primary_transcript.startswith('ENST'):
            if primary_transcript not in transcript_to_mutations:
                transcript_to_mutations[primary_transcript] = []
            
            protein_change = mutation.get('HGVSp_Short', '') or mutation.get('HGVSp', '')
            if not protein_change:
                # Check if protein position information is available
                protein_position = mutation.get('Protein_position', '')
                amino_acids = mutation.get('Amino_acids', '')
                if protein_position and amino_acids:
                    protein_change = f"p.{amino_acids}{protein_position}"
            
            if protein_change:
                transcript_to_mutations[primary_transcript].append({
                    'gene': mutation.get('Hugo_Symbol', ''),
                    'protein_change': protein_change,
                    'variant_classification': mutation.get('Variant_Classification', ''),
                    'hgvsc': mutation.get('HGVSc', ''),
                })
        
        # Process additional transcripts from all_effects field
        if all_effects:
            effects = all_effects.split(';')
            for effect in effects:
                if not effect:
                    continue
                
                effect_parts = effect.split(',')
                if len(effect_parts) >= 4:  # Ensure there are enough parts
                    gene = effect_parts[0]
                    effect_type = effect_parts[1]
                    protein_change = effect_parts[2]
                    transcript_id = ''
                    
                    # Find the ENST ID in the effect parts
                    for part in effect_parts:
                        if part.startswith('ENST'):
                            transcript_id = part
                            break
                    
                    if transcript_id.startswith('ENST'):
                        if transcript_id not in transcript_to_mutations:
                            transcript_to_mutations[transcript_id] = []
                        
                        transcript_to_mutations[transcript_id].append({
                            'gene': gene,
                            'protein_change': protein_change,
                            'effect_type': effect_type,
                        })
    
    return transcript_to_mutations

def get_transcripts_with_protein_mutations(file_path):
    """
    Extract ENST IDs and their corresponding protein mutations from a MAF file.
    
    Args:
        file_path (str): Path to the MAF file.
        
    Returns:
        dict: Dictionary mapping ENST IDs to their protein mutations.
    """
    mutations = parse_maf_file(file_path)
    return extract_transcripts_and_mutations(mutations)

def display_transcripts_and_mutations(transcript_to_mutations):
    """
    Display the ENST IDs and their corresponding protein mutations.
    
    Args:
        transcript_to_mutations (dict): Dictionary mapping ENST IDs to their protein mutations.
    """
    print(f"Found {len(transcript_to_mutations)} transcripts with mutations:")
    
    for transcript_id, mutations in transcript_to_mutations.items():
        print(f"\n{transcript_id}:")
        for mutation in mutations:
            gene = mutation.get('gene', 'N/A')
            protein_change = mutation.get('protein_change', 'N/A')
            effect_type = mutation.get('effect_type', mutation.get('variant_classification', 'N/A'))
            
            print(f"  - Gene: {gene}, Protein change: {protein_change}, Effect: {effect_type}")



# Example to parse a specific MAF file and save the results to a CSV
def save_results_to_csv(transcript_to_mutations, output_file):
    """
    Save the transcript and mutation data to a CSV file.
    
    Args:
        transcript_to_mutations (dict): Dictionary mapping ENST IDs to their protein mutations.
        output_file (str): Path to the output CSV file.
    """
    import csv
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Transcript_ID', 'Gene', 'Protein_Change', 'Effect_Type'])
        
        for transcript_id, mutations in transcript_to_mutations.items():
            for mutation in mutations:
                writer.writerow([
                    transcript_id,
                    mutation.get('gene', 'N/A'),
                    mutation.get('protein_change', 'N/A'),
                    mutation.get('effect_type', mutation.get('variant_classification', 'N/A'))
                ])
    
    print(f"Results saved to {output_file}")

# Additional usage example:
# Parse a MAF file and save the results to a CSV
# Save results to CSV
# maf_file = "data/e0ad030c-034a-4315-a8a8-fa5be7feaac4.wxs.aliquot_ensemble_masked.maf"
# output_file = "data/out_transcript_mutations.csv"
# transcript_to_mutations = get_transcripts_with_protein_mutations(maf_file)
# # display_transcripts_and_mutations(transcript_to_mutations)
# save_results_to_csv(transcript_to_mutations, output_file)
