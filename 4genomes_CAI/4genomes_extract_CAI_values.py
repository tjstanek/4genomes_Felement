import os

# Define the directory where the normalized CAI output files are located
input_directory = 'normalized_CAI_output'
output_file = 'normalized_CAI_summary.txt'

# Open the output file in write mode
with open(output_file, mode='w') as output_txt:
    # Write the header to the output file
    output_txt.write('transcript\ttrue_CAI\tnorm_CAI\tspecies\telement\n')

    # Iterate through all files in the normalized_CAI_output directory
    for filename in os.listdir(input_directory):
        if filename.endswith(".txt"):  # Only process .txt files
            # Extract species and element information from the filename
            species = filename.split('_')[0]  # First string before the first underscore
            element = f"Muller_{filename.split('_')[1]}"  # String between the second and third underscore
            
            # Construct the full path to the current file
            file_path = os.path.join(input_directory, filename)
            
            # Open the current file and process it
            with open(file_path, mode='r') as file:
                lines = file.readlines()
                
                # Process each line in the file
                for line in lines:
                    # Skip empty lines and header lines
                    if line.startswith("seq_name") or not line.strip():
                        continue
                    
                    # Split the line by tabs
                    columns = line.strip().split('\t')
                    
                    # The seq_name is the first column and normalised_cai is the 7th column
                    seq_name = columns[0]
                    true_cai = columns[3] # The true_cai is at index 4
                    normalised_cai = columns[6]  # The normalised_cai is at index 6
                    
                    # Write the seq_name, normalised_cai, species, and element to the output file
                    output_txt.write(f"{seq_name}\t{true_cai}\t{normalised_cai}\t{species}\t{element}\n")

print("Processing complete. The summary has been written to", output_file)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_principal_isoforms(file_path):
    """
    Parse the APPRIS principal isoform table and return a dictionary
    with gene_id as keys and a list of (isoform_id, isoform_type) tuples as values.
    """
    gene_data = {}
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            gene_id = parts[5].rsplit('-', 1)[0]
            isoform_id = parts[5]
            isoform_type = parts[4]
            if gene_id not in gene_data:
                gene_data[gene_id] = []
            gene_data[gene_id].append((isoform_id, isoform_type))
    return gene_data

def find_longest_principal_isoform_from_gff(gene_data, gff_file_path):
    """
    Find the longest principal isoform for each gene using the GFF file.
    For genes in the APPRIS data, select the longest PRINCIPAL isoform.
    For novel genes not in the APPRIS data, select the longest isoform.
    """
    # Track all principal isoforms from APPRIS
    principal_isoforms = {}
    for gene_id, isoforms in gene_data.items():
        for isoform_id, isoform_type in isoforms:
            if 'PRINCIPAL' in isoform_type:
                if gene_id not in principal_isoforms:
                    principal_isoforms[gene_id] = []
                principal_isoforms[gene_id].append(isoform_id)

    # Track all isoform lengths, including those for novel genes
    isoform_lengths = {}
    gene_isoforms = {}

    with open(gff_file_path, "r") as f:
        for line in f:
            if not line.startswith('#'):
                parts = line.strip().split('\t')
                feature_type = parts[2]
                attributes = dict(item.split('=') for item in parts[8].split(';') if item)
                if feature_type == 'CDS' and 'Parent' in attributes:
                    isoform_id = attributes['Parent']
                    gene_id = isoform_id.rsplit('-', 1)[0]
                    
                    # Track all isoforms for the gene
                    if gene_id not in gene_isoforms:
                        gene_isoforms[gene_id] = set()
                    gene_isoforms[gene_id].add(isoform_id)
                    
                    # Track lengths
                    if isoform_id not in isoform_lengths:
                        isoform_lengths[isoform_id] = 0
                    isoform_lengths[isoform_id] += int(parts[4]) - int(parts[3]) + 1

    # Identify the longest isoform for each gene
    longest_isoforms = {}
    for gene_id, isoforms in gene_isoforms.items():
        if gene_id in principal_isoforms:
            # For genes in the APPRIS data, select the longest PRINCIPAL isoform
            longest_isoform = None
            max_length = 0
            for isoform_id in principal_isoforms[gene_id]:
                length = isoform_lengths.get(isoform_id, 0)
                if (length > max_length) or (length == max_length and (longest_isoform is None or isoform_id < longest_isoform)):
                    max_length = length
                    longest_isoform = isoform_id
        else:
            # For novel genes, select the longest isoform
            longest_isoform = None
            max_length = 0
            for isoform_id in isoforms:
                length = isoform_lengths.get(isoform_id, 0)
                if (length > max_length) or (length == max_length and (longest_isoform is None or isoform_id < longest_isoform)):
                    max_length = length
                    longest_isoform = isoform_id
        
        if longest_isoform is not None:
            longest_isoforms[gene_id] = longest_isoform

    return longest_isoforms


def extract_related_gff_entries(gff_file, longest_isoforms):
    """
    Extract related GFF entries for the longest principal isoforms.
    """
    related_entries = []
    intergenic_entries = []
    gene_ids_to_include = set(longest_isoforms.keys())
    isoform_ids_to_include = set(longest_isoforms.values())
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            feature_type = parts[2]
            attributes = parts[8]
            if feature_type == 'intergenic':
                intergenic_entries.append(line)
            elif feature_type == 'gene':
                attributes_dict = dict(item.split('=') for item in attributes.split(';') if item)
                if 'ID' in attributes_dict:
                    feature_id = attributes_dict['ID']
                    if feature_id in gene_ids_to_include:
                        related_entries.append(line)
            elif feature_type == 'mRNA':
                attributes_dict = dict(item.split('=') for item in attributes.split(';') if item)
                if 'ID' in attributes_dict:
                    feature_id = attributes_dict['ID']
                    if feature_id in isoform_ids_to_include:
                        related_entries.append(line)
            else:
                attributes_dict = dict(item.split('=') for item in attributes.split(';') if item)
                if 'Parent' in attributes_dict:
                    parent_id = attributes_dict['Parent']
                    if parent_id in isoform_ids_to_include or parent_id in gene_ids_to_include:
                        related_entries.append(line)
                elif 'ID' in attributes_dict:
                    feature_id = attributes_dict['ID']
                    if feature_id in isoform_ids_to_include or feature_id in gene_ids_to_include:
                        related_entries.append(line)
    
    return related_entries, intergenic_entries

def write_gff_entries(output_gff_file, related_entries, intergenic_entries):
    """
    Write the related and intergenic GFF entries to a new file.
    """
    with open(output_gff_file, 'w') as f:
        for entry in related_entries:
            f.write(entry)
        for entry in intergenic_entries:
            f.write(entry)

def extract_and_write_cds_sequences(longest_isoforms, gff_file, genome_fasta_file, output_nucleotide_fasta_D, output_peptide_fasta_D, output_nucleotide_fasta_F, output_peptide_fasta_F):
    """
    Extract the full-length CDS sequences for the longest principal isoforms and write to nucleotide and peptide FASTA files.
    """
    fasta_sequences = {record.id: record for record in SeqIO.parse(genome_fasta_file, "fasta")}
    print(fasta_sequences)
    nucleotide_records_D = []
    peptide_records_D = []
    nucleotide_records_F = []
    peptide_records_F = []

    for gene_id, isoform_id in longest_isoforms.items():
        cds_sequences = []
        chromosome = None
        with open(gff_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if parts[2] == 'CDS':
                    attributes = dict(item.split('=') for item in parts[8].split(';') if item)
                    parent_id = attributes['Parent']
                    if parent_id == isoform_id:
                        start = int(parts[3]) - 1
                        end = int(parts[4])
                        strand = parts[6]
                        chromosome = parts[0]
                        seq = fasta_sequences[chromosome].seq[start:end]
                        #print(seq)
                        if strand == '-':
                            seq = seq.reverse_complement()
                        cds_sequences.append(seq)
        
        if cds_sequences:
            full_cds_seq = Seq('').join(cds_sequences)
            nucleotide_record = SeqRecord(full_cds_seq, id=isoform_id, description='')
            peptide_record = SeqRecord(full_cds_seq.translate(to_stop=True), id=isoform_id, description='')
            if chromosome == "Muller_D":
                nucleotide_records_D.append(nucleotide_record)
                peptide_records_D.append(peptide_record)
            elif chromosome == "Muller_F":
                nucleotide_records_F.append(nucleotide_record)
                peptide_records_F.append(peptide_record)
    
    SeqIO.write(nucleotide_records_D, output_nucleotide_fasta_D, "fasta")
    SeqIO.write(peptide_records_D, output_peptide_fasta_D, "fasta")
    SeqIO.write(nucleotide_records_F, output_nucleotide_fasta_F, "fasta")
    SeqIO.write(peptide_records_F, output_peptide_fasta_F, "fasta")

# Main script execution
appris_file = "Dmel_appris_data.principal_isoforms2.txt"
input_gff_file = "DtakHiC1_Annotations-2024_06_29.clean.gff3"
genome_fasta_file = "DtakHiC1_Muller_DF_original.fasta"
output_gff_file = "Dtak-principal-r6.55.gff3"
output_nucleotide_fasta_D = "Dtak_D_principal.nuc.fna"
output_peptide_fasta_D = "Dtak_D_principal.peptides.faa"
output_nucleotide_fasta_F = "Dtak_F_principal.nuc.fna"
output_peptide_fasta_F = "Dtak_F_principal.peptides.faa"

# Step 1: Parse the APPRIS principal isoform table
gene_data = parse_principal_isoforms(appris_file)

# Step 2: Find the longest principal isoform for each gene using the GFF file
longest_isoforms = find_longest_principal_isoform_from_gff(gene_data, input_gff_file)

# Step 3: Extract related GFF entries for the longest principal isoforms
related_entries, intergenic_entries = extract_related_gff_entries(input_gff_file, longest_isoforms)

# Step 4: Write the related and intergenic GFF entries to a new file
write_gff_entries(output_gff_file, related_entries, intergenic_entries)

# Step 5: Extract and write the CDS sequences to nucleotide and peptide FASTA files
extract_and_write_cds_sequences(longest_isoforms, input_gff_file, genome_fasta_file, output_nucleotide_fasta_D, output_peptide_fasta_D, output_nucleotide_fasta_F, output_peptide_fasta_F)

print('Done!')
