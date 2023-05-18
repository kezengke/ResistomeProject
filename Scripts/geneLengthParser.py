import csv

fasta_file = 'nucleotide_fasta_protein_homolog_model.fasta'  
csv_file = 'geneLength.csv'     

sequences = []  # List to store the parsed header and length information

# Open the FASTA file and read its contents
with open(fasta_file, 'r') as file:
    header = ''
    sequence = ''
    
    # Iterate through each line in the file
    for line in file:
        line = line.strip()  # Remove leading/trailing whitespace
        
        if line.startswith('>'):  # Header line
            # If a sequence has been collected, add the header and its length to the list
            if sequence:
                sequences.append([header, len(sequence)])
            
            # Reset header and sequence for the new entry
            header = line[1:]  # Exclude the ">" symbol
            sequence = ''
        else:  # Sequence line
            sequence += line
    
    # Add the last sequence after reaching the end of the file
    if sequence:
        sequences.append([header, len(sequence)])

# Write the sequences list to a CSV file
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Header', 'Length'])  # Write the header row
    writer.writerows(sequences)  # Write the sequence data

