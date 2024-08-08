# Okay, the thing is we have to update the VirStrain database every once in a while
# I get the input fasta from GISAID and ENA to make sure that all lineages are represented
    # GISAID : 
    # ENA : 

# For now I'm doing this on googleColab, okok
# I just need to save the script in the repo so I won't forget

# Install biopython
pip install biopython

# Import modules
import os
import re
import pandas as pd
import Bio
from Bio import SeqIO

# Directory containing the metadata and fasta files
directory_path = "/content/drive/MyDrive/Coinfection_Paper/VirStrain_DatabaseFiles"

# I got these metadata by running the input files in Pangolin, to get the accurate lineage assignment for each sequence
# Input metadata containing lineage assignment of fasta from GISAID
df_GISAID = pd.read_csv('/content/drive/MyDrive/Coinfection_Paper/VirStrain_DatabaseFiles/lineage_report_VOI_RepresentativeGenomes_20220815.csv')
# Input metadata containing lineage assignment of fasta from ENA
df_ENA_original = pd.read_csv('/content/drive/MyDrive/Coinfection_Paper/VirStrain_DatabaseFiles/lineage_report_ena_data_20240808.csv')


# I process the input files from GISAID and ENA differently cause their format are different from each other


### GISAID input data wrangling ###
# multifasta from GISAID, input 
input_fasta_GISAID = "/content/drive/MyDrive/Coinfection_Paper/VirStrain_DatabaseFiles/VOI_RepresentativeGenomes_CodonAligned_2022-08-15.fasta"
# and output FASTA (renamed headers)
output_fasta_GISAID = "/content/drive/MyDrive/Coinfection_Paper/VirStrain_DatabaseFiles/fromGISAID.fasta"
# Creating an dictionary to map the 'taxon' values to the 'lineage' values
taxon_to_lineage_GISAID = dict(zip(df_GISAID['taxon'], df_GISAID['lineage']))
# Renaming the FASTA headers by lineage assignment
with open(output_fasta_GISAID, "w") as output_handle:
  for record in SeqIO.parse(input_fasta_GISAID, "fasta"):
    # Get the new header using the directory
    new_header = taxon_to_lineage_GISAID.get(record.id, record.id)
    # Rename the record
    record.id = new_header
    record.description = new_header
    # Write the modified record to the output file
    SeqIO.write(record, output_handle, "fasta")


### ENA input data wrangling ###
# multifasta from ENA, original input file
input_fasta_ENA_original = "/content/drive/MyDrive/Coinfection_Paper/VirStrain_DatabaseFiles/ena_data_20240808-0546.fasta"
# processed input file
input_fasta_ENA = "/content/drive/MyDrive/Coinfection_Paper/VirStrain_DatabaseFiles/fromENA_editedHeaders.fasta"
# and output FASTA (renamed headers)
output_fasta_ENA = "/content/drive/MyDrive/Coinfection_Paper/VirStrain_DatabaseFiles/fromENA.fasta"
# Shortening ENA FASTA headers
with open(input_fasta_ENA, "w") as output_handle:
  for record in SeqIO.parse(input_fasta_ENA_original, "fasta"):
    # Modifying the headers cause they are too looong, retaning only the first part before the space
    new_header = record.id.split()[0]
    # Update the record's ID and description with the new header
    record.id = new_header
    record.description = new_header
    # Write the modified record to the output file
    SeqIO.write(record, output_handle, "fasta")
# Shortening ENA FASTA headers on the metadata csv
# Modify the 'taxon' column to retain only the first part before the space
df_ENA_original['taxon'] = df_ENA_original['taxon'].apply(lambda x: x.split()[0])
# Save the modified metadata to a new CSV file
output_metadata = '/content/drive/MyDrive/Coinfection_Paper/VirStrain_DatabaseFiles/fromENA_editedHeaders.csv'
df_ENA_original.to_csv(output_metadata, index=False)
df_ENA = pd.read_csv(output_metadata)
# Creating an dictionary to map the 'taxon' values to the 'lineage' values
taxon_to_lineage_ENA = dict(zip(df_ENA['taxon'], df_ENA['lineage']))
# Renaming the FASTA headers by lineage assignment
with open(output_fasta_ENA, "w") as output_handle:
  for record in SeqIO.parse(input_fasta_ENA, "fasta"):
    # Get the new header using the directory
    new_header = taxon_to_lineage_ENA.get(record.id, record.id)
    # Rename the record
    record.id = new_header
    record.description = new_header
    # Write the modified record to the output file
    SeqIO.write(record, output_handle, "fasta")


# Then I save the output fasta files and manually concatenate
# After concatenating, I make a new VirStrain Database