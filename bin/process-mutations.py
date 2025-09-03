#!/usr/bin/env python3

# written by Adeliza Realingo (04 March 2025)

# Import modules
import pandas as pd
import os
import re
import sys

# Input mutations tsv from process freyja_get_lineage_def
# Input: Folder path containing two files
sample = sys.argv[1]

# Assign each file to A_lineage_mutations_path and B_lineage_mutations_path
A_lineage_mutations_path = sys.argv[2]
B_lineage_mutations_path = sys.argv[3]


# Get value for first (A) and second (B) lineage, then use the value as column name
A_lineage = (os.path.basename(A_lineage_mutations_path).split('_')[-1]).rsplit('.', 1)[0]
A_lineage_mutations = pd.read_csv( A_lineage_mutations_path, sep='\t', header=None)
A_lineage_mutations.columns = [A_lineage]
# ~ #
B_lineage = (os.path.basename(B_lineage_mutations_path).split('_')[-1]).rsplit('.', 1)[0]
B_lineage_mutations = pd.read_csv( B_lineage_mutations_path, sep='\t', header=None)
B_lineage_mutations.columns = [B_lineage]


# Only retain nonsynonymous mutations
A_filter_nonsynonymous_mutations = A_lineage_mutations.loc[A_lineage_mutations[A_lineage].astype(str).str.match(r'[A-Za-z0-9]+\([^)]*\)')]
A_filter_nonsynonymous_mutations.reset_index(drop=True, inplace=True)
# ~ #
B_filter_nonsynonymous_mutations = B_lineage_mutations.loc[B_lineage_mutations[B_lineage].astype(str).str.match(r'[A-Za-z0-9]+\([^)]*\)')]
B_filter_nonsynonymous_mutations.reset_index(drop=True, inplace=True)


## Create mutations dataframe for both lineages

# Data wrangling to fill the positions and RNA_change columns
A_positions = A_filter_nonsynonymous_mutations[A_lineage].str.extract(r'([0-9]+)')[0]
A_RNA_change = A_filter_nonsynonymous_mutations[A_lineage].str.extract(r'([A-Za-z0-9]+)')[0]
# ~ #
B_positions = B_filter_nonsynonymous_mutations[B_lineage].str.extract(r'([0-9]+)')[0]
B_RNA_change = B_filter_nonsynonymous_mutations[B_lineage].str.extract(r'([A-Za-z0-9]+)')[0]

# Function to transform RNA change to VCF
def transform_variant(variant):
    match = re.match(r"([A-Z])(\d+)([A-Z])", variant)
    if match:
        return (match.group(2), match.group(1), match.group(3))
    return None

# Apply the function to the DataFrame
A_VCF = A_RNA_change.apply(transform_variant)
B_VCF = B_RNA_change.apply(transform_variant)


A_gene_prot_change = A_filter_nonsynonymous_mutations[A_lineage].str.extract(r'(\([A-Za-z0-9]+:[A-Za-z0-9]+\))')
A_gene = A_gene_prot_change[0].str.extract(r'([A-Za-z0-9]+)')[0]
A_protein_change = A_gene_prot_change[0].str.extract(r'(:[A-Za-z0-9]+)')[0]
A_protein_change = A_protein_change.str.replace(':', '')
# ~ #
B_gene_prot_change = B_filter_nonsynonymous_mutations[B_lineage].str.extract(r'(\([A-Za-z0-9]+:[A-Za-z0-9]+\))')
B_gene = B_gene_prot_change[0].str.extract(r'([A-Za-z0-9]+)')[0]
B_protein_change = B_gene_prot_change[0].str.extract(r'(:[A-Za-z0-9]+)')[0]
B_protein_change = B_protein_change.str.replace(':', '')

# Finally create the mutations dataframe following the mutations.tsv format
A_df_mutations = pd.DataFrame({'Position':A_positions, 'RNA_change':A_RNA_change, 'VCF':A_VCF, 'gene':A_gene, 'protein_change':A_protein_change, A_lineage:'yes'})
B_df_mutations = pd.DataFrame({'Position':B_positions, 'RNA_change':B_RNA_change, 'VCF':B_VCF, 'gene':B_gene, 'protein_change':B_protein_change, B_lineage:'yes'})

# Merge the dataframes for both mutations
df_mutations = pd.merge(A_df_mutations, B_df_mutations, on=['Position', 'RNA_change', 'VCF', 'gene', 'protein_change'], how='outer')
common_rows = df_mutations[A_lineage].eq('yes') & df_mutations[B_lineage].eq('yes')
df_mutations.loc[common_rows, [A_lineage, B_lineage]] = 'yes'
df_mutations['Position'] = df_mutations['Position'].astype(int)
df_mutations = df_mutations.sort_values(by=['Position'])

# Export the dataframe to tsv file
mutations_tsv_path = f'{sample}_mutations.tsv'
df_mutations.to_csv( mutations_tsv_path, sep='\t', index=False)