"""Look at conservation in the EL4 flexible loop of Fks1."""

import pandas as pd

# Import fasta
with open('../ortholog_data/EL4.fasta', 'r') as file:
    file_content = file.read()
    
    lines = file_content.split('\n')

headers = [h[1:] for h in lines if h[0]=='>']
sequences = [h for h in lines if h[0]!='>']        

df = pd.DataFrame(list(zip(headers,sequences)), columns=['header','sequence'])

# Check that EL4 is delimited by the two conserved cysteines
df['C_both_sides'] = df.sequence.apply(lambda x: True if (x[0]=='C') & (x[-1]=='C') else False)

# Filter
df['length'] = df.sequence.apply(lambda x: len(x) - x.count('-'))
df_filtered = df[(df.C_both_sides) 
                 #& (df.length == 18)
                 ].reset_index(drop=True)

# Count arginines and lysines
df_filtered['R_count'] = df_filtered.sequence.str.count('R')
df_filtered['K_count'] = df_filtered.sequence.str.count('K')
df_filtered['RK_count'] = df_filtered.R_count + df_filtered.K_count
gby = df_filtered.groupby('sequence')[['header','RK_count']].agg(n_homologues = ('header','nunique'),
                                                           RK_count = ('RK_count', lambda x: x.astype(bool).sum())
                                                          )
