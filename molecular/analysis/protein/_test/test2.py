
import pandas as pd

# Read in secondary structure and convert to correct format
df = pd.read_csv('test.csv').set_index('structure').T.rename_axis('residue')

# Convert from condensed to long format
df = df.stack().to_frame('code').reset_index()

# Cross-tabulate
df['count'] = 1
print(df.pivot_table(index='residue', columns='code', values='count', aggfunc='sum').fillna(0))
