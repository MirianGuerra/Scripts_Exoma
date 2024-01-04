import pyoncoprint
import pandas as pd
import gget

df = pd.read_csv('mafzao_sem_variant_dbscan.csv', sep=';')

df['gene'] = df['ENSP'].apply(lambda x: gget.info(x))
print(df)
df.to_csv('df_com_gene.csv', sep=';')

df['track_name'] = df['Reference_Allele'] + '>' + df['Tumor_Seq_Allele2']

agrupar_track_name = df.groupby('track_name')['AD']
oncoplot = df[['AD', 'gene', 'Variant_Classification']]
new_df = oncoplot.pivot(index='gene', columns='AD', values='Variant_Classification')
new_df = new_df.fillna('')