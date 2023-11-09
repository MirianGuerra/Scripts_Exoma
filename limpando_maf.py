import pandas as pd
import os

itens = os.listdir('.')
maf_files = [item for item in itens if item.endswith('.maf')]
print(maf_files)

relacao_id_ad = pd.read_csv('Id_Ad.csv', sep=';')
relacao_id_ad['ID']=relacao_id_ad['Codigo'].apply(lambda x: x.split('-')[0])
relacao_id_ad.drop(columns='Codigo', inplace=True)
print(relacao_id_ad)

mafzao =  pd.DataFrame()
#mafzao_sem_variant = pd.DataFrame()
#mafzao_sem_variant_dbscan = pd.DataFrame()

for file in maf_files:
    df = pd.read_csv(f'{file}', sep='\t')
    df['ID'] = file.split('-')[0]
    #print(df)
    mafzao = pd.concat([mafzao, df], ignore_index = True)
    
mafzao = mafzao.merge(relacao_id_ad, how = 'left', on='ID')
mafzao_sem_variant = mafzao[~mafzao['Variant_Classification'].isin(['intron_variant', 'synonymous_variant', '3_prime_UTR_variant'])]    
mafzao_sem_variant_dbscan = mafzao_sem_variant[mafzao_sem_variant['dbSNP_RS'].isin(['.'])]

mafzao.to_csv('mafzao.csv', sep=';', index=False)
mafzao_sem_variant.to_csv('mafzao_sem_variant.csv', sep=';', index=False)
mafzao_sem_variant_dbscan.to_csv('mafzao_sem_variant_dbscan.csv', sep=';', index=False)

print(mafzao_sem_variant_dbscan)

print(f'Mafzao --> {len(mafzao)}')
print(f'Mafzao Variant--> {len(mafzao_sem_variant)}')
print(f'Mafzao Variant e DbScan --> {len(mafzao_sem_variant_dbscan)}')

