import pandas as pd
import os
from comut import comut
from comut import fileparsers
import json
import tarfile
import re
"""
itens = os.listdir('Scripts/maf/')
maf_files = [item for item in itens if item.endswith('.maf')]
#print(maf_files)

relacao_id_ad = pd.read_csv('Scripts/Id_Ad.csv', sep=';')
relacao_id_ad['ID']=relacao_id_ad['Codigo'].apply(lambda x: x.split('-')[0])
relacao_id_ad.drop(columns='Codigo', inplace=True)
#print(relacao_id_ad)

mafzao =  pd.DataFrame()
#mafzao_sem_variant = pd.DataFrame()
#mafzao_sem_variant_dbscan = pd.DataFrame()

for file in maf_files:
    df = pd.read_csv(f'Scripts/maf/{file}', sep='\t')
    df['ID'] = file.split('-')[0]
    #print(df)
    mafzao = pd.concat([mafzao, df], ignore_index = True)
    
mafzao = mafzao.merge(relacao_id_ad, how = 'left', on='ID')
mafzao_sem_variant = mafzao[~mafzao['Variant_Classification'].isin(['intron_variant', 'synonymous_variant', '3_prime_UTR_variant'])]    
mafzao_sem_variant_dbscan = mafzao_sem_variant[mafzao_sem_variant['dbSNP_RS'].isin(['.'])]
linhas = []
mafzao_sem_variant_dbscan.reset_index(drop=True, inplace=True)
for i in range(0,len(mafzao_sem_variant_dbscan)-1):
    #print(i)
    if mafzao_sem_variant_dbscan['Variant_Classification'][i]=='5_prime_UTR_variant':
        #print(mafzao_sem_variant_dbscan['Variant_Classification'][i])
        try:
            hgvsc_string = mafzao_sem_variant_dbscan['HGVSc'][i].split('.')[2]
            #print(hgvsc_string)
            numeros = re.findall(r'\d+', hgvsc_string)[0]
            numeros_int = int(numeros)
            #print(f'{i} - {numeros_int}')
            if numeros_int>50:
                linhas.append(i)
        except:
            continue
        
mafzao_sem_variant_dbscan_5UTR = mafzao_sem_variant_dbscan.drop(linhas)
mafzao_sem_variant_dbscan_5UTR.to_csv('Scripts/mafzao_sem_variant_dbscan_5UTR.maf', sep='\t', index=False)
mafzao_sem_variant_dbscan_5UTR_Depth50 = mafzao_sem_variant_dbscan_5UTR[mafzao_sem_variant_dbscan_5UTR['Total_Depth']>50]
mafzao_sem_variant_dbscan_5UTR_Depth50.to_csv('Scripts/mafzao_sem_variant_dbscan_5UTR_Depth50.maf', sep='\t', index=False)

#print(linhas)

#mafzao_sem_variant_dbscan_Hugo_com_ponto = mafzao_sem_variant_dbscan[~mafzao_sem_variant_dbscan['Hugo_Symbol'].isin(['.'])]
#mafzao_limpo = mafzao_sem_variant_dbscan_Hugo_com_ponto[["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele","Variant_Classification", "Variant_Type", "Tumor_Seq_Allele2","Tumor_Sample_Barcode",'AD']]

mafzao.to_csv('mafzao.maf', sep='\t', index=False)
mafzao_sem_variant.to_csv('mafzao_sem_variant.maf', sep=';', index=False)
mafzao_sem_variant_dbscan.to_csv('mafzao_sem_variant_dbscan.maf', sep=';', index=False)
#mafzao_sem_variant_dbscan_Hugo_com_ponto.to_csv('mafzao_sem_variant_dbscan_Hugo_com_ponto.maf', sep='\t', index=False)
#mafzao_limpo.to_csv('mafzao_limpo.maf', sep='\t', index=False)

#print(mafzao_sem_variant_dbscan)

print(f'Mafzao --> {len(mafzao)}')
print(f'Mafzao Variant--> {len(mafzao_sem_variant)}')
print(f'Mafzao Variant e DbScan --> {len(mafzao_sem_variant_dbscan)}')
print(f'Mafzao Variant e DbScan e Hugo sem ponto --> {len(mafzao_sem_variant_dbscan_Hugo_com_ponto)}')

#maf_pronto_comut = mafzao_limpo[['AD', 'Hugo_Symbol','Variant_Classification']]
#maf_pronto_comut.rename(columns={'AD':'sample', 'Hugo_Symbol': 'category', 'Variant_Classification': 'value'}, inplace=True)
#print(maf_pronto_comut)
#toy_comut = comut.CoMut()
#data = pd.read_csv('tutorial_data/tutorial_mutation_data.tsv', sep = '\t')


# Caminho para o arquivo JSON
arquivo_json = "Scripts/Json/BPP237-N_BPP237-T.hard-filtered.vcf.annotated"


fontes = ['clingen', 'clinvar', 'variants', 'gnomAD-preview', 'clingenDosageSensitivityMap' ]
#clingen_data = []

# Decodificar o conteúdo do arquivo JSON diretamente
with open(f'{arquivo_json}.json', 'r') as file:
    json_data = json.load(file)
    
# Escrever o arquivo JSON com identação
with open(f'{arquivo_json}_identado.json', 'w') as file:
    json.dump(json_data, file, indent=4)
    
for fonte in fontes:
    print(fonte)
    clingen_data = []
    for position in json_data["positions"]:
        if fonte in position:
            for entry in position[fonte]:
                clingen_data.append(entry)
                
    df_clingen = pd.DataFrame(clingen_data)
    #print(df_clingen)
    df_clingen.to_csv(f'{fonte}.csv', sep=';')
    #print(df_clingen)


arquivo_json = "Scripts/Json/BPP237-N_BPP237-T.hard-filtered.vcf.annotated.json"
with open(arquivo_json, 'r') as arquivo:
    # Carregar os dados do JSON para uma variável em Python
    dados = json.load(arquivo)

#data = json.loads(arquivo_json)

filtered_data = []

for posicao in dados['positions']:
    for variante in posicao['variants']:
        if variante['variantType'] in ['SNV','insertion']:
            filtered_data.append(variante)
            
df = pd.DataFrame(filtered_data)
df.to_csv('Scripts/nirvana_hard_filteres.csv', sep=';')
#print(df)



arquivo_json = "Scripts/Json/BPP237-N_BPP237-T.hard-filtered.vcf.annotated.json"
with open(arquivo_json, 'r') as arquivo:
    # Carregar os dados do JSON para uma variável em Python
    dados = json.load(arquivo)
    
clinvar_data =[]
# Iterar através de cada variante e extrair informações de clinvar
for posicao in dados['positions']:
    for variant in posicao['variants']:
        for clinvar_entry in variant.get('clinvar', []):
            clinvar_data.append(clinvar_entry)

# Criar um DataFrame com os dados de clinvar
clinvar_df = pd.DataFrame(clinvar_data)

# Mostrar o DataFrame
print(clinvar_df)
clinvar_df.to_csv('Scripts/clinvar_nirvana.csv', sep=';')
"""

#Juntando todos os JSONs do nirvana
diretorio = '/home/mirian/Documentos/Mirian/CPOM/Leticia_Pulmao/Scripts/nirvana/'

# Lista para armazenar os nomes dos arquivos
lista_arquivos = []

# Percorre todos os arquivos no diretório
for nome_arquivo in os.listdir(diretorio):
    caminho_completo = os.path.join(diretorio, nome_arquivo)
    lista_arquivos.append(caminho_completo)

#print(lista_arquivos)
# Agora 'lista_arquivos' contém todos os nomes dos arquivos no diretório especificado
df_nirvana = pd.DataFrame()
for arquivo in lista_arquivos:
    with open(arquivo, 'r') as arquivo:
    # Carregar os dados do JSON para uma variável em Python
        dados = json.load(arquivo)
    
    padrao = r"([^/]+)-N"
    # Buscar pela correspondência
    resultado = re.search(padrao, str(arquivo))
    
    if resultado:
        id_arquivo = resultado.group(1)
    else:
        id_arquivo = "Desconhecido"
    #print(id_arquivo)
    
    clinvar_data =[]
    # Iterar através de cada variante e extrair informações de clinvar
    for posicao in dados['positions']:
        for variant in posicao['variants']:
            for clinvar_entry in variant.get('clinvar', []):
                clinvar_data.append(clinvar_entry)
        
    # Criar um DataFrame com os dados de clinvar
    clinvar_df = pd.DataFrame(clinvar_data)
    clinvar_df['ID'] = id_arquivo
    #print(clinvar_df)
    df_nirvana = pd.concat([df_nirvana, clinvar_df], ignore_index=True)
    #df_nirvana['ID']=resultado
 
relacao_id_ad = pd.read_csv('Scripts/Id_Ad.csv', sep=';')
relacao_id_ad['ID']=relacao_id_ad['Codigo'].apply(lambda x: x.split('-')[0])
relacao_id_ad.drop(columns='Codigo', inplace=True)

df_nirvana = df_nirvana.merge(relacao_id_ad, how = 'left', on='ID')

#print(df_nirvana)  
#df_nirvana.to_csv('Scripts/nirvanao.csv', sep = ';')
df_nirvana['significance2']=df_nirvana['significance'].apply(lambda x: str(x[0]))
#print(df_nirvana)
#print(df_nirvana['significance2'].unique())

df_nirvana_filtrado = df_nirvana[~df_nirvana['significance2'].isin(['benign', 'likely benign'])]
#print(df_nirvana_filtrado['significance2'].unique())

df_nirvana_filtrado.to_csv('Scripts/nirvanao_filtrado.csv', sep = ';')
