
import os
import gzip
import shutil

#Descompactando e colocando na pasta nirvana

# Define o diretório raiz onde as pastas estão localizadas
diretorio_raiz = '/run/user/1002/gvfs/smb-share:server=nas-cpom.local,share=dados/Leticia/sophia/fastq_sophia_bioinfogpt/'
# Define o diretório onde você deseja extrair os arquivos
diretorio_destino = '/home/mirian/Documentos/Mirian/CPOM/Leticia_Pulmao/Scripts/nirvana/'

# Caminhar por todos os subdiretórios e arquivos no diretório raiz
for raiz, pastas, arquivos in os.walk(diretorio_raiz):
    for nome_arquivo in arquivos:
        # Verificar se o arquivo termina com .hard-filtered.vcf.annotated.json.tar.gz
        if nome_arquivo.endswith('.hard-filtered.vcf.annotated.json.gz'):
            caminho_arquivo = os.path.join(raiz, nome_arquivo)

            # Definir o nome do arquivo descompactado (removendo a extensão .gz)
            nome_arquivo_descompactado = nome_arquivo[:-3]
            caminho_arquivo_descompactado = os.path.join(diretorio_destino, nome_arquivo_descompactado)

            # Descompactar o arquivo
            with gzip.open(caminho_arquivo, 'rb') as f_in:
                with open(caminho_arquivo_descompactado, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            print(f"Arquivo {nome_arquivo} descompactado para {caminho_arquivo_descompactado}")
