
import pandas as pd
import json
import os
def MANE_trans(MANE_zip_file):
    
    gene_transcript = pd.read_csv(MANE_zip_file, sep='\t', compression='gzip')


    gene_transcript = gene_transcript[['symbol', 'RefSeq_nuc']].copy()
    gene_transcript['RefSeq_nuc'] = gene_transcript['RefSeq_nuc'].str.replace('.[0-9]$', '',regex=True)
    gene_transcript_kv = gene_transcript.set_index('symbol').to_dict(orient='dict')['RefSeq_nuc']

    return gene_transcript_kv

if  __name__ == '__main__':
    gtdict = MANE_trans('../source/MANE.GRCh38.v1.0.summary.txt.gz')
    if not os.path.exists('../json'):  #判断所在目录下有该文件名的文件夹
        os.makedirs('../json')

    with open('../json/MANE_gene_transcript.json', 'w') as fp:
        json.dump(gtdict, fp)