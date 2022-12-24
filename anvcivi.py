import pandas as pd

#https://civicdb.org/downloads/01-Dec-2022/01-Dec-2022-ClinicalEvidenceSummaries.tsv
res=pd.read_csv('01-Dec-2022-ClinicalEvidenceSummaries.tsv',sep='\t')
res=res[['gene','variant','disease','drugs','drug_interaction_type','evidence_statement','citation_id','source_type','variant_summary']]
res= res[~((res['drugs'].isna()) ) ]  #只管有药  (res['disease'].isna()) & 
res['variant']=res['variant'].str.split('/')
res=res.explode('variant')
res['variant']=res['variant'].str.replace('*','',regex=False)

#res.to_csv('var_medic.tsv',sep='\t')

#MANE转录本
#https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.summary.txt.gz

gene_tcr=pd.read_csv('MANE.GRCh38.v1.0.summary.txt.gz',sep='\t',compression='gzip')
gene_tcr=gene_tcr[['symbol','RefSeq_nuc']]
gene_tcr['RefSeq_nuc'] = gene_tcr['RefSeq_nuc'].str.replace('.[0-9]$','',regex=True)
gene_tcr_kv=gene_tcr.set_index('symbol').to_dict(orient='dict')['RefSeq_nuc']

#print(gene_tcr)
#gene_tcr.to_csv('gene_nm.tsv',sep='\t',index=None)


##annovar hg19 results
annovar_out=pd.read_csv('anovar_add_af.tsv',sep='\t')
#print(annovar_out)

def A_ref_chose(str_aac,gene_tcr_dic):
    rt_a = ''
    for aa in str_aac.split(','):
        alist=aa.split(':')
        if alist[0] in gene_tcr_dic.keys() and alist[1] == gene_tcr_dic[alist[0]]:
            rt_a=aa
            break
    if rt_a == '' :
        rt_a = str_aac.split(',')[0]

    return rt_a

annovar_out['AAChange_MANE']=annovar_out['AAChange.refGene'].apply(A_ref_chose,args=(gene_tcr_kv,))

def connect_medic(row,pd_var_medic):

    alist=row['AAChange_MANE'].split(':')
    extyp=row['ExonicFunc.refGene']
    tmppd=pd_var_medic[pd_var_medic['gene']==alist[0]]
    pval=alist[-1].replace('p.', '')  #有的无氨基酸变化，利用这点可以比对上c.***

    if  pval in  tmppd['variant'].values: #氨基酸变化完全匹配,返回row['AAChange_MANE']；如果是indel,则可以返回修改后的
        
        return row['AAChange_MANE'],tmppd[tmppd['variant']==pval][['disease','drugs','evidence_statement']].to_dict('records')
        
    else:
        return row['AAChange_MANE'],None
    
annovar_out[['AA_match','INFO']]=annovar_out.apply(connect_medic,args=(res,),axis=1,result_type='expand')
#annovar_out=annovar_out.explode('INFO')

table_var=annovar_out[~(annovar_out['INFO'].isna())][['AA_match','af','INFO']]
table_var=table_var.explode('INFO')

#x=table_var['INFO'].apply(pd.Series)
# https://datascientyst.com/normalize-json-dict-new-columns-pandas/
from flatten_json import flatten
table_var= pd.concat([table_var, table_var['INFO'].apply(flatten).apply(pd.Series)], axis=1)

table_var['af']

table_var_af=table_var[['AA_match','af']].drop_duplicates().reset_index(drop=True)

table_var_af.to_markdown('a.md')

table_var_af_md=table_var[['AA_match','af','disease','drugs']].reset_index(drop=True)

table_var_af_md.to_markdown('x.md')

table_var_ev=table_var[['AA_match','disease','drugs','evidence_statement']].reset_index(drop=True)

table_var_ev.to_markdown('y.md')
