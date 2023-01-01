import pandas as pd
import yaml
import os
import argparse

script_path=os.path.split(os.path.realpath(__file__))[0]



parser = argparse.ArgumentParser(description='anvcivi argparse')
parser.add_argument('-i','--input', type=str, required=True,help='输入文件，annovar输出调整得到')
parser.add_argument('-c','--config', type=str, default=script_path+'/anvcivi.config.yml',help='config路径，默认位于 $scritpt_path/anvcivi.config.yml')
parser.add_argument('-o','--out_dir',type=str,required=True,help='输出文件，tsv格式文本' )
args = parser.parse_args()



with open(args.config,'r',encoding='utf8') as conf:
    path_dic = yaml.safe_load(conf)


#https://civicdb.org/downloads/01-Dec-2022/01-Dec-2022-ClinicalEvidenceSummaries.tsv
res=pd.read_csv(path_dic['civic'],sep='\t')
#print(set(res['clinical_significance'].values.tolist()))
res=res[['gene','variant','clinical_significance','disease','drugs','drug_interaction_type','evidence_statement','citation_id','source_type','variant_summary']]
res= res[~((res['drugs'].isna()) ) ]  #只管有药  (res['disease'].isna()) & 
res['variant']=res['variant'].str.split('/')
res=res.explode('variant')

res.to_csv('tmp.csv')
#初步判断，civic是根据hgvs标准命名，但是annovar的氨基酸水平命名并不遵守此规则。
# 初步设想是:先根据annovar  -hgvs 获得满足hgvs命名规则的cDNA突变命名，再借由python hgvs包，解析出氨基酸水平突变命名。最终将其与用药信息关联。


#MANE转录本
#https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.summary.txt.gz

gene_tcr=pd.read_csv(path_dic['mane'],sep='\t',compression='gzip')
gene_tcr=gene_tcr[['symbol','RefSeq_nuc']]
gene_tcr['RefSeq_nuc'] = gene_tcr['RefSeq_nuc'].str.replace('.[0-9]$','',regex=True)
gene_tcr_kv=gene_tcr.set_index('symbol').to_dict(orient='dict')['RefSeq_nuc']

#print(gene_tcr)
#gene_tcr.to_csv('gene_nm.tsv',sep='\t',index=None)


##annovar hg19 results
annovar_out=pd.read_csv(args.input,sep='\t')
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

    if  pval in  tmppd['variant'].values: #氨基酸变化完全匹配,返回row['AAChange_MANE']
        
        return row['AAChange_MANE'],tmppd[tmppd['variant']==pval][['clinical_significance','disease','drugs','evidence_statement']].to_dict('records')
    
    elif pval[:-1] in  tmppd['variant'].values: #判断仅关注氨基酸位置而非完整变异
        return row['AAChange_MANE'],tmppd[tmppd['variant']==pval[:-1]][['clinical_significance','disease','drugs','evidence_statement']].to_dict('records')
    
    #elif alist[3] in   #exon n#直接判断 是否包含exon indel字符较危险，最好还是将数据库文件分列，彻底清洗干净，搁置
    #氨基酸变化不完全匹配，例如:如果是indel,则可以返回EGFR exon19 del之类的自定义字符串


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


table_var_af=table_var[['AA_match','af']].drop_duplicates().reset_index(drop=True)

table_var_af.to_csv(args.out_dir+'/var_af.tsv',sep='\t')

table_var_af_md=table_var[['AA_match','af','clinical_significance','disease','drugs']].reset_index(drop=True)

table_var_af_md.to_csv(args.out_dir+'/var_medic.tsv',sep='\t')

table_var_ev=table_var[['AA_match','disease','drugs','evidence_statement']].reset_index(drop=True)

table_var_ev.to_csv(args.out_dir+'/var_evidence.tsv',sep='\t')
