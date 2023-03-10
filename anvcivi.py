import pandas as pd
import yaml
import os
import json
import argparse
from flatten_json import flatten


def A_ref_chose(str_aac, gene_transcript_dict):
    rt_a = ''
    for aa in str_aac.split(','):
        alist = aa.split(':')
        if alist[0] in gene_transcript_dict.keys() and alist[1] == gene_transcript_dict[alist[0]]:
            rt_a = aa
            break
    if rt_a == '':
        rt_a = str_aac.split(',')[0]
        waringtxt = '警告：输入基因{}与转录本{}与MANE转录本{}冲突'.format( rt_a.split(':')[0],rt_a.split(':')[1],gene_transcript_dict.get(rt_a.split(':')[0])) 
        return rt_a , waringtxt
    else:

        return rt_a , None

def connect_medic(row, pd_var_medic):

    alist = row['AAChange_MANE'].split(':')

    tmppd = pd_var_medic[pd_var_medic['gene'] == alist[0]] 

    #特殊情况之annovar不给p.var，或者格式不对劲
    if len(alist) == 4  and alist[3].upper() in tmppd['variant'].values:
        return   alist[3].upper() , tmppd[tmppd['variant'] == cval][[
            'clinical_significance', 'disease', 'drugs', 'evidence_level','citation_id',
            'evidence_statement'
        ]].to_dict('records')
    elif len(alist) != 5 :
        return row['AAChange_MANE'], None


    pval = alist[4].replace('p.', '').replace('*', 'X').upper()  #氨基酸变化
    cval = alist[3].upper()   #cDNA变化

    if pval in tmppd['variant'].values:  #蛋白质变化完全匹配,返回row['AAChange_MANE']

        return row['AAChange_MANE'], tmppd[tmppd['variant'] == pval][[
            'clinical_significance', 'disease', 'drugs', 'evidence_level','citation_id','evidence_statement']].to_dict('records')

    elif pval[:-1] in tmppd['variant'].values:  #判断仅关注蛋白质位置而非完整变异
        return   pval[:-1] , tmppd[tmppd['variant'] == pval[:-1]][[
            'clinical_significance', 'disease', 'drugs', 'evidence_level','citation_id',
            'evidence_statement'
        ]].to_dict('records')

    elif cval in tmppd['variant'].values :   # 判断 cDNA变化
        return   pval[:-1] , tmppd[tmppd['variant'] == cval][[
            'clinical_significance', 'disease', 'drugs', 'evidence_level','citation_id',
            'evidence_statement'
        ]].to_dict('records')

    elif 'MUTATION' in tmppd['variant'].values:
        return 'MUTATION', tmppd[tmppd['variant'] == 'MUTATION'][['clinical_significance', 'disease', 'drugs','evidence_level', 'citation_id','evidence_statement']].to_dict('records')
    
    else:
        return row['AAChange_MANE'], None



def var_civic(annovar_variant,civic_data,gene_transcript_kv ):
    '''
    输入annovar变体列表 如 AAChange.refGene 
    '''
    #print(annovar_variant['AAChange.refGene'].apply(A_ref_chose,result_type="expand" ,args=(gene_transcript_kv,)))
    annovar_variant['AAChange_MANE'],annovar_variant['select_waring'] = zip(*annovar_variant['AAChange.refGene'].apply(A_ref_chose, args=(gene_transcript_kv,) ))

    annovar_variant[['AA_match', 'INFO']] = annovar_variant.apply(connect_medic,
                                                      args=(civic_data, ),
                                                      axis=1,
                                                      result_type='expand')

    # https://datascientyst.com/normalize-json-dict-new-columns-pandas/   
    # json 再explot
    table_var = annovar_variant[~(annovar_variant['INFO'].isna())][['AAChange_MANE', 'INFO']]
    table_var = table_var.explode('INFO')
    table_var = pd.concat([table_var, table_var['INFO'].apply(flatten).apply(pd.Series)], axis=1)

    waring_pd= annovar_variant[~(annovar_variant['select_waring'].isna()) ][['select_waring']].drop_duplicates().reset_index(drop=True)
    annovar_variant['isannotate']=(~annovar_variant['INFO'].isna())


    var_af = annovar_variant[['AAChange_MANE','isannotate']]#.drop_duplicates().reset_index(drop=True)


    table_var_af_medic = table_var[[
    'AAChange_MANE', 'citation_id', 'clinical_significance', 'disease',
    'drugs','evidence_level'
    ]]
    table_var_af_medic = table_var_af_medic.groupby(
    ['AAChange_MANE',
     'clinical_significance', 'disease', 'drugs','evidence_level'])['citation_id'].apply(
         lambda x: ','.join(set(x.values))).reset_index()
    
    #print(table_var)

    table_ev = table_var[['citation_id', 'evidence_statement']].reset_index(drop=True)

    return var_af, table_var_af_medic , table_ev ,waring_pd

#https://civicdb.org/downloads/01-Dec-2022/01-Dec-2022-ClinicalEvidenceSummaries.tsv

def main():
    
    script_path = os.path.split(os.path.realpath(__file__))[0]

    parser = argparse.ArgumentParser(description='anvcivi argparse')
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        required=True,
                        help='输入文件，annovar输出调整得到')
    parser.add_argument('-d',
                        '--disease',
                        type=str,
                        required=True,
                        help='选择只注释的疾病,需要civic本身支持，ALL注释全疾病')
    parser.add_argument('-c',
                        '--config',
                        type=str,
                        default=script_path + '/anvcivi.config.yml',
                        help='config路径，默认位于 $scritpt_path/anvcivi.config.yml')
    parser.add_argument('-o',
                        '--out_dir',
                        type=str,
                        required=True,
                        help='输出文件，tsv格式文本')
    args = parser.parse_args()

    with open(args.config, 'r', encoding='utf8') as conf:
        path_dic = yaml.safe_load(conf)


    
    civic_data = pd.read_csv(path_dic['civic'], sep='\t')
    with open(path_dic['mane'], 'r') as fp:
        gene_transcript_dict = json.load(fp)

    annovar_variant = pd.read_csv(args.input, sep='\t')
    

    with open(path_dic['sig_term'], 'r', encoding='utf8') as fp:
        sig_term_dict = json.load(fp)

    civic_data = civic_data.replace({'clinical_significance':sig_term_dict})
    #civic_data=civic_data.replace({"disease": disease_dict})
    civic_data['variant']=civic_data['variant'].str.replace('*', 'X',regex=False)



    if  args.disease != 'ALL':
        civic_data=civic_data[civic_data['disease']==args.disease]

    if not os.path.exists(args.out_dir):  #判断所在目录下有该文件名的文件夹
        os.makedirs(args.out_dir)


    table_var_af ,table_var_af_medic,table_ev ,waring_pd=  var_civic(annovar_variant,civic_data,gene_transcript_dict)



    table_var_af.to_csv(args.out_dir + '/var_af.tsv', sep='\t', index=None)
    table_var_af_medic.to_csv(args.out_dir + '/var_medic.tsv', sep='\t', index=None)
    table_ev.to_csv(args.out_dir + '/evidence.tsv', sep='\t', index=None)
    print(waring_pd)


    
if __name__ == '__main__':
    
    main()

