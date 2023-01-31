import pandas as pd
import re
import os

def var_add_col(df, colname, pattern):
    df[colname] = df['variant'].str.contains(pattern, regex=True)
    return df
def list_change(list):
    rt_list=[]
    start=re.findall('^[A-Z]*[0-9]*', list[0])[0]
    if len(list) >1:
        for x in list:    
            if x !=start and len(x) ==1:  # D835H/Y
                rt_list.append(start+x)
            elif x =='K and Amplification'  :    #  V600E/K and Amplification
                rt_list.append(start+x)
            elif x =='V600E':
                rt_list.append('V600E and Amplification')
            else:                           # G12/G13
                rt_list.append(x)
    else:
        return list

    return rt_list
def str_cut(string):
    pt1 = re.compile(r"[(](C{1}.*?)[)]", re.S)   #最小匹配 匹配括号内 c.*(大写C)
    
    str = re.findall(pt1, string)

    
    if len(str) != 0:
        print(str)
        return str[0]
    else:
        return string

res = pd.read_csv('source/01-Dec-2022-ClinicalEvidenceSummaries.tsv', sep='\t')

res=res[res['evidence_type']=='Predictive'] # 用药预测
res['citation_id']=res['citation_id'].astype(str)+'('+res['source_type']+')'


res['variant'] = res['variant'].str.upper()  #全转大写
res['drugs'] = res['drugs'].str.split(',')
res=res.explode('drugs')
res['variant'] = res['variant'].str.split('/')
res['variant'] = res['variant'].apply(list_change)
res=res.explode('variant')
res['variant'] = res['variant'].apply(str_cut)

res = var_add_col(res, 'isppos', '^[A-Z]*[0-9]*[0-9]$')
res = var_add_col(res, 'ispvar', '^[A-Z][0-9]+[A-Z*]$|^[A-Z][0-9]+_{1}|^[A-Z][0-9]+FS|^[A-Z][0-9]+DEL')  # smp|T123_|T123FS
res = var_add_col(res, 'isexonvar', 'EXON\s{1}[0-9]+\s{1}[A-Z]+')
res = var_add_col(res, 'iscvar', 'C\.')
res = var_add_col(res, 'isfusion', ':{2}[A-Z0-9]+[A-Z0-9]$')
res = var_add_col(res, 'ismutation', '^MUTATION$')
#
res['匹配次数']=res[['isppos', 'ispvar','isexonvar','iscvar','isfusion','ismutation']].sum(axis=1)

not_single_pattern_var=res[['gene','variant','isexonvar','isppos', 'ispvar','iscvar','isfusion','ismutation','匹配次数']].query('匹配次数>1')
if not_single_pattern_var.shape[0] !=0:
    print(not_single_pattern_var)

 #有关分子谱对 治疗 反应的影响的证据,即用药相关部分
#res=res[res['evidence_type']=='Oncogenic'] 

if not os.path.exists('civic_clean/'):  #判断所在目录下有该文件名的文件夹
    os.makedirs('civic_clean/')
#------------------------------------#
no_clas=res.query('(isexonvar == False) & (isppos == False) & (ispvar == False) & (iscvar == False) & (isfusion == False) & (ismutation == False)' )
no_clas=no_clas.reset_index(drop=True)
no_clas=no_clas[['gene','variant','evidence_level','drugs','clinical_significance','evidence_statement']].sort_values(by='evidence_level', ascending=True)
no_clas.to_csv('civic_clean/unaccept_var.tsv', sep='\t')
#------------------------#
# acclas[['gene','variant','isexonvar','isppos', 'ispvar','iscvar','isfusion']].to_csv('aceptvar.tsv', sep='\t')

acclas=res.query('~( (isexonvar == False) & (isppos == False) & (ispvar == False) & (iscvar == False) & (isfusion == False) & (ismutation == False) )' )
acclas=acclas.reset_index(drop=True)
acclas.to_csv('civic_clean/accept_var_info.tsv', sep='\t')

#---------------------------------------#

#print(acclas['drugs'].value_counts())



'''
import sweetviz as sv
# read dataset

my_report = sv.analyze(acclas)
my_report.show_html()

'''
