from anvcivi import var_civic
import yaml
import pandas as pd
import json
from pywebio.input import textarea ,radio
from pywebio.output import put_html,put_text #,put_row,put_table,span
from pywebio import start_server

def web_anvcivi():

    with open('web_anvcivi.yaml' ,'r', encoding='utf8') as conf:
        path_dic = yaml.safe_load(conf)

    civic_data = pd.read_csv(path_dic['civic'], sep='\t')
    with open(path_dic['mane'], 'r') as fp:
        gene_transcript_dict = json.load(fp)

    with open(path_dic['sig_term'], 'r', encoding='utf8') as fp:
        sig_term_dict = json.load(fp)

    civic_data = civic_data.replace({'clinical_significance':sig_term_dict})
    #civic_data=civic_data.replace({"disease": disease_dict})
    civic_data['variant']=civic_data['variant'].str.replace('*', 'X',regex=False)


    input_txt = textarea('annovar格式变体',value='''
EGFR:NM_001346898:exon21:c.T2573G:p.L858R,EGFR:NM_001346900:exon21:c.T2414G:p.L805R,EGFR:NM_005228:exon21:c.T2573G:p.L858R
EGFR:NM_001346941:exon18:c.2154_2155delinsAA:p.G719S
KRAS:NM_033360:exon2:c.33_34delinsAT:p.G12C''')


    disease_list=list(set(civic_data['disease'].values)) #.update('all')
    disease_list.append('ALL')
    #print(disease_list)
    disease_chose=radio(options=disease_list,required=True,value='ALL')

    if disease_chose != 'ALL':
        civic_data=civic_data[civic_data['disease']==disease_chose]
    
    annovar_variant = pd.DataFrame({ 'AAChange.refGene' : input_txt.split("\n")})
    
    try :
        table_var_af ,table_var_af_medic,table_ev,warning =  var_civic(annovar_variant,civic_data,gene_transcript_dict)
        put_html(table_var_af.to_html(border=0)) 
        put_html(table_var_af_medic.to_html(border=0))
        put_html(table_ev.to_html(border=0))
        #print(table_var_af ,table_var_af_medic,table_ev,warning)
        put_html(warning.to_html(border=0))
    except:
        put_text('意外错误，可能是输入突变与选择疾病无关联用药，可尝试选择ALL')

    

    
#annovar_variant = pd.read_csv(args.input, sep='\t')
    
if __name__ == '__main__':
    start_server(web_anvcivi,port=8800)